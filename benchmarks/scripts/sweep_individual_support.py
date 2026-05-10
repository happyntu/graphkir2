"""Sweep individual-support filters for allele candidate selection."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PANELS: tuple[tuple[str, str, str], ...] = (
    (
        "synthetic-functional8",
        "benchmarks/configs/synthetic-functional8-likelihood-exon2-minfrac07.json",
        "benchmarks/configs/synthetic-functional8.json",
    ),
    (
        "synthetic-functional8x6",
        "benchmarks/configs/synthetic-functional8x6-likelihood-exon2-minfrac07.json",
        "benchmarks/configs/synthetic-functional8x6.json",
    ),
    (
        "synthetic-difficult5",
        "benchmarks/configs/synthetic-difficult5-likelihood-exon2-minfrac07.json",
        "benchmarks/configs/synthetic-difficult5.json",
    ),
)


def build_parser() -> argparse.ArgumentParser:
    """Build the sweep parser."""
    parser = argparse.ArgumentParser(
        description="Sweep individual-support candidate filters.",
    )
    parser.add_argument(
        "--thresholds",
        default="0.02,0.05,0.1,0.2,0.4",
        help="Comma-separated min individual support ratios.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=600,
        help="Typing top-n parameter.",
    )
    return parser


def select_with_ratio(
    result: object,
    min_fraction_ratio: float,
    min_individual_ratio: float,
) -> list[str]:
    """Select the first candidate passing fraction and individual-support filters."""
    candidate_ids = []
    cn = result.n  # type: ignore[attr-defined]
    expected_fraction = 1 / cn
    for index in range(len(result.fraction)):  # type: ignore[attr-defined]
        fraction = result.fraction[index]  # type: ignore[attr-defined]
        if not all(value >= expected_fraction * min_fraction_ratio for value in fraction):
            continue
        individual_values = np.abs(result.value_sum_indv[index])  # type: ignore[attr-defined]
        if (
            len(individual_values)
            and individual_values.max() > 0
            and individual_values.min() / individual_values.max() < min_individual_ratio
        ):
            continue
        candidate_ids.append(index)
    best_id = (candidate_ids or [0])[0]
    return result.allele_name[best_id]  # type: ignore[attr-defined]


def main() -> None:
    """Run the individual-support sweep."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir.kir_cn import loadCN
    from graphkir.kir_typing import TypingWithPosNegAllele, selectKirTypingModel
    from graphkir2.benchmark import load_benchmark_preset
    from graphkir2.core.pipeline import GraphKir2Pipeline
    from graphkir2.io.manifest import load_sample_manifest

    args = build_parser().parse_args()
    thresholds = [float(value) for value in args.thresholds.split(",") if value.strip()]
    result_rows: list[dict[str, str]] = []

    for threshold in thresholds:
        for label, preset_config, eval_config in PANELS:
            preset = load_benchmark_preset(preset_config)
            plan = GraphKir2Pipeline(preset.to_run_config()).build_typing_plan()
            manifest = load_sample_manifest(preset.input_csv)
            cn_hints = {
                sample.sample_id: sample.copy_number_hint
                for sample in manifest.samples
            }
            rows: list[dict[str, str]] = []
            for sample in plan.samples:
                model = selectKirTypingModel(
                    plan.strategy,
                    sample.variant_json,
                    top_n=args.top_n,
                    variant_correction=True,
                    multi_map_mode="likelihood",
                    exon_weight=2.0,
                    ambiguity_neutral_prob=getattr(
                        preset,
                        "allele_ambiguity_neutral_prob",
                        0.999,
                    ),
                    select_min_fraction_ratio=0.7,
                )
                if not isinstance(model, TypingWithPosNegAllele):
                    raise TypeError("individual-support sweep requires full allele typing")
                cn = model._normalizeGeneCopyNumbers(
                    loadCN(cn_hints.get(sample.sample_id, "") or sample.cn_tsv)
                )
                alleles: list[str] = []
                warnings: list[str] = []
                for gene, copy_number in cn.items():
                    if not copy_number:
                        continue
                    _, reads_num = model.typingPerGene(gene, copy_number)
                    result = model._result[gene][-1]
                    chosen = select_with_ratio(result, 0.7, threshold)
                    pure_gene = gene.split("*")[0]
                    alleles.extend(
                        allele if allele != "fail" else f"{pure_gene}*"
                        for allele in chosen
                    )
                    if reads_num < 100:
                        warnings.append(gene)
                rows.append(
                    {
                        "name": sample.output_prefix,
                        "alleles": "_".join(alleles),
                        "warnings": "_".join(warnings),
                    }
                )

            threshold_label = str(threshold).replace(".", "p")
            pred_tsv = Path(
                f"benchmarks/results/{label}/tmp_top{args.top_n}_ratio_{threshold_label}.allele.tsv"
            )
            pred_tsv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(rows).to_csv(pred_tsv, sep="\t", index=False)
            bundle_dir = Path(
                f"benchmarks/results/{label}-tmp-top{args.top_n}-ratio-{threshold_label}-bundle"
            )
            subprocess.run(
                [
                    sys.executable,
                    "benchmarks/scripts/run_compare.py",
                    "--config",
                    eval_config,
                    "--evaluate",
                    "--pred-tsv",
                    str(pred_tsv),
                    "--multi-map-policy-override",
                    "likelihood",
                    "--output-dir",
                    str(bundle_dir),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            summary = pd.read_csv(bundle_dir / "summary.tsv", sep="\t").iloc[0]
            result_rows.append(
                {
                    "threshold": str(threshold),
                    "panel": label,
                    "three_digit_f1": str(summary["three_digit_f1"]),
                    "five_digit_f1": str(summary["five_digit_f1"]),
                    "seven_digit_f1": str(summary["seven_digit_f1"]),
                    "pred_tsv": str(pred_tsv),
                    "bundle_dir": str(bundle_dir),
                }
            )

    output_path = Path("benchmarks/results/individual-support-sweep.tsv")
    pd.DataFrame(result_rows).to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {output_path}")
    print(pd.DataFrame(result_rows).to_string(index=False))


if __name__ == "__main__":
    main()
