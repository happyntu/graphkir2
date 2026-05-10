"""Re-run the legacy typing stage from an existing benchmark preset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    """Build the rerun typing parser."""
    parser = argparse.ArgumentParser(
        description="Re-run Graph-KIR typing from an existing benchmark preset.",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to a benchmark preset JSON file.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path to write the merged allele TSV.",
    )
    parser.add_argument(
        "--multi-map-mode",
        default="discard",
        choices=["discard", "keep-all", "best-only", "weighted", "margin", "likelihood"],
        help="How to handle multi-mapped reads during the re-run typing stage.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=600,
        help="Typing top-n parameter passed to the legacy model.",
    )
    parser.add_argument(
        "--exon-weight",
        type=float,
        help="Extra likelihood weight applied to exon variants.",
    )
    parser.add_argument(
        "--gene-exon-weights",
        default="",
        help="Optional `GENE:WEIGHT,...` overrides for exon weighting.",
    )
    parser.add_argument(
        "--confusion-gene-groups",
        default="",
        help="Optional `GENE1/GENE2,...` groups for confusion-aware best-only scoring.",
    )
    parser.add_argument(
        "--margin-gap",
        type=float,
        help="Top1-top2 score gap threshold below which near-tie alignments share weight.",
    )
    parser.add_argument(
        "--margin-scale",
        type=float,
        help="Softmax score scale used when `--multi-map-mode margin` keeps near-tie alignments.",
    )
    parser.add_argument(
        "--ambiguity-neutral-prob",
        type=float,
        help="Neutral allele probability used for non-target read mass in likelihood mode.",
    )
    parser.add_argument(
        "--ambiguity-target-weight-power",
        type=float,
        help="Power applied to ambiguous target alignment weights in likelihood mode.",
    )
    parser.add_argument(
        "--select-min-fraction-ratio",
        type=float,
        help="Minimum allele fraction as a ratio of the expected 1/CN fraction.",
    )
    return parser


def main() -> None:
    """Re-run typing for every sample in a benchmark preset."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir.kir_cn import loadCN
    from graphkir.kir_typing import (
        parseConfusionGeneGroups,
        parseGeneWeightSpec,
        selectKirTypingModel,
    )
    from graphkir2.benchmark import load_benchmark_preset
    from graphkir2.core.pipeline import GraphKir2Pipeline
    from graphkir2.io.manifest import load_sample_manifest

    args = build_parser().parse_args()
    preset = load_benchmark_preset(args.config)
    pipeline = GraphKir2Pipeline(preset.to_run_config())
    typing_plan = pipeline.build_typing_plan()
    manifest = load_sample_manifest(preset.input_csv)
    cn_hints = {sample.sample_id: sample.copy_number_hint for sample in manifest.samples}
    exon_weight = args.exon_weight if args.exon_weight is not None else preset.allele_exon_weight
    gene_weight_spec = args.gene_exon_weights or getattr(preset, "allele_gene_exon_weights", "")
    gene_exon_weights = parseGeneWeightSpec(gene_weight_spec)
    confusion_gene_groups = parseConfusionGeneGroups(
        args.confusion_gene_groups or getattr(preset, "allele_confusion_groups", "")
    )
    margin_gap = (
        args.margin_gap
        if args.margin_gap is not None
        else getattr(preset, "allele_margin_gap", 5.0)
    )
    margin_scale = (
        args.margin_scale
        if args.margin_scale is not None
        else getattr(preset, "allele_margin_scale", 2.0)
    )
    ambiguity_neutral_prob = (
        args.ambiguity_neutral_prob
        if args.ambiguity_neutral_prob is not None
        else getattr(preset, "allele_ambiguity_neutral_prob", 0.999)
    )
    ambiguity_target_weight_power = (
        args.ambiguity_target_weight_power
        if args.ambiguity_target_weight_power is not None
        else getattr(preset, "allele_ambiguity_target_weight_power", 1.0)
    )
    select_min_fraction_ratio = (
        args.select_min_fraction_ratio
        if args.select_min_fraction_ratio is not None
        else getattr(preset, "allele_select_min_fraction_ratio", 0.5)
    )

    rows: list[dict[str, str]] = []
    for sample in typing_plan.samples:
        model = selectKirTypingModel(
            typing_plan.strategy,
            sample.variant_json,
            top_n=args.top_n,
            variant_correction=True,
            multi_map_mode=args.multi_map_mode,
            exon_weight=exon_weight,
            gene_exon_weights=gene_exon_weights,
            confusion_gene_groups=confusion_gene_groups,
            margin_score_gap=margin_gap,
            margin_score_scale=margin_scale,
            ambiguity_neutral_prob=ambiguity_neutral_prob,
            ambiguity_target_weight_power=ambiguity_target_weight_power,
            select_min_fraction_ratio=select_min_fraction_ratio,
        )
        cn_path = cn_hints.get(sample.sample_id, "") or sample.cn_tsv
        cn = loadCN(cn_path)
        called_alleles, warning_genes = model.typing(cn)
        rows.append(
            {
                "name": sample.output_prefix,
                "alleles": "_".join(called_alleles),
                "warnings": "_".join(warning_genes),
            }
        )

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    print(f"Wrote typing results to {output_path}")


if __name__ == "__main__":
    main()
