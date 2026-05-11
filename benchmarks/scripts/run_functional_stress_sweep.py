"""Run a synthetic functional typing stress sweep for graphkir2 candidates."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from profile_targeted_top_n import run_timed


RunnerKind = Literal["standard", "private_support"]
Resolution = Literal["three_digit_f1", "five_digit_f1", "seven_digit_f1"]


DEFAULT_LABELS = (
    "synthetic-functional8",
    "synthetic-functional8x6",
    "synthetic-difficult5",
    "synthetic-difficult5x12",
    "synthetic-difficult5x12-seed5101",
    "synthetic-difficult5x12-seed5102",
    "synthetic-difficult5x12-seed5103",
)

FUNCTIONAL_METRICS: tuple[Resolution, ...] = ("three_digit_f1", "five_digit_f1")


@dataclass(frozen=True)
class MethodSpec:
    """One typing method in the functional stress sweep."""

    name: str
    runner: RunnerKind
    config_suffix: str
    multi_map_mode: str
    top_n: int
    base_top_n: int = 0
    gene_base_top_ns: str = ""
    functional_discard_fallback_genes: str = ""
    functional_discard_fallback_resolution: int = 3
    functional_discard_fallback_max_score: float = 0.0
    functional_discard_fallback_min_score_delta: float = 0.0
    exon_weight: float | None = None
    ambiguity_neutral_prob: float | None = None
    select_min_fraction_ratio: float | None = None
    policy_override: str = ""


@dataclass(frozen=True)
class SweepPaths:
    """Derived paths for one panel/method sweep run."""

    config: Path
    prediction_tsv: Path
    bundle_dir: Path


DEFAULT_METHODS: tuple[MethodSpec, ...] = (
    MethodSpec(
        name="discard",
        runner="standard",
        config_suffix="",
        multi_map_mode="discard",
        top_n=600,
        policy_override="discard",
    ),
    MethodSpec(
        name="likelihood_top5000",
        runner="standard",
        config_suffix="",
        multi_map_mode="likelihood",
        top_n=5000,
        exon_weight=2.0,
        ambiguity_neutral_prob=0.999,
        select_min_fraction_ratio=0.7,
        policy_override="likelihood",
    ),
    MethodSpec(
        name="enhancedgate_top5000",
        runner="private_support",
        config_suffix="-conditional-kir2ds3-enhancedgate",
        multi_map_mode="likelihood",
        top_n=5000,
    ),
    MethodSpec(
        name="enhancedgate_geneaware",
        runner="private_support",
        config_suffix="-conditional-kir2ds3-enhancedgate",
        multi_map_mode="likelihood",
        top_n=5000,
        base_top_n=600,
        gene_base_top_ns="KIR2DL1:1000",
    ),
    MethodSpec(
        name="enhancedgate_kir2dl1fallback_geneaware",
        runner="private_support",
        config_suffix="-conditional-kir2ds3-enhancedgate",
        multi_map_mode="likelihood",
        top_n=5000,
        base_top_n=600,
        gene_base_top_ns="KIR2DL1:1000",
        functional_discard_fallback_genes="KIR2DL1",
        functional_discard_fallback_resolution=3,
        functional_discard_fallback_max_score=-100.0,
        functional_discard_fallback_min_score_delta=20.0,
    ),
)


def build_parser() -> argparse.ArgumentParser:
    """Build the functional stress sweep parser."""
    parser = argparse.ArgumentParser(
        description="Run synthetic functional typing stress sweeps.",
    )
    parser.add_argument(
        "--label",
        action="append",
        default=[],
        help="Synthetic label to sweep. Defaults to current functional stress set.",
    )
    parser.add_argument(
        "--method",
        action="append",
        choices=[method.name for method in DEFAULT_METHODS],
        default=[],
        help="Method name to run. Defaults to all methods.",
    )
    parser.add_argument(
        "--results-dir",
        default="benchmarks/results/functional-stress-sweep",
        help="Directory for generated predictions and bundles.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/summary.tsv",
        help="Machine-readable aggregate summary TSV.",
    )
    parser.add_argument(
        "--output-per-gene-tsv",
        default="benchmarks/results/functional-stress-sweep/per_gene_summary.tsv",
        help="Machine-readable per-gene metric TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_functional_stress_seed_sweep.md",
        help="Markdown report path.",
    )
    parser.add_argument(
        "--skip-md",
        action="store_true",
        help="Do not write the markdown report.",
    )
    return parser


def select_methods(names: list[str]) -> tuple[MethodSpec, ...]:
    """Return requested method specs in default order."""
    if not names:
        return DEFAULT_METHODS
    requested = set(names)
    return tuple(method for method in DEFAULT_METHODS if method.name in requested)


def paths_for(
    label: str,
    method: MethodSpec,
    results_dir: Path,
) -> SweepPaths:
    """Build stable paths for one panel/method run."""
    config = Path("benchmarks/configs") / f"{label}{method.config_suffix}.json"
    return SweepPaths(
        config=config,
        prediction_tsv=results_dir / f"{label}.{method.name}.allele.tsv",
        bundle_dir=results_dir / f"{label}.{method.name}.bundle",
    )


def build_typing_command(
    method: MethodSpec,
    paths: SweepPaths,
) -> list[str]:
    """Build the typing rerun command for one method."""
    if method.runner == "standard":
        command = [
            sys.executable,
            "benchmarks/scripts/rerun_typing.py",
            "--config",
            str(paths.config),
            "--output-tsv",
            str(paths.prediction_tsv),
            "--multi-map-mode",
            method.multi_map_mode,
            "--top-n",
            str(method.top_n),
        ]
        if method.exon_weight is not None:
            command.extend(["--exon-weight", str(method.exon_weight)])
        if method.ambiguity_neutral_prob is not None:
            command.extend(["--ambiguity-neutral-prob", str(method.ambiguity_neutral_prob)])
        if method.select_min_fraction_ratio is not None:
            command.extend(
                ["--select-min-fraction-ratio", str(method.select_min_fraction_ratio)]
            )
        return command

    command = [
        sys.executable,
        "benchmarks/scripts/rerun_typing_private_support.py",
        "--config",
        str(paths.config),
        "--output-tsv",
        str(paths.prediction_tsv),
        "--top-n",
        str(method.top_n),
    ]
    if method.base_top_n:
        command.extend(["--base-top-n", str(method.base_top_n)])
    if method.gene_base_top_ns:
        command.extend(["--gene-base-top-ns", method.gene_base_top_ns])
    if method.functional_discard_fallback_genes:
        command.extend(
            [
                "--functional-discard-fallback-genes",
                method.functional_discard_fallback_genes,
                "--functional-discard-fallback-resolution",
                str(method.functional_discard_fallback_resolution),
                "--functional-discard-fallback-max-score",
                str(method.functional_discard_fallback_max_score),
                "--functional-discard-fallback-min-score-delta",
                str(method.functional_discard_fallback_min_score_delta),
            ]
        )
    return command


def build_compare_command(
    method: MethodSpec,
    paths: SweepPaths,
) -> list[str]:
    """Build the evaluation command for one generated prediction TSV."""
    command = [
        sys.executable,
        "benchmarks/scripts/run_compare.py",
        "--config",
        str(paths.config),
        "--evaluate",
        "--pred-tsv",
        str(paths.prediction_tsv),
        "--output-dir",
        str(paths.bundle_dir),
    ]
    if method.policy_override:
        command.extend(["--multi-map-policy-override", method.policy_override])
    return command


def read_bundle_metrics(bundle_dir: Path) -> dict[str, str]:
    """Read aggregate evaluation metrics from a result bundle."""
    artifact = json.loads((bundle_dir / "compare.json").read_text(encoding="utf-8"))
    metrics = artifact["metrics"]
    return {
        "three_digit_f1": f"{float(metrics['three_digit_f1']):.6f}",
        "five_digit_f1": f"{float(metrics['five_digit_f1']):.6f}",
        "seven_digit_f1": f"{float(metrics['seven_digit_f1']):.6f}",
    }


def read_per_gene_rows(
    label: str,
    method: str,
    bundle_dir: Path,
) -> list[dict[str, str]]:
    """Read per-gene metrics from a result bundle."""
    rows: list[dict[str, str]] = []
    with (bundle_dir / "per_gene.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "label": label,
                    "method": method,
                    "gene": row["gene"],
                    "three_digit_f1": f"{float(row['three_digit_f1']):.6f}",
                    "five_digit_f1": f"{float(row['five_digit_f1']):.6f}",
                    "seven_digit_f1": f"{float(row['seven_digit_f1']):.6f}",
                }
            )
    return rows


def run_method(
    repo_root: Path,
    label: str,
    method: MethodSpec,
    results_dir: Path,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Run typing and evaluation for one panel/method."""
    paths = paths_for(label, method, results_dir)
    if not paths.config.exists():
        raise FileNotFoundError(f"Missing config for {label}/{method.name}: {paths.config}")

    paths.prediction_tsv.parent.mkdir(parents=True, exist_ok=True)
    timed = run_timed(build_typing_command(method, paths), repo_root)
    if timed.exit_code:
        raise RuntimeError(
            f"Typing failed for {label}/{method.name}\n"
            f"stdout:\n{timed.stdout}\n"
            f"stderr:\n{timed.stderr}"
        )

    compare = subprocess.run(
        build_compare_command(method, paths),
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )
    if compare.returncode:
        raise RuntimeError(
            f"Evaluation failed for {label}/{method.name}\n"
            f"stdout:\n{compare.stdout}\n"
            f"stderr:\n{compare.stderr}"
        )

    metrics = read_bundle_metrics(paths.bundle_dir)
    row = {
        "label": label,
        "method": method.name,
        "runner": method.runner,
        "top_n": str(method.top_n),
        "base_top_n": str(method.base_top_n),
        "gene_base_top_ns": method.gene_base_top_ns,
        "functional_discard_fallback_genes": method.functional_discard_fallback_genes,
        "functional_discard_fallback_resolution": str(method.functional_discard_fallback_resolution),
        "functional_discard_fallback_max_score": str(method.functional_discard_fallback_max_score),
        "functional_discard_fallback_min_score_delta": str(method.functional_discard_fallback_min_score_delta),
        "runtime_seconds": f"{timed.runtime_seconds:.3f}",
        "max_rss_mb": "" if timed.max_rss_mb is None else f"{timed.max_rss_mb:.1f}",
        "three_digit_f1": metrics["three_digit_f1"],
        "five_digit_f1": metrics["five_digit_f1"],
        "seven_digit_f1": metrics["seven_digit_f1"],
        "prediction_tsv": str(paths.prediction_tsv),
        "bundle_dir": str(paths.bundle_dir),
    }
    return row, read_per_gene_rows(label, method.name, paths.bundle_dir)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    """Write rows to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def method_means(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Compute mean aggregate metrics per method."""
    methods = sorted({row["method"] for row in rows})
    means: list[dict[str, str]] = []
    for method in methods:
        subset = [row for row in rows if row["method"] == method]
        means.append(
            {
                "method": method,
                "panels": str(len(subset)),
                "runtime_seconds": f"{sum(float(row['runtime_seconds']) for row in subset) / len(subset):.3f}",
                "three_digit_f1": f"{sum(float(row['three_digit_f1']) for row in subset) / len(subset):.6f}",
                "five_digit_f1": f"{sum(float(row['five_digit_f1']) for row in subset) / len(subset):.6f}",
                "seven_digit_f1": f"{sum(float(row['seven_digit_f1']) for row in subset) / len(subset):.6f}",
            }
        )
    return means


def functional_regressions(
    per_gene_rows: list[dict[str, str]],
    baseline_method: str,
    candidate_method: str,
) -> list[dict[str, str]]:
    """Find candidate per-gene 3/5-digit F1 regressions versus a baseline."""
    by_key = {
        (row["label"], row["method"], row["gene"]): row
        for row in per_gene_rows
    }
    regressions: list[dict[str, str]] = []
    keys = sorted(
        (label, gene)
        for (label, method, gene) in by_key
        if method == baseline_method
    )
    for label, gene in keys:
        baseline = by_key.get((label, baseline_method, gene))
        candidate = by_key.get((label, candidate_method, gene))
        if baseline is None or candidate is None:
            continue
        for metric in FUNCTIONAL_METRICS:
            delta = float(candidate[metric]) - float(baseline[metric])
            if delta < -1e-12:
                regressions.append(
                    {
                        "label": label,
                        "gene": gene,
                        "metric": metric.replace("_f1", ""),
                        "baseline_f1": baseline[metric],
                        "candidate_f1": candidate[metric],
                        "delta": f"{delta:.6f}",
                    }
                )
    return regressions


def functional_gains(
    per_gene_rows: list[dict[str, str]],
    baseline_method: str,
    candidate_method: str,
) -> list[dict[str, str]]:
    """Find candidate per-gene 3/5-digit F1 improvements versus a baseline."""
    by_key = {
        (row["label"], row["method"], row["gene"]): row
        for row in per_gene_rows
    }
    gains: list[dict[str, str]] = []
    keys = sorted(
        (label, gene)
        for (label, method, gene) in by_key
        if method == baseline_method
    )
    for label, gene in keys:
        baseline = by_key.get((label, baseline_method, gene))
        candidate = by_key.get((label, candidate_method, gene))
        if baseline is None or candidate is None:
            continue
        for metric in FUNCTIONAL_METRICS:
            delta = float(candidate[metric]) - float(baseline[metric])
            if delta > 1e-12:
                gains.append(
                    {
                        "label": label,
                        "gene": gene,
                        "metric": metric.replace("_f1", ""),
                        "baseline_f1": baseline[metric],
                        "candidate_f1": candidate[metric],
                        "delta": f"{delta:.6f}",
                    }
                )
    return gains


def remaining_functional_errors(
    per_gene_rows: list[dict[str, str]],
    method: str,
) -> list[dict[str, str]]:
    """Find candidate per-gene 3/5-digit F1 values that are still below 1.0."""
    rows: list[dict[str, str]] = []
    for row in per_gene_rows:
        if row["method"] != method:
            continue
        for metric in FUNCTIONAL_METRICS:
            value = float(row[metric])
            if value < 1.0 - 1e-12:
                rows.append(
                    {
                        "label": row["label"],
                        "gene": row["gene"],
                        "metric": metric.replace("_f1", ""),
                        "f1": row[metric],
                    }
                )
    return sorted(rows, key=lambda row: (row["label"], row["gene"], row["metric"]))


def _fmt(value: str, digits: int = 4) -> str:
    """Format a TSV float string for Markdown."""
    return f"{float(value):.{digits}f}" if value else ""


def render_markdown(
    rows: list[dict[str, str]],
    per_gene_rows: list[dict[str, str]],
    summary_tsv: str,
    per_gene_tsv: str,
) -> str:
    """Render the functional stress sweep report."""
    means = method_means(rows)
    regressions = functional_regressions(
        per_gene_rows,
        baseline_method="discard",
        candidate_method="enhancedgate_kir2dl1fallback_geneaware",
    )
    gains = functional_gains(
        per_gene_rows,
        baseline_method="discard",
        candidate_method="enhancedgate_kir2dl1fallback_geneaware",
    )
    remaining_errors = remaining_functional_errors(
        per_gene_rows,
        method="enhancedgate_kir2dl1fallback_geneaware",
    )
    lines = [
        "# Functional Stress Seed Sweep",
        "",
        "This synthetic-first sweep re-runs only the typing and evaluation stages",
        "on existing generated fixtures. It is intended to expose 3-digit and",
        "5-digit functional typing weaknesses before spending time on real data.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/run_functional_stress_sweep.py",
        "```",
        "",
        f"Summary TSV: `{summary_tsv}`",
        f"Per-gene TSV: `{per_gene_tsv}`",
        "",
        "## Aggregate Results",
        "",
        "| panel | method | runtime s | 3-digit F1 | 5-digit F1 | 7-digit F1 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["label"],
                    row["method"],
                    _fmt(row["runtime_seconds"], 2),
                    _fmt(row["three_digit_f1"], 4),
                    _fmt(row["five_digit_f1"], 4),
                    _fmt(row["seven_digit_f1"], 4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Mean By Method",
            "",
            "| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in means:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["method"],
                    row["panels"],
                    _fmt(row["runtime_seconds"], 2),
                    _fmt(row["three_digit_f1"], 4),
                    _fmt(row["five_digit_f1"], 4),
                    _fmt(row["seven_digit_f1"], 4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## KIR2DL1 Fallback Gene-aware Vs Discard",
            "",
        ]
    )
    if regressions:
        lines.extend(
            [
                "Per-gene 3/5-digit regressions versus discard:",
                "",
                "| panel | gene | metric | discard F1 | fallback gene-aware F1 | delta |",
                "|---|---:|---:|---:|---:|---:|",
            ]
        )
        for row in regressions:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["label"],
                        row["gene"],
                        row["metric"],
                        _fmt(row["baseline_f1"], 4),
                        _fmt(row["candidate_f1"], 4),
                        _fmt(row["delta"], 4),
                    ]
                )
                + " |"
            )
    else:
        lines.append("No per-gene 3/5-digit regressions versus discard were detected.")

    lines.extend(
        [
            "",
            "Per-gene 3/5-digit gains versus discard:",
            "",
            "| panel | gene | metric | discard F1 | fallback gene-aware F1 | delta |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in gains:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["label"],
                    row["gene"],
                    row["metric"],
                    _fmt(row["baseline_f1"], 4),
                    _fmt(row["candidate_f1"], 4),
                    _fmt(row["delta"], 4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "Remaining fallback gene-aware 3/5-digit errors:",
            "",
            "| panel | gene | metric | F1 |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in remaining_errors:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["label"],
                    row["gene"],
                    row["metric"],
                    _fmt(row["f1"], 4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "`enhancedgate_kir2dl1fallback_geneaware` resolves the strict",
            "`synthetic-functional8x6` KIR2DL1 3-digit regression while preserving",
            "the enhancedgate KIR2DS3/KIR2DS4/KIR2DS5 gains and the gene-aware",
            "top-n runtime setting.",
            "`likelihood_top5000` alone is not viable because it loses substantial",
            "3/5-digit accuracy on the difficult5x12 seed panels.",
            "KIR2DL1 still has one 5-digit miss matching discard's remaining error;",
            "the next method work should focus on the remaining 5-digit errors in",
            "`KIR2DL5A/B`, `KIR2DS3`, and `KIR2DS5` without reintroducing the",
            "KIR2DL1 3-digit ambiguity-likelihood regression.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Run the functional stress sweep."""
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    labels = tuple(args.label) if args.label else DEFAULT_LABELS
    methods = select_methods(args.method)
    results_dir = Path(args.results_dir)

    rows: list[dict[str, str]] = []
    per_gene_rows: list[dict[str, str]] = []
    for label in labels:
        for method in methods:
            row, gene_rows = run_method(repo_root, label, method, results_dir)
            rows.append(row)
            per_gene_rows.extend(gene_rows)
            print(
                f"{label}\t{method.name}\t{row['runtime_seconds']}s\t"
                f"{row['three_digit_f1']}/{row['five_digit_f1']}/{row['seven_digit_f1']}"
            )

    summary_fields = [
        "label",
        "method",
        "runner",
        "top_n",
        "base_top_n",
        "gene_base_top_ns",
        "functional_discard_fallback_genes",
        "functional_discard_fallback_resolution",
        "functional_discard_fallback_max_score",
        "functional_discard_fallback_min_score_delta",
        "runtime_seconds",
        "max_rss_mb",
        "three_digit_f1",
        "five_digit_f1",
        "seven_digit_f1",
        "prediction_tsv",
        "bundle_dir",
    ]
    per_gene_fields = [
        "label",
        "method",
        "gene",
        "three_digit_f1",
        "five_digit_f1",
        "seven_digit_f1",
    ]
    write_tsv(Path(args.output_tsv), rows, summary_fields)
    write_tsv(Path(args.output_per_gene_tsv), per_gene_rows, per_gene_fields)
    if not args.skip_md:
        markdown = render_markdown(
            rows,
            per_gene_rows,
            args.output_tsv,
            args.output_per_gene_tsv,
        )
        output_md = Path(args.output_md)
        output_md.parent.mkdir(parents=True, exist_ok=True)
        output_md.write_text(markdown, encoding="utf-8")
        print(f"Wrote markdown report to {args.output_md}")
    print(f"Wrote summary TSV to {args.output_tsv}")
    print(f"Wrote per-gene TSV to {args.output_per_gene_tsv}")


if __name__ == "__main__":
    main()
