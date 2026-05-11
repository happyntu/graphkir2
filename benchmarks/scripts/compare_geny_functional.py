"""Normalize and evaluate Geny HPRC predictions against graphkir2 metrics."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any

DEFAULT_TRUTH_TSV = "data/groundtruth/hprc_summary_v1_2_e.tsv"
DEFAULT_GENY_OUTPUT = "data/geny_hprc44.txt"
DEFAULT_OUTPUT_DIR = "benchmarks/results/geny-functional-comparison"


def build_parser() -> argparse.ArgumentParser:
    """Build the Geny functional comparison parser."""
    parser = argparse.ArgumentParser(
        description="Normalize Geny output and evaluate it at 3/5/7-digit resolution.",
    )
    parser.add_argument(
        "--truth-tsv",
        default=DEFAULT_TRUTH_TSV,
        help="Ground-truth allele TSV in graphkir merged-allele format.",
    )
    parser.add_argument(
        "--geny-output",
        default=DEFAULT_GENY_OUTPUT,
        help="Raw Geny output text file.",
    )
    parser.add_argument(
        "--graphkir-tsv",
        default="",
        help="Optional graphkir/graphkir2 prediction TSV to evaluate beside Geny.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for normalized Geny TSV and comparison reports.",
    )
    parser.add_argument(
        "--normalized-geny-tsv",
        default="",
        help="Optional explicit output path for normalized Geny allele TSV.",
    )
    return parser


def normalize_geny_allele(gene: str, raw_allele: str) -> str:
    """Normalize one Geny allele token into `GENE*fields` format."""
    clean_gene = gene.strip()
    token = raw_allele.strip().strip(",;")
    token = re.sub(r"[\$\#\+\=\-]+$", "", token)
    if not token:
        return ""
    if token.startswith("*") and clean_gene:
        token = f"{clean_gene}{token}"
    elif "*" not in token and clean_gene:
        token = f"{clean_gene}*{token}"
    if "*" not in token:
        return ""
    allele_gene, field = token.split("*", 1)
    field = field.split("e", 1)[0]
    field = re.sub(r"[^0-9A-Za-z]+$", "", field)
    if not allele_gene or not field:
        return ""
    return f"{allele_gene}*{field}"


def parse_geny_output(path: Path) -> dict[str, list[str]]:
    """Parse Geny text output into sample -> allele list."""
    results: dict[str, list[str]] = {}
    current_sample = ""
    current_alleles: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("Sample:"):
                if current_sample:
                    results[current_sample] = current_alleles
                current_sample = line.split("Sample:", 1)[1].strip()
                current_alleles = []
                continue
            parts = line.split()
            if (
                current_sample
                and len(parts) >= 3
                and parts[0] == "[ilp]"
                and parts[1].startswith("KIR")
            ):
                allele = normalize_geny_allele(parts[1], parts[2])
                if allele:
                    current_alleles.append(allele)
        if current_sample:
            results[current_sample] = current_alleles
    return results


def write_allele_tsv(sample_alleles: dict[str, list[str]], path: Path) -> None:
    """Write sample alleles in graphkir-compatible merged allele TSV format."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "alleles"], delimiter="\t")
        writer.writeheader()
        for sample_id in sorted(sample_alleles):
            writer.writerow(
                {
                    "id": sample_id,
                    "alleles": "_".join(sample_alleles[sample_id]),
                }
            )


def input_status_rows(
    truth_tsv: Path,
    geny_output: Path,
    graphkir_tsv: Path | None,
) -> list[dict[str, str]]:
    """Build input availability rows."""
    rows = [
        {
            "input": "truth_tsv",
            "path": str(truth_tsv),
            "required": "yes",
            "exists": str(truth_tsv.exists()).lower(),
            "action": "" if truth_tsv.exists() else "provide ground-truth TSV",
        },
        {
            "input": "geny_output",
            "path": str(geny_output),
            "required": "yes",
            "exists": str(geny_output.exists()).lower(),
            "action": "" if geny_output.exists() else "provide raw Geny output",
        },
    ]
    if graphkir_tsv is not None:
        rows.append(
            {
                "input": "graphkir_tsv",
                "path": str(graphkir_tsv),
                "required": "no",
                "exists": str(graphkir_tsv.exists()).lower(),
                "action": (
                    "" if graphkir_tsv.exists() else "optional comparison skipped"
                ),
            }
        )
    return rows


def write_input_status(rows: list[dict[str, str]], path: Path) -> None:
    """Write input availability rows as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["input", "path", "required", "exists", "action"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def has_missing_required(rows: list[dict[str, str]]) -> bool:
    """Return whether any required input is missing."""
    return any(row["required"] == "yes" and row["exists"] != "true" for row in rows)


def metric_row(tool: str, pred_path: Path, evaluation: Any) -> dict[str, str]:
    """Convert an AlleleEvaluationSummary-like object to a flat TSV row."""
    return {
        "tool": tool,
        "status": "evaluated",
        "pred_path": str(pred_path),
        "truth_samples": str(evaluation.truth_samples),
        "pred_samples": str(evaluation.pred_samples),
        "matched_samples": str(evaluation.matched_samples),
        "three_digit_precision": f"{evaluation.three_digit.precision:.6f}",
        "three_digit_recall": f"{evaluation.three_digit.recall:.6f}",
        "three_digit_f1": f"{evaluation.three_digit.f1:.6f}",
        "five_digit_precision": f"{evaluation.five_digit.precision:.6f}",
        "five_digit_recall": f"{evaluation.five_digit.recall:.6f}",
        "five_digit_f1": f"{evaluation.five_digit.f1:.6f}",
        "seven_digit_precision": f"{evaluation.seven_digit.precision:.6f}",
        "seven_digit_recall": f"{evaluation.seven_digit.recall:.6f}",
        "seven_digit_f1": f"{evaluation.seven_digit.f1:.6f}",
        "missing_truth_samples": ",".join(evaluation.missing_truth_samples),
        "missing_pred_samples": ",".join(evaluation.missing_pred_samples),
    }


def per_gene_rows(tool: str, evaluation: Any) -> list[dict[str, str]]:
    """Convert per-gene evaluation metrics to flat TSV rows."""
    rows: list[dict[str, str]] = []
    for gene in evaluation.per_gene:
        rows.append(
            {
                "tool": tool,
                "gene": gene.gene,
                "three_digit_precision": f"{gene.three_digit.precision:.6f}",
                "three_digit_recall": f"{gene.three_digit.recall:.6f}",
                "three_digit_f1": f"{gene.three_digit.f1:.6f}",
                "five_digit_precision": f"{gene.five_digit.precision:.6f}",
                "five_digit_recall": f"{gene.five_digit.recall:.6f}",
                "five_digit_f1": f"{gene.five_digit.f1:.6f}",
                "seven_digit_precision": f"{gene.seven_digit.precision:.6f}",
                "seven_digit_recall": f"{gene.seven_digit.recall:.6f}",
                "seven_digit_f1": f"{gene.seven_digit.f1:.6f}",
            }
        )
    return rows


def write_tsv(rows: list[dict[str, str]], fieldnames: list[str], path: Path) -> None:
    """Write rows as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def render_notes(
    input_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    normalized_geny_tsv: Path,
) -> str:
    """Render a short Markdown report."""
    lines = [
        "# Geny Functional Comparison",
        "",
        "This report normalizes raw Geny calls into graphkir-compatible allele TSV",
        "format and evaluates 3/5/7-digit precision, recall, and F1 using the",
        "same graphkir2 benchmark evaluator.",
        "",
        f"Normalized Geny TSV: `{normalized_geny_tsv}`",
        "",
        "## Input Status",
        "",
        "| input | required | exists | path | action |",
        "|---|---:|---:|---|---|",
    ]
    for row in input_rows:
        lines.append(
            f"| {row['input']} | {row['required']} | {row['exists']} | {row['path']} | {row['action']} |"
        )
    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| tool | matched samples | 3-digit F1 | 5-digit F1 | 7-digit F1 |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary_rows:
        lines.append(
            f"| {row['tool']} | {row['matched_samples']} | {row['three_digit_f1']} | {row['five_digit_f1']} | {row['seven_digit_f1']} |"
        )
    if not summary_rows:
        lines.extend(
            [
                "",
                "No comparison was evaluated because required inputs are missing.",
            ]
        )
    return "\n".join(lines) + "\n"


def run_comparison(
    truth_tsv: Path,
    geny_output: Path,
    graphkir_tsv: Path | None,
    output_dir: Path,
    normalized_geny_tsv: Path,
) -> int:
    """Run Geny normalization and optional graphkir comparison."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))
    from graphkir2.benchmark import evaluate_allele_calls  # type: ignore[import-not-found]

    output_dir.mkdir(parents=True, exist_ok=True)
    input_rows = input_status_rows(truth_tsv, geny_output, graphkir_tsv)
    write_input_status(input_rows, output_dir / "missing_inputs.tsv")
    summary_rows: list[dict[str, str]] = []
    gene_rows: list[dict[str, str]] = []
    if not has_missing_required(input_rows):
        geny_calls = parse_geny_output(geny_output)
        write_allele_tsv(geny_calls, normalized_geny_tsv)
        geny_eval = evaluate_allele_calls(str(truth_tsv), str(normalized_geny_tsv))
        summary_rows.append(metric_row("Geny", normalized_geny_tsv, geny_eval))
        gene_rows.extend(per_gene_rows("Geny", geny_eval))
        if graphkir_tsv is not None and graphkir_tsv.exists():
            graphkir_eval = evaluate_allele_calls(str(truth_tsv), str(graphkir_tsv))
            summary_rows.append(metric_row("graphkir2", graphkir_tsv, graphkir_eval))
            gene_rows.extend(per_gene_rows("graphkir2", graphkir_eval))
    summary_fields = [
        "tool",
        "status",
        "pred_path",
        "truth_samples",
        "pred_samples",
        "matched_samples",
        "three_digit_precision",
        "three_digit_recall",
        "three_digit_f1",
        "five_digit_precision",
        "five_digit_recall",
        "five_digit_f1",
        "seven_digit_precision",
        "seven_digit_recall",
        "seven_digit_f1",
        "missing_truth_samples",
        "missing_pred_samples",
    ]
    per_gene_fields = [
        "tool",
        "gene",
        "three_digit_precision",
        "three_digit_recall",
        "three_digit_f1",
        "five_digit_precision",
        "five_digit_recall",
        "five_digit_f1",
        "seven_digit_precision",
        "seven_digit_recall",
        "seven_digit_f1",
    ]
    write_tsv(summary_rows, summary_fields, output_dir / "summary.tsv")
    write_tsv(gene_rows, per_gene_fields, output_dir / "per_gene.tsv")
    (output_dir / "notes.md").write_text(
        render_notes(input_rows, summary_rows, normalized_geny_tsv),
        encoding="utf-8",
    )
    return 0


def main() -> None:
    """Run the Geny functional comparison adapter."""
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    normalized_geny_tsv = (
        Path(args.normalized_geny_tsv)
        if args.normalized_geny_tsv
        else output_dir / "geny.normalized.allele.tsv"
    )
    graphkir_tsv = Path(args.graphkir_tsv) if args.graphkir_tsv else None
    run_comparison(
        Path(args.truth_tsv),
        Path(args.geny_output),
        graphkir_tsv,
        output_dir,
        normalized_geny_tsv,
    )
    print(f"Wrote Geny comparison reports to {output_dir}")


if __name__ == "__main__":
    main()
