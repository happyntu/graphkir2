"""Summarize synthetic method comparisons and wrong calls."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MethodBundle:
    """One evaluated method bundle for a synthetic panel."""

    panel: str
    method: str
    bundle_dir: str


@dataclass(frozen=True)
class CandidatePrediction:
    """Truth/prediction pair used for wrong-call reporting."""

    panel: str
    method: str
    truth_tsv: str
    pred_tsv: str


DEFAULT_METHODS: tuple[MethodBundle, ...] = (
    MethodBundle("synthetic-functional8", "discard", "benchmarks/results/synthetic-functional8-discard-bundle"),
    MethodBundle("synthetic-functional8", "best-only+exon2", "benchmarks/results/synthetic-functional8-bestonly-exon2-bundle"),
    MethodBundle("synthetic-functional8", "margin+exon2", "benchmarks/results/synthetic-functional8-margin-g5-s2-bundle"),
    MethodBundle("synthetic-functional8", "likelihood+exon2", "benchmarks/results/synthetic-functional8-likelihood-exon2-bundle"),
    MethodBundle("synthetic-functional8", "likelihood+exon2+minfrac0.7+top5000", "benchmarks/results/synthetic-functional8-likelihood-exon2-minfrac07-top5000-bundle"),
    MethodBundle("synthetic-functional8", "functional-target+kir2ds3-private", "benchmarks/results/synthetic-functional8-private-support-kir2ds3-group_neutral-bundle"),
    MethodBundle("synthetic-functional8x6", "discard", "benchmarks/results/synthetic-functional8x6-discard-bundle"),
    MethodBundle("synthetic-functional8x6", "best-only+exon2", "benchmarks/results/synthetic-functional8x6-bestonly-exon2-bundle"),
    MethodBundle("synthetic-functional8x6", "margin+exon2", "benchmarks/results/synthetic-functional8x6-margin-g5-s2-bundle"),
    MethodBundle("synthetic-functional8x6", "likelihood+exon2", "benchmarks/results/synthetic-functional8x6-likelihood-exon2-bundle"),
    MethodBundle("synthetic-functional8x6", "likelihood+exon2+minfrac0.7+top5000", "benchmarks/results/synthetic-functional8x6-likelihood-exon2-minfrac07-top5000-bundle"),
    MethodBundle("synthetic-functional8x6", "functional-target+kir2ds3-private", "benchmarks/results/synthetic-functional8x6-private-support-kir2ds3-group_neutral-bundle"),
    MethodBundle("synthetic-difficult5", "discard", "benchmarks/results/synthetic-difficult5-discard-bundle"),
    MethodBundle("synthetic-difficult5", "best-only+exon2", "benchmarks/results/synthetic-difficult5-bestonly-exon2-bundle"),
    MethodBundle("synthetic-difficult5", "margin+exon2", "benchmarks/results/synthetic-difficult5-margin-g5-s2-bundle"),
    MethodBundle("synthetic-difficult5", "likelihood+exon2", "benchmarks/results/synthetic-difficult5-likelihood-exon2-bundle"),
    MethodBundle("synthetic-difficult5", "likelihood+exon2+minfrac0.7+top5000", "benchmarks/results/synthetic-difficult5-likelihood-exon2-minfrac07-top5000-bundle"),
    MethodBundle("synthetic-difficult5", "functional-target+kir2ds3-private", "benchmarks/results/synthetic-difficult5-private-support-kir2ds3-group_neutral-bundle"),
)

DEFAULT_CANDIDATES: tuple[CandidatePrediction, ...] = (
    CandidatePrediction(
        "synthetic-functional8",
        "likelihood+exon2+minfrac0.7+top5000",
        "benchmarks/generated/synthetic-functional8/synthetic-functional8_truth.tsv",
        "benchmarks/results/synthetic-functional8/likelihood_exon2_minfrac07_top5000.allele.tsv",
    ),
    CandidatePrediction(
        "synthetic-functional8x6",
        "likelihood+exon2+minfrac0.7+top5000",
        "benchmarks/generated/synthetic-functional8x6/synthetic-functional8x6_truth.tsv",
        "benchmarks/results/synthetic-functional8x6/likelihood_exon2_minfrac07_top5000.allele.tsv",
    ),
    CandidatePrediction(
        "synthetic-difficult5",
        "likelihood+exon2+minfrac0.7+top5000",
        "benchmarks/generated/synthetic-difficult5/synthetic-difficult5_truth.tsv",
        "benchmarks/results/synthetic-difficult5/likelihood_exon2_minfrac07_top5000.allele.tsv",
    ),
)


def build_parser() -> argparse.ArgumentParser:
    """Build the summary report parser."""
    parser = argparse.ArgumentParser(
        description="Build synthetic method comparison and wrong-call reports.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_synthetic_method_comparison.md",
        help="Markdown report path.",
    )
    parser.add_argument(
        "--output-comparison-tsv",
        default="benchmarks/results/synthetic-method-comparison.tsv",
        help="Machine-readable method comparison TSV path.",
    )
    parser.add_argument(
        "--output-wrong-calls-tsv",
        default="benchmarks/results/synthetic-likelihood-wrong-calls.tsv",
        help="Wrong-call TSV path for the current candidate baseline.",
    )
    return parser


def read_summary(bundle_dir: str) -> dict[str, str]:
    """Read a bundle summary TSV and return its single data row."""
    path = Path(bundle_dir) / "summary.tsv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one summary row in {path}, got {len(rows)}")
    return rows[0]


def read_comparison_rows(methods: tuple[MethodBundle, ...]) -> list[dict[str, str]]:
    """Collect selected metrics from method bundle summaries."""
    rows: list[dict[str, str]] = []
    for method in methods:
        summary = read_summary(method.bundle_dir)
        rows.append(
            {
                "panel": method.panel,
                "method": method.method,
                "three_digit_f1": summary["three_digit_f1"],
                "five_digit_f1": summary["five_digit_f1"],
                "seven_digit_f1": summary["seven_digit_f1"],
                "runtime_seconds": summary["runtime_seconds"],
                "max_rss_mb": summary["max_rss_mb"],
                "bundle_dir": method.bundle_dir,
            }
        )
    return rows


def write_tsv(path: str, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    """Write dictionary rows as TSV."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _format_metric(value: str) -> str:
    """Format a metric string for compact Markdown display."""
    if not value:
        return ""
    return f"{float(value):.5g}"


def _extract_sample_id(raw: str) -> str:
    """Match graphkir2 eval sample-id extraction for merged prediction names."""
    name = Path(raw).name
    parts = name.split(".")
    if len(parts) >= 2 and parts[1].isdigit():
        return parts[1]
    return raw.strip()


def _normalize_allele(raw: str) -> str:
    """Normalize an allele token enough for wrong-call comparison."""
    import re

    token = raw.strip()
    if not token or "*" not in token:
        return ""
    token = re.sub(r"[\$\#\+\=\-]+$", "", token)
    gene, field = token.split("*", 1)
    field = field.split("e", 1)[0]
    field = re.sub(r"[^0-9A-Za-z]+$", "", field)
    if not gene or not field:
        return ""
    return f"{gene}*{field}"


def _limit_allele(allele: str, resolution: int) -> str:
    """Limit allele field length to the requested KIR resolution."""
    gene, field = allele.split("*", 1)
    return f"{gene}*{field[:resolution]}"


def _gene_name(allele: str) -> str:
    """Return the gene portion of an allele name."""
    return allele.split("*", 1)[0]


def read_allele_rows(path: str) -> dict[str, list[str]]:
    """Read truth or prediction TSV into sample -> alleles."""
    rows: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Allele file is missing a header: {path}")
        for row in reader:
            sample_id = (
                row.get("id")
                or row.get("sample_id")
                or row.get("name")
                or ""
            ).strip()
            if not sample_id:
                continue
            if "name" in row and not row.get("id") and not row.get("sample_id"):
                sample_id = _extract_sample_id(sample_id)
            allele_field = row.get("alleles", "").strip()
            alleles = [_normalize_allele(token) for token in allele_field.split("_")]
            rows[sample_id] = [allele for allele in alleles if allele]
    return rows


def build_gene_counter(alleles: list[str], resolution: int) -> dict[str, Counter[str]]:
    """Build per-gene allele counters at one resolution."""
    counters: dict[str, Counter[str]] = defaultdict(Counter)
    for allele in alleles:
        counters[_gene_name(allele)][_limit_allele(allele, resolution)] += 1
    return counters


def build_wrong_call_rows(
    candidates: tuple[CandidatePrediction, ...],
) -> list[dict[str, str]]:
    """Build per-sample, per-gene wrong-call rows for candidate predictions."""
    rows: list[dict[str, str]] = []
    for candidate in candidates:
        truth = read_allele_rows(candidate.truth_tsv)
        pred = read_allele_rows(candidate.pred_tsv)
        sample_ids = sorted(set(truth) | set(pred))
        for sample_id in sample_ids:
            for resolution in (3, 5, 7):
                truth_genes = build_gene_counter(truth.get(sample_id, []), resolution)
                pred_genes = build_gene_counter(pred.get(sample_id, []), resolution)
                for gene in sorted(set(truth_genes) | set(pred_genes)):
                    truth_counter = truth_genes.get(gene, Counter())
                    pred_counter = pred_genes.get(gene, Counter())
                    if truth_counter == pred_counter:
                        continue
                    missing = truth_counter - pred_counter
                    extra = pred_counter - truth_counter
                    rows.append(
                        {
                            "panel": candidate.panel,
                            "method": candidate.method,
                            "sample_id": sample_id,
                            "gene": gene,
                            "resolution": str(resolution),
                            "truth": "_".join(truth_counter.elements()),
                            "pred": "_".join(pred_counter.elements()),
                            "missing": "_".join(missing.elements()),
                            "extra": "_".join(extra.elements()),
                        }
                    )
    return rows


def summarize_wrong_calls(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize wrong calls by panel/gene/resolution."""
    counter: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        counter[(row["panel"], row["gene"], row["resolution"])] += 1
    return [
        {
            "panel": panel,
            "gene": gene,
            "resolution": resolution,
            "wrong_call_rows": str(count),
        }
        for (panel, gene, resolution), count in sorted(counter.items())
    ]


def render_markdown(
    comparison_rows: list[dict[str, str]],
    wrong_rows: list[dict[str, str]],
    comparison_tsv: str,
    wrong_calls_tsv: str,
) -> str:
    """Render the synthetic method report as Markdown."""
    lines = [
        "# 2026-05-11 Synthetic Method Comparison",
        "",
        "## Scope",
        "",
        "This report compares the current synthetic multi-map candidates across the",
        "`synthetic-functional8`, `synthetic-functional8x6`, and `synthetic-difficult5`",
        "panels. It is generated from existing benchmark bundles and can be rebuilt with:",
        "",
        "```bash",
        "python benchmarks/scripts/summarize_synthetic_methods.py",
        "```",
        "",
        "## Method Comparison",
        "",
        "| panel | method | 3-digit F1 | 5-digit F1 | 7-digit F1 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in comparison_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["method"],
                    _format_metric(row["three_digit_f1"]),
                    _format_metric(row["five_digit_f1"]),
                    _format_metric(row["seven_digit_f1"]),
                ]
            )
            + " |"
        )

    wrong_summary = summarize_wrong_calls(wrong_rows)
    lines.extend(
        [
            "",
            "## Candidate Wrong-Call Summary",
            "",
            "The wrong-call table uses the current candidate baseline:",
            "`likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 + top_n=5000`.",
            "",
            f"Full comparison TSV: `{comparison_tsv}`",
            f"Full wrong-call TSV: `{wrong_calls_tsv}`",
            "",
            "| panel | gene | resolution | wrong-call rows |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in wrong_summary:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["gene"],
                    row["resolution"],
                    row["wrong_call_rows"],
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Current Decision",
            "",
            "`likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 + top_n=5000` is the current synthetic candidate baseline.",
            "It is the first tested multi-map method that improves `5-digit` and `7-digit`",
            "on `synthetic-functional8x6` while preserving the discard baseline `3-digit`.",
            "Raising `top_n` from `600` to `5000` improves `synthetic-functional8x6`",
            "without regressing `synthetic-functional8` or `synthetic-difficult5`.",
            "On `synthetic-difficult5`, the same candidate now also normalizes",
            "`KIR2DL5A/B` CN hints onto the merged `KIR2DL5*BACKBONE` graph key,",
            "which removes the previous missing-call failure mode.",
            "",
            "A separate functional-target candidate adds targeted `KIR2DS3` private-support",
            "reranking plus directional `KIR2DS3/KIR2DS5` cross-gene neutralization.",
            "It preserves `synthetic-functional8` and `synthetic-functional8x6` at",
            "the same `3/5-digit` scores as the top5000 baseline and raises",
            "`synthetic-difficult5` to `3/5/7-digit = 1.0` when combined with",
            "the `KIR2DS4` same-5-digit highest-suffix tie-break. Treat it as the",
            "current functional target.",
            "",
            "Remaining failure focus:",
            "",
            "* balanced baseline: `KIR2DL1` on `synthetic-functional8x6`",
            "* balanced baseline: `KIR2DS3` functional typing on `synthetic-difficult5`",
            "* functional target: no remaining synthetic-difficult5 error in this panel",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build comparison and wrong-call reports."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))
    args = build_parser().parse_args()

    comparison_rows = read_comparison_rows(DEFAULT_METHODS)
    wrong_rows = build_wrong_call_rows(DEFAULT_CANDIDATES)

    comparison_fields = [
        "panel",
        "method",
        "three_digit_f1",
        "five_digit_f1",
        "seven_digit_f1",
        "runtime_seconds",
        "max_rss_mb",
        "bundle_dir",
    ]
    wrong_fields = [
        "panel",
        "method",
        "sample_id",
        "gene",
        "resolution",
        "truth",
        "pred",
        "missing",
        "extra",
    ]
    write_tsv(args.output_comparison_tsv, comparison_rows, comparison_fields)
    write_tsv(args.output_wrong_calls_tsv, wrong_rows, wrong_fields)

    markdown = render_markdown(
        comparison_rows,
        wrong_rows,
        args.output_comparison_tsv,
        args.output_wrong_calls_tsv,
    )
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(markdown, encoding="utf-8")
    print(f"Wrote {args.output_md}")
    print(f"Wrote {args.output_comparison_tsv}")
    print(f"Wrote {args.output_wrong_calls_tsv}")


if __name__ == "__main__":
    main()
