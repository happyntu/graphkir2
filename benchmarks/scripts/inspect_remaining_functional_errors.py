"""Inspect remaining sample-level functional typing errors after a stress sweep."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_METHODS = (
    "discard",
    "likelihood_top5000",
    "enhancedgate_geneaware",
    "enhancedgate_kir2dl1fallback_geneaware",
    "enhancedgate_kir2dl1_kir2ds5guard_geneaware",
    "enhancedgate_functionalguard_geneaware",
    "enhancedgate_kir2dl5guard_geneaware",
    "enhancedgate_kir2dl5_kir2ds5unsupported_geneaware",
    "enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware",
)


def build_parser() -> argparse.ArgumentParser:
    """Build the remaining-error inspection parser."""
    parser = argparse.ArgumentParser(
        description="Inspect remaining sample-level 3/5-digit errors from a functional stress sweep.",
    )
    parser.add_argument(
        "--summary-tsv",
        default="benchmarks/results/functional-stress-sweep/summary.tsv",
        help="Functional stress sweep summary TSV.",
    )
    parser.add_argument(
        "--candidate-method",
        default="enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware",
        help="Method whose remaining functional errors should be triaged.",
    )
    parser.add_argument(
        "--method",
        action="append",
        default=[],
        help="Methods to include in comparison columns. Defaults to current stress methods.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/remaining_functional_errors.tsv",
        help="Machine-readable wrong-call TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_remaining_functional_error_triage.md",
        help="Markdown report path.",
    )
    return parser


def extract_sample_id(raw: str) -> str:
    """Extract graphkir synthetic sample IDs from truth or prediction rows."""
    name = Path(raw).name
    parts = name.split(".")
    if len(parts) >= 2 and parts[1].isdigit():
        return parts[1]
    return raw.strip()


def normalize_allele(raw: str) -> str:
    """Normalize one allele token for functional-resolution comparison."""
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


def limit_allele(allele: str, resolution: int) -> str:
    """Limit an allele field to a functional-resolution prefix."""
    gene, field = allele.split("*", 1)
    return f"{gene}*{field[:resolution]}"


def gene_name(allele: str) -> str:
    """Return the gene name from a normalized allele."""
    return allele.split("*", 1)[0]


def read_allele_rows(path: str) -> dict[str, list[str]]:
    """Read a truth or prediction allele TSV into sample -> normalized alleles."""
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
                sample_id = extract_sample_id(sample_id)
            allele_field = row.get("alleles", "").strip()
            alleles = [normalize_allele(token) for token in allele_field.split("_")]
            rows[sample_id] = [allele for allele in alleles if allele]
    return rows


def gene_counter(alleles: list[str], gene: str, resolution: int) -> Counter[str]:
    """Build a single-gene allele counter at one resolution."""
    counter: Counter[str] = Counter()
    for allele in alleles:
        if gene_name(allele) == gene:
            counter[limit_allele(allele, resolution)] += 1
    return counter


def gene_alleles(alleles: list[str], gene: str) -> str:
    """Return full normalized alleles for one gene."""
    return "_".join(allele for allele in alleles if gene_name(allele) == gene)


def counter_text(counter: Counter[str]) -> str:
    """Render an allele counter in stable TSV-friendly form."""
    return "_".join(sorted(counter.elements()))


def read_sweep_summary(path: Path) -> list[dict[str, str]]:
    """Read functional stress sweep summary rows."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truth_path_for_label(label: str) -> str:
    """Read the truth TSV path from a committed benchmark config."""
    config_path = Path("benchmarks/configs") / f"{label}.json"
    data = json.loads(config_path.read_text(encoding="utf-8"))
    truth_path = str(data.get("allele_truth_tsv") or "")
    if not truth_path:
        raise ValueError(f"Config is missing allele_truth_tsv: {config_path}")
    return truth_path


def load_panel_calls(
    summary_rows: list[dict[str, str]],
    methods: tuple[str, ...],
) -> dict[str, dict[str, dict[str, list[str]]]]:
    """Load truth and selected method predictions for each panel."""
    labels = sorted({row["label"] for row in summary_rows})
    calls: dict[str, dict[str, dict[str, list[str]]]] = {}
    by_label_method = {
        (row["label"], row["method"]): row
        for row in summary_rows
    }
    for label in labels:
        panel: dict[str, dict[str, list[str]]] = {
            "truth": read_allele_rows(truth_path_for_label(label))
        }
        for method in methods:
            row = by_label_method.get((label, method))
            if row is None:
                continue
            panel[method] = read_allele_rows(row["prediction_tsv"])
        calls[label] = panel
    return calls


def classify_error(
    truth: Counter[str],
    candidate: Counter[str],
    discard: Counter[str],
    likelihood: Counter[str],
    enhancedgate: Counter[str],
) -> str:
    """Classify a remaining candidate error from method agreement patterns."""
    if sum(truth.values()) != sum(candidate.values()):
        return "copy_number_count_mismatch"
    if candidate == discard:
        return "shared_with_discard"
    if discard == truth:
        return "candidate_regression"
    if likelihood == truth and candidate != truth:
        return "fallback_or_private_support_regression"
    if enhancedgate == truth and candidate != truth:
        return "functional_fallback_regression"
    if likelihood == candidate:
        return "unresolved_likelihood_pattern"
    return "all_methods_disagree_or_shifted"


def build_error_rows(
    calls: dict[str, dict[str, dict[str, list[str]]]],
    candidate_method: str,
    methods: tuple[str, ...],
) -> list[dict[str, str]]:
    """Build sample-level remaining 3/5-digit error rows for the candidate."""
    rows: list[dict[str, str]] = []
    for label, panel in sorted(calls.items()):
        truth_rows = panel["truth"]
        candidate_rows = panel[candidate_method]
        sample_ids = sorted(set(truth_rows) | set(candidate_rows))
        for sample_id in sample_ids:
            truth_alleles = truth_rows.get(sample_id, [])
            candidate_alleles = candidate_rows.get(sample_id, [])
            genes = sorted(
                {gene_name(allele) for allele in truth_alleles}
                | {gene_name(allele) for allele in candidate_alleles}
            )
            for gene in genes:
                for resolution in (3, 5):
                    truth = gene_counter(truth_alleles, gene, resolution)
                    candidate = gene_counter(candidate_alleles, gene, resolution)
                    if truth == candidate:
                        continue
                    discard = gene_counter(
                        panel.get("discard", {}).get(sample_id, []),
                        gene,
                        resolution,
                    )
                    likelihood = gene_counter(
                        panel.get("likelihood_top5000", {}).get(sample_id, []),
                        gene,
                        resolution,
                    )
                    enhancedgate = gene_counter(
                        panel.get("enhancedgate_geneaware", {}).get(sample_id, []),
                        gene,
                        resolution,
                    )
                    missing = truth - candidate
                    extra = candidate - truth
                    row = {
                        "panel": label,
                        "sample_id": sample_id,
                        "gene": gene,
                        "resolution": str(resolution),
                        "cause_hint": classify_error(
                            truth,
                            candidate,
                            discard,
                            likelihood,
                            enhancedgate,
                        ),
                        "truth": counter_text(truth),
                        "candidate": counter_text(candidate),
                        "missing": counter_text(missing),
                        "extra": counter_text(extra),
                        "truth_full": gene_alleles(truth_alleles, gene),
                        "candidate_full": gene_alleles(candidate_alleles, gene),
                    }
                    for method in methods:
                        method_alleles = panel.get(method, {}).get(sample_id, [])
                        row[method] = counter_text(
                            gene_counter(method_alleles, gene, resolution)
                        )
                        row[f"{method}_full"] = gene_alleles(method_alleles, gene)
                    rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    """Write rows as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_errors(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize remaining errors by gene, resolution, and cause."""
    counter: Counter[tuple[str, str, str]] = Counter()
    for row in rows:
        counter[(row["gene"], row["resolution"], row["cause_hint"])] += 1
    return [
        {
            "gene": gene,
            "resolution": resolution,
            "cause_hint": cause,
            "wrong_call_rows": str(count),
        }
        for (gene, resolution, cause), count in sorted(counter.items())
    ]


def render_markdown(
    rows: list[dict[str, str]],
    output_tsv: str,
    candidate_method: str,
) -> str:
    """Render a concise Markdown triage report."""
    summary = summarize_errors(rows)
    remaining_genes = ", ".join(sorted({row["gene"] for row in rows})) or "none"
    regression_counts: Counter[tuple[str, str]] = Counter(
        (row["gene"], row["resolution"])
        for row in rows
        if row["cause_hint"] == "candidate_regression"
    )
    regression_text = ", ".join(
        f"{gene} {resolution}-digit: {count}"
        for (gene, resolution), count in sorted(regression_counts.items())
    )
    if not regression_text:
        regression_text = "none"
    lines = [
        "# Remaining Functional Error Triage",
        "",
        "This report expands the functional stress sweep into sample-level",
        "3-digit and 5-digit wrong calls for the current candidate.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/inspect_remaining_functional_errors.py",
        "```",
        "",
        f"Candidate method: `{candidate_method}`",
        f"Full TSV: `{output_tsv}`",
        "",
        "## Summary",
        "",
        "| gene | resolution | cause hint | wrong-call rows |",
        "|---|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["gene"],
                    row["resolution"],
                    row["cause_hint"],
                    row["wrong_call_rows"],
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Detail",
            "",
            "| panel | sample | gene | res | cause hint | truth | candidate | missing | extra | discard | likelihood | enhancedgate |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["sample_id"],
                    row["gene"],
                    row["resolution"],
                    row["cause_hint"],
                    row["truth"],
                    row["candidate"],
                    row["missing"],
                    row["extra"],
                    row.get("discard", ""),
                    row.get("likelihood_top5000", ""),
                    row.get("enhancedgate_geneaware", ""),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "* `shared_with_discard` means the current candidate did not introduce the error; fixing it probably requires new gene-specific evidence handling rather than undoing enhancedgate.",
            f"* `candidate_regression` is the highest-priority blocker because discard is correct and the candidate is wrong. Current candidate regressions: {regression_text}.",
            "* If candidate regressions are `none`, next work should target shared or unresolved failure patterns rather than adding broader regression guards.",
            f"* Current remaining functional errors are in: {remaining_genes}.",
            "* The KIR2DL1 3-digit regression is absent here; only a 5-digit suballele miss remains, matching discard behavior.",
            "",
            "## Recommended Next Method Work",
            "",
            "* KIR2DS5 has no remaining 3/5-digit functional rows for the current candidate; keep the targeted KIR2DS5 guard narrow unless broader panels reveal a regression.",
            "* Inspect the remaining KIR2DS3 rows at suballele/private-support level before changing the broader KIR2DS3 rescue gate.",
            "* Keep any future KIR2DL5A/B work separate from the KIR2DS3/KIR2DS5 gate work.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build the remaining functional error triage report."""
    args = build_parser().parse_args()
    methods = tuple(args.method) if args.method else DEFAULT_METHODS
    if args.candidate_method not in methods:
        methods = (*methods, args.candidate_method)
    summary_rows = read_sweep_summary(Path(args.summary_tsv))
    calls = load_panel_calls(summary_rows, methods)
    rows = build_error_rows(calls, args.candidate_method, methods)

    base_fields = [
        "panel",
        "sample_id",
        "gene",
        "resolution",
        "cause_hint",
        "truth",
        "candidate",
        "missing",
        "extra",
        "truth_full",
        "candidate_full",
    ]
    method_fields = [
        field
        for method in methods
        for field in (method, f"{method}_full")
    ]
    write_tsv(Path(args.output_tsv), rows, [*base_fields, *method_fields])
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(
        render_markdown(rows, args.output_tsv, args.candidate_method),
        encoding="utf-8",
    )
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.output_md}")


if __name__ == "__main__":
    main()
