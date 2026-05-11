"""Audit remaining KIR2DS3 functional typing failures."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from inspect_kir2dl5_discriminating_variants import (  # noqa: E402
    build_allele_variant_index,
    carriers_of_interest,
    compute_variant_evidence,
    discriminating_variant_ids,
    support_class,
    support_direction,
    variant_detail_sort_key,
)
from inspect_kir2dl5_remaining_failures import (  # noqa: E402
    allele_text,
    count_variant_observations,
    format_float,
    has_read_evidence,
    load_typing_context,
    rank_evidence,
    read_remaining_rows,
)
from inspect_kir2ds5_remaining_failures import (  # noqa: E402
    VARIANT_FIELDNAMES,
    summarize_variant_rows,
    write_tsv,
    write_variant_tsv,
)
from inspect_remaining_functional_errors import (  # noqa: E402
    DEFAULT_METHODS,
    gene_name,
    load_panel_calls,
    read_sweep_summary,
)


KIR2DS3_GENE = "KIR2DS3"
KIR2DS3_BACKBONE = "KIR2DS3*BACKBONE"


def build_parser() -> argparse.ArgumentParser:
    """Build the KIR2DS3 failure audit parser."""
    parser = argparse.ArgumentParser(
        description="Audit KIR2DS3 remaining functional typing failures.",
    )
    parser.add_argument(
        "--summary-tsv",
        default="benchmarks/results/functional-stress-sweep/summary.tsv",
        help="Functional stress sweep summary TSV.",
    )
    parser.add_argument(
        "--remaining-tsv",
        default="benchmarks/results/functional-stress-sweep/remaining_functional_errors.tsv",
        help="Sample-level remaining functional error TSV.",
    )
    parser.add_argument(
        "--candidate-method",
        default="enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware",
        help="Candidate method to audit.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2ds3_failure_audit.tsv",
        help="Machine-readable KIR2DS3 audit TSV.",
    )
    parser.add_argument(
        "--variant-output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2ds3_discriminating_variants.tsv",
        help="Machine-readable KIR2DS3 discriminating-variant TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_kir2ds3_failure_audit.md",
        help="Markdown audit report path.",
    )
    parser.add_argument(
        "--max-detail-rows",
        type=int,
        default=24,
        help="Maximum discriminating variants to show in Markdown detail.",
    )
    return parser


def gene_alleles_list(
    alleles: Sequence[str],
    gene: str = KIR2DS3_GENE,
) -> list[str]:
    """Return full allele names belonging to one gene."""
    return [allele for allele in alleles if gene_name(allele) == gene]


def gene_copy_text(alleles: Sequence[str]) -> str:
    """Render the KIR2DS3 copy count."""
    return f"{KIR2DS3_GENE}={len(gene_alleles_list(alleles))}"


def classify_kir2ds3_cause(
    truth: Sequence[str],
    candidate: Sequence[str],
    discard: Sequence[str] | None = None,
    likelihood: Sequence[str] | None = None,
) -> str:
    """Classify remaining KIR2DS3 failures at single-gene scope."""
    truth_key = sorted(truth)
    candidate_key = sorted(candidate)
    discard_key = sorted(discard) if discard is not None else None
    likelihood_key = sorted(likelihood) if likelihood is not None else None
    if len(truth_key) != len(candidate_key):
        return "copy_count_mismatch"
    if truth_key == candidate_key:
        return "resolved"
    if discard_key == truth_key:
        return "candidate_regression"
    if discard_key == candidate_key:
        return "shared_with_discard"
    if likelihood_key == candidate_key:
        return "unresolved_likelihood_pattern"
    return "all_methods_disagree_or_shifted"


def resolve_kir2ds3_top_n(summary_row: dict[str, str]) -> int:
    """Resolve the effective KIR2DS3 top-n for the private-support candidate."""
    return int(summary_row.get("top_n") or "5000")


def _variant_sort_key(variant: object) -> tuple[int, str]:
    """Sort variants by genomic position then ID."""
    return int(getattr(variant, "pos", 0)), str(getattr(variant, "id"))


def build_variant_rows_from_context(
    label: str,
    sample_id: str,
    root_cause: str,
    truth: Sequence[str],
    candidate: Sequence[str],
    context: dict[str, Any],
) -> list[dict[str, str]]:
    """Build discriminating-variant rows for one KIR2DS3 failure case."""
    reads_by_gene: dict[str, list[Any]] = context["likelihood_gene_reads"]
    variants_by_gene: dict[str, list[Any]] = context["gene_variants"]
    reads = reads_by_gene.get(KIR2DS3_BACKBONE, [])
    variants = variants_by_gene.get(KIR2DS3_BACKBONE, [])
    variant_map = {str(variant.id): variant for variant in variants}
    allele_variants = build_allele_variant_index(variants)
    truth_only, candidate_only = discriminating_variant_ids(
        truth,
        candidate,
        allele_variants,
    )

    rows: list[dict[str, str]] = []
    for side, variant_ids in (
        ("truth_only", truth_only),
        ("candidate_only", candidate_only),
    ):
        for variant_id in sorted(variant_ids, key=lambda key: _variant_sort_key(variant_map[key])):
            variant = variant_map[variant_id]
            evidence = compute_variant_evidence(variant_id, reads)
            truth_carriers = carriers_of_interest(variant, truth)
            candidate_carriers = carriers_of_interest(variant, candidate)
            rows.append(
                {
                    "panel": label,
                    "sample_id": sample_id,
                    "root_cause": root_cause,
                    "truth_full": allele_text(truth),
                    "candidate_full": allele_text(candidate),
                    "side": side,
                    "variant_id": variant_id,
                    "position": str(variant.pos),
                    "variant_type": str(variant.typ),
                    "variant_value": str(variant.val),
                    "in_exon": "yes" if variant.in_exon else "no",
                    "carrier_allele_count": str(len(variant.allele)),
                    "truth_carriers": allele_text(truth_carriers),
                    "candidate_carriers": allele_text(candidate_carriers),
                    "positive_weight": format_float(evidence.positive_weight),
                    "negative_weight": format_float(evidence.negative_weight),
                    "net_weight": format_float(evidence.net_weight),
                    "positive_reads": str(evidence.positive_reads),
                    "negative_reads": str(evidence.negative_reads),
                    "ambiguous_positive_weight": format_float(
                        evidence.ambiguous_positive_weight
                    ),
                    "ambiguous_negative_weight": format_float(
                        evidence.ambiguous_negative_weight
                    ),
                    "ambiguous_positive_ratio": format_float(
                        evidence.ambiguous_positive_ratio
                    ),
                    "ambiguous_negative_ratio": format_float(
                        evidence.ambiguous_negative_ratio
                    ),
                    "support_class": support_class(evidence),
                    "support_direction": support_direction(side, evidence),
                }
            )
    return rows


def inspect_kir2ds3_sample(
    label: str,
    sample_id: str,
    summary_row: dict[str, str],
    root_cause: str,
    truth: Sequence[str],
    candidate: Sequence[str],
) -> tuple[dict[str, str], list[dict[str, str]]]:
    """Inspect read/variant/rank diagnostics for one remaining KIR2DS3 sample."""
    context, normalized_cn, config_path = load_typing_context(label, sample_id)
    raw_gene_reads: dict[str, list[Any]] = context["raw_gene_reads"]
    likelihood_gene_reads: dict[str, list[Any]] = context["likelihood_gene_reads"]
    discard_gene_reads: dict[str, list[Any]] = context["discard_gene_reads"]
    gene_variants: dict[str, list[Any]] = context["gene_variants"]
    raw_cn: dict[str, int] = context["raw_cn"]

    raw_reads = raw_gene_reads.get(KIR2DS3_BACKBONE, [])
    likelihood_reads = likelihood_gene_reads.get(KIR2DS3_BACKBONE, [])
    discard_reads = discard_gene_reads.get(KIR2DS3_BACKBONE, [])
    variants = gene_variants.get(KIR2DS3_BACKBONE, [])
    copy_number = normalized_cn.get(KIR2DS3_BACKBONE, 0)
    top_n = resolve_kir2ds3_top_n(summary_row)

    likelihood_evidence = rank_evidence(
        likelihood_reads,
        variants,
        copy_number,
        top_n,
        exon_weight=2.0,
        ambiguity_likelihood=True,
        ambiguity_neutral_prob=float(context["ambiguity_neutral_prob"]),
        min_fraction_ratio=float(context["select_min_fraction_ratio"]),
        truth=truth,
        candidate=candidate,
    )
    discard_evidence = rank_evidence(
        discard_reads,
        variants,
        copy_number,
        top_n=600,
        exon_weight=1.0,
        ambiguity_likelihood=False,
        ambiguity_neutral_prob=0.999,
        min_fraction_ratio=0.5,
        truth=truth,
        candidate=candidate,
    )
    positive_variants, negative_variants = count_variant_observations(likelihood_reads)
    variant_rows = build_variant_rows_from_context(
        label,
        sample_id,
        root_cause,
        truth,
        candidate,
        context,
    )

    diagnostics = {
        "config": str(config_path),
        "typing_backbone": KIR2DS3_BACKBONE,
        "raw_cn": str(raw_cn.get(KIR2DS3_BACKBONE, 0)),
        "normalized_cn": str(copy_number),
        "current_top_n": str(top_n),
        "raw_reads": str(len(raw_reads)),
        "raw_informative_reads": str(sum(1 for read in raw_reads if has_read_evidence(read))),
        "likelihood_reads": str(len(likelihood_reads)),
        "likelihood_informative_reads": str(
            sum(1 for read in likelihood_reads if has_read_evidence(read))
        ),
        "discard_reads": str(len(discard_reads)),
        "discard_informative_reads": str(
            sum(1 for read in discard_reads if has_read_evidence(read))
        ),
        "variants": str(len(variants)),
        "exon_variants": str(sum(1 for variant in variants if variant.in_exon)),
        "candidate_alleles": str(
            len({allele for variant in variants for allele in variant.allele})
        ),
        "positive_variant_ids": str(positive_variants),
        "negative_variant_ids": str(negative_variants),
        **{f"likelihood_{key}": value for key, value in likelihood_evidence.items()},
        **{f"discard_{key}": value for key, value in discard_evidence.items()},
        **summarize_variant_rows(variant_rows),
    }
    return diagnostics, variant_rows


def remaining_key(row: dict[str, str]) -> tuple[str, str]:
    """Build a sample key from a remaining-error row."""
    return row["panel"], row["sample_id"]


def build_audit_outputs(
    summary_rows: list[dict[str, str]],
    remaining_rows: list[dict[str, str]],
    candidate_method: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build KIR2DS3 sample-level and variant-level audit rows."""
    methods: tuple[str, ...] = DEFAULT_METHODS
    if candidate_method not in methods:
        methods = (*methods, candidate_method)
    calls = load_panel_calls(summary_rows, methods)
    by_label_method = {
        (row["label"], row["method"]): row
        for row in summary_rows
    }
    rows_by_sample: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in remaining_rows:
        rows_by_sample.setdefault(remaining_key(row), []).append(row)

    audit_rows: list[dict[str, str]] = []
    variant_rows: list[dict[str, str]] = []
    for (label, sample_id), sample_errors in sorted(rows_by_sample.items()):
        summary_row = by_label_method[(label, candidate_method)]
        panel = calls[label]
        truth = gene_alleles_list(panel["truth"].get(sample_id, []))
        candidate = gene_alleles_list(panel[candidate_method].get(sample_id, []))
        discard = gene_alleles_list(panel.get("discard", {}).get(sample_id, []))
        likelihood = gene_alleles_list(
            panel.get("likelihood_top5000", {}).get(sample_id, [])
        )
        cause = classify_kir2ds3_cause(truth, candidate, discard, likelihood)
        diagnostics, sample_variant_rows = inspect_kir2ds3_sample(
            label,
            sample_id,
            summary_row,
            cause,
            truth,
            candidate,
        )
        variant_rows.extend(sample_variant_rows)
        audit_rows.append(
            {
                "panel": label,
                "sample_id": sample_id,
                "root_cause": cause,
                "remaining_error_rows": ";".join(
                    sorted(
                        f"{row['gene']}:{row['resolution']}:{row['cause_hint']}"
                        for row in sample_errors
                    )
                ),
                "truth_full": allele_text(truth),
                "candidate_full": allele_text(candidate),
                "discard_full": allele_text(discard),
                "likelihood_full": allele_text(likelihood),
                "truth_gene_count": gene_copy_text(truth),
                "candidate_gene_count": gene_copy_text(candidate),
                "discard_gene_count": gene_copy_text(discard),
                "likelihood_gene_count": gene_copy_text(likelihood),
                **diagnostics,
            }
        )
    return audit_rows, variant_rows


def render_markdown(
    audit_rows: list[dict[str, str]],
    variant_rows: list[dict[str, str]],
    output_tsv: str,
    variant_output_tsv: str,
    max_detail_rows: int,
) -> str:
    """Render the KIR2DS3 audit Markdown report."""
    cause_counts = Counter(row["root_cause"] for row in audit_rows)
    signal_counts = Counter(row["signal"] for row in audit_rows)
    direction_counts = Counter(row["support_direction"] for row in variant_rows)
    lines = [
        "# KIR2DS3 Remaining Failure Audit",
        "",
        "This report audits the remaining KIR2DS3 functional errors for the",
        "current functional-stress lead and expands each miss into likelihood",
        "rank/support diagnostics plus truth-only vs candidate-only variant evidence.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/inspect_kir2ds3_remaining_failures.py",
        "```",
        "",
        f"Sample TSV: `{output_tsv}`",
        f"Variant TSV: `{variant_output_tsv}`",
        "",
        "## Summary",
        "",
        "| root cause | samples |",
        "|---|---:|",
    ]
    for cause, count in sorted(cause_counts.items()):
        lines.append(f"| {cause} | {count} |")

    lines.extend(
        [
            "",
            "| variant signal | samples |",
            "|---|---:|",
        ]
    )
    for signal, count in sorted(signal_counts.items()):
        lines.append(f"| {signal} | {count} |")

    lines.extend(
        [
            "",
            "## Sample Detail",
            "",
            "| panel | sample | root cause | truth | candidate | discard | likelihood selected | truth rank | candidate rank | truth gap | candidate gap | truth support | candidate support | truth-only variants | candidate-only variants | truth-candidate margin | ambiguous positive ratio | signal |",
            "|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in audit_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["sample_id"],
                    row["root_cause"],
                    row["truth_full"],
                    row["candidate_full"],
                    row["discard_full"],
                    row["likelihood_selected"],
                    row["likelihood_truth_rank"] or "not in top-n",
                    row["likelihood_candidate_rank"] or "not in top-n",
                    row["likelihood_truth_value_gap"] or "",
                    row["likelihood_candidate_value_gap"] or "",
                    row["likelihood_truth_support"] or "",
                    row["likelihood_candidate_support"] or "",
                    row["truth_only_variants"],
                    row["candidate_only_variants"],
                    row["truth_minus_candidate_net"],
                    row["ambiguous_positive_ratio"],
                    row["signal"],
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Direction Counts",
            "",
            "| support direction | variants |",
            "|---|---:|",
        ]
    )
    for direction, count in sorted(direction_counts.items()):
        lines.append(f"| {direction} | {count} |")

    lines.extend(
        [
            "",
            "## Top Discriminating Variants",
            "",
            "| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |",
            "|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in sorted(variant_rows, key=variant_detail_sort_key)[:max_detail_rows]:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["sample_id"],
                    row["side"],
                    row["variant_id"],
                    row["position"],
                    row["in_exon"],
                    row["support_class"],
                    row["support_direction"],
                    row["positive_weight"],
                    row["negative_weight"],
                    row["net_weight"],
                    row["ambiguous_positive_ratio"],
                    row["truth_carriers"],
                    row["candidate_carriers"],
                ]
            )
            + " |"
        )

    lines.extend(["", "## Interpretation", ""])
    if not audit_rows:
        lines.extend(
            [
                "* No KIR2DS3 3-digit or 5-digit remaining errors are present for the current candidate.",
                "* Keep KIR2DS3 private-support and fallback gates unchanged unless broader stress panels reveal a candidate regression.",
            ]
        )
    else:
        candidate_signal = signal_counts.get("variant_signal_favors_candidate", 0)
        truth_signal = signal_counts.get("variant_signal_favors_truth", 0)
        mixed_signal = signal_counts.get("variant_signal_mixed", 0)
        candidate_regressions = cause_counts.get("candidate_regression", 0)
        unsupported_candidate_kir2ds3_002 = sum(
            1
            for row in variant_rows
            if row["side"] == "candidate_only"
            and row["support_class"] == "unsupported"
            and "KIR2DS3*002" in row["candidate_carriers"]
        )
        unsupported_truth_only = sum(
            1
            for row in variant_rows
            if row["side"] == "truth_only" and row["support_class"] == "unsupported"
        )
        lines.extend(
            [
                f"* Variant evidence favors truth in {truth_signal} sample(s), candidate in {candidate_signal} sample(s), and is mixed in {mixed_signal} sample(s).",
                f"* Candidate-regression KIR2DS3 samples: {candidate_regressions}.",
                f"* Unsupported candidate-only `KIR2DS3*002` variants appear in {unsupported_candidate_kir2ds3_002} row(s).",
                f"* Unsupported truth-only variants appear in {unsupported_truth_only} row(s), so the residual set is not uniformly truth-rescuable.",
                "* Treat `shared_with_discard` rows separately from all-method shifted rows; shared rows need new evidence, not rollback of current guards.",
                "* Any next KIR2DS3 method should be narrow and preserve the current KIR2DS5/KIR2DL5 gains.",
            ]
        )
        if truth_signal and not candidate_signal:
            lines.append(
                "* A truth-favoring KIR2DS3 guard may be testable, but only after checking which allele prefix and discriminating variants recur across the rows above."
            )
        else:
            lines.append(
                "* Do not add a KIR2DS3 guard from aggregate F1 alone; the sample-level variant signal is not uniformly truth-favoring."
            )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build the KIR2DS3 remaining-failure audit report."""
    args = build_parser().parse_args()
    summary_rows = read_sweep_summary(Path(args.summary_tsv))
    remaining_rows = read_remaining_rows(Path(args.remaining_tsv), [KIR2DS3_GENE])
    audit_rows, variant_rows = build_audit_outputs(
        summary_rows,
        remaining_rows,
        args.candidate_method,
    )
    write_tsv(Path(args.output_tsv), audit_rows)
    write_variant_tsv(Path(args.variant_output_tsv), variant_rows)
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(
        render_markdown(
            audit_rows,
            variant_rows,
            args.output_tsv,
            args.variant_output_tsv,
            args.max_detail_rows,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.variant_output_tsv}")
    print(f"Wrote {args.output_md}")


if __name__ == "__main__":
    main()
