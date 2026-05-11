"""Inspect why one KIR2DS3 truth genotype is ranked far below the selected call."""

from __future__ import annotations

import argparse
import sys
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
)
from inspect_kir2ds3_remaining_failures import (  # noqa: E402
    KIR2DS3_BACKBONE,
    KIR2DS3_GENE,
    build_variant_rows_from_context,
    gene_alleles_list,
    resolve_kir2ds3_top_n,
)
from inspect_kir2ds5_remaining_failures import (  # noqa: E402
    summarize_variant_rows,
    write_tsv,
    write_variant_tsv,
)
from inspect_remaining_functional_errors import (  # noqa: E402
    DEFAULT_METHODS,
    load_panel_calls,
    read_sweep_summary,
)


RANK_FIELDNAMES = [
    "panel",
    "sample_id",
    "rank",
    "kind",
    "genotype",
    "three_digit",
    "five_digit",
    "value",
    "gap",
    "private_support",
    "has_truth_allele",
    "has_selected_allele",
    "has_kir2ds3_016",
    "has_kir2ds3_002",
    "truth_only_variants",
    "candidate_only_variants",
    "truth_only_exon",
    "candidate_only_exon",
    "truth_only_supported",
    "truth_only_unsupported",
    "candidate_only_supported",
    "candidate_only_unsupported",
    "truth_only_net",
    "candidate_only_net",
    "truth_minus_candidate_net",
    "ambiguous_positive_ratio",
    "signal",
]


def build_parser() -> argparse.ArgumentParser:
    """Build the KIR2DS3 rank-failure parser."""
    parser = argparse.ArgumentParser(
        description="Inspect KIR2DS3 sample-03 rank failure under likelihood typing.",
    )
    parser.add_argument(
        "--summary-tsv",
        default="benchmarks/results/functional-stress-sweep/summary.tsv",
        help="Functional stress sweep summary TSV.",
    )
    parser.add_argument(
        "--label",
        default="synthetic-difficult5x12",
        help="Panel label to inspect.",
    )
    parser.add_argument(
        "--sample-id",
        default="03",
        help="Sample ID to inspect.",
    )
    parser.add_argument(
        "--candidate-method",
        default="enhancedgate_kir2dl5_kir2ds5unsupported_geneaware",
        help="Candidate method whose selected call should be compared.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2ds3_rank_failure.tsv",
        help="Machine-readable genotype-rank TSV.",
    )
    parser.add_argument(
        "--variant-output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2ds3_rank_failure_variants.tsv",
        help="Machine-readable selected-vs-truth variant TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_kir2ds3_rank_failure.md",
        help="Markdown report path.",
    )
    parser.add_argument(
        "--top-detail-rows",
        type=int,
        default=30,
        help="Number of top likelihood genotypes to include before requested calls.",
    )
    parser.add_argument(
        "--max-variant-rows",
        type=int,
        default=30,
        help="Maximum selected-vs-truth variants to show in Markdown.",
    )
    return parser


def limit_allele_prefix(allele: str, resolution: int) -> str:
    """Limit an allele name to a fixed functional prefix."""
    if "*" not in allele:
        return allele
    gene, field = allele.split("*", 1)
    return f"{gene}*{field[:resolution]}"


def functional_key(alleles: Sequence[str], resolution: int) -> str:
    """Render a genotype at one functional resolution while preserving copies."""
    return allele_text([limit_allele_prefix(allele, resolution) for allele in alleles])


def genotype_key(alleles: Sequence[str]) -> tuple[str, ...]:
    """Return a stable exact genotype key."""
    return tuple(sorted(alleles))


def contains_allele_prefix(alleles: Sequence[str], prefix: str) -> bool:
    """Return whether any allele begins with a requested prefix."""
    return any(allele.startswith(prefix) for allele in alleles)


def rank_kind(
    alleles: Sequence[str],
    truth: Sequence[str],
    current: Sequence[str],
    discard: Sequence[str],
    likelihood: Sequence[str],
    rank: int,
    top_detail_rows: int,
) -> str:
    """Label why one genotype was included in the rank report."""
    allele_key = genotype_key(alleles)
    labels: list[str] = []
    if rank <= top_detail_rows:
        labels.append("top")
    if allele_key == genotype_key(truth):
        labels.append("truth")
    if allele_key == genotype_key(current):
        labels.append("current")
    if allele_key == genotype_key(discard):
        labels.append("discard")
    if allele_key == genotype_key(likelihood):
        labels.append("likelihood")
    return ";".join(labels) if labels else "requested"


def find_result_index(result: Any, alleles: Sequence[str]) -> int | None:
    """Find the zero-based exact rank index for a genotype."""
    if not alleles or result.isFail():
        return None
    target = genotype_key(alleles)
    for index, allele_names in enumerate(result.allele_name):
        if genotype_key(allele_names) == target:
            return index
    return None


def select_result_indexes(
    result: Any,
    requested_genotypes: Sequence[Sequence[str]],
    top_detail_rows: int,
) -> list[int]:
    """Return top indexes plus exact indexes for truth/current/context calls."""
    indexes = set(range(min(top_detail_rows, len(result.allele_name))))
    for genotype in requested_genotypes:
        index = find_result_index(result, genotype)
        if index is not None:
            indexes.add(index)
    return sorted(indexes)


def _variant_sort_key(variant: object) -> tuple[int, str]:
    """Sort variants by genomic position then ID."""
    return int(getattr(variant, "pos", 0)), str(getattr(variant, "id"))


def build_pair_variant_rows(
    truth: Sequence[str],
    candidate: Sequence[str],
    reads: Sequence[object],
    variants: Sequence[object],
) -> list[dict[str, str]]:
    """Build compact truth-vs-candidate variant rows for rank summaries."""
    variant_map = {str(getattr(variant, "id")): variant for variant in variants}
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
            rows.append(
                {
                    "side": side,
                    "in_exon": "yes" if getattr(variant, "in_exon") else "no",
                    "support_class": support_class(evidence),
                    "support_direction": support_direction(side, evidence),
                    "positive_weight": format_float(evidence.positive_weight),
                    "net_weight": format_float(evidence.net_weight),
                    "ambiguous_positive_weight": format_float(
                        evidence.ambiguous_positive_weight
                    ),
                }
            )
    return rows


def support_score_text(
    alleles: Sequence[str],
    allele_variants: dict[str, set[str]],
    positive: Any,
    negative: Any,
) -> str:
    """Compute private-support score when all alleles are model candidates."""
    if not alleles or any(allele not in allele_variants for allele in alleles):
        return ""
    from graphkir2.typing.private_support import private_support_score  # type: ignore[import-not-found]

    return format_float(private_support_score(list(alleles), allele_variants, positive, negative))


def run_likelihood_model(
    reads: list[Any],
    variants: list[Any],
    copy_number: int,
    top_n: int,
    ambiguity_neutral_prob: float,
) -> tuple[Any, Any, Any, dict[str, set[str]]]:
    """Run the KIR2DS3 likelihood model and return support context."""
    from graphkir.typing_mulit_allele import AlleleTyping
    from graphkir2.typing.private_support import collect_variant_support  # type: ignore[import-not-found]

    model = AlleleTyping(
        reads,
        variants,
        force_homo=False,
        top_n=top_n,
        variant_correction=True,
        exon_weight=2.0,
        ambiguity_likelihood=True,
        ambiguity_neutral_prob=ambiguity_neutral_prob,
    )
    result = model.typing(copy_number)
    positive, negative = collect_variant_support(model.reads)
    allele_variants = {
        allele: {
            str(variant.id)
            for variant in variants
            if allele in variant.allele
        }
        for allele in model.allele_to_id
    }
    return result, positive, negative, allele_variants


def build_rank_rows(
    label: str,
    sample_id: str,
    result: Any,
    truth: Sequence[str],
    current: Sequence[str],
    discard: Sequence[str],
    likelihood: Sequence[str],
    reads: Sequence[object],
    variants: Sequence[object],
    allele_variants: dict[str, set[str]],
    positive: Any,
    negative: Any,
    top_detail_rows: int,
) -> list[dict[str, str]]:
    """Build genotype-rank rows for top calls and requested context calls."""
    if result.isFail():
        return []
    indexes = select_result_indexes(
        result,
        [truth, current, discard, likelihood],
        top_detail_rows,
    )
    top_value = float(result.value[0])
    rows: list[dict[str, str]] = []
    for index in indexes:
        alleles = list(result.allele_name[index])
        summary = summarize_variant_rows(
            build_pair_variant_rows(truth, alleles, reads, variants)
        )
        rank = index + 1
        value = float(result.value[index])
        rows.append(
            {
                "panel": label,
                "sample_id": sample_id,
                "rank": str(rank),
                "kind": rank_kind(
                    alleles,
                    truth,
                    current,
                    discard,
                    likelihood,
                    rank,
                    top_detail_rows,
                ),
                "genotype": allele_text(alleles),
                "three_digit": functional_key(alleles, 3),
                "five_digit": functional_key(alleles, 5),
                "value": format_float(value),
                "gap": format_float(top_value - value),
                "private_support": support_score_text(
                    alleles,
                    allele_variants,
                    positive,
                    negative,
                ),
                "has_truth_allele": "yes"
                if any(allele in alleles for allele in truth)
                else "no",
                "has_selected_allele": "yes"
                if any(allele in alleles for allele in current)
                else "no",
                "has_kir2ds3_016": "yes"
                if contains_allele_prefix(alleles, "KIR2DS3*016")
                else "no",
                "has_kir2ds3_002": "yes"
                if contains_allele_prefix(alleles, "KIR2DS3*002")
                else "no",
                **summary,
            }
        )
    return rows


def selected_variant_rows(
    label: str,
    sample_id: str,
    truth: Sequence[str],
    current: Sequence[str],
    context: dict[str, Any],
) -> list[dict[str, str]]:
    """Build full selected-vs-truth variant rows for the Markdown detail."""
    return build_variant_rows_from_context(
        label,
        sample_id,
        "rank_failure",
        truth,
        current,
        context,
    )


def load_sample_calls(
    summary_rows: list[dict[str, str]],
    label: str,
    sample_id: str,
    candidate_method: str,
) -> tuple[list[str], list[str], list[str], list[str], dict[str, str]]:
    """Load truth/current/discard/likelihood calls for one sample."""
    methods: tuple[str, ...] = DEFAULT_METHODS
    if candidate_method not in methods:
        methods = (*methods, candidate_method)
    calls = load_panel_calls(summary_rows, methods)
    panel = calls[label]
    by_label_method = {
        (row["label"], row["method"]): row
        for row in summary_rows
    }
    truth = gene_alleles_list(panel["truth"].get(sample_id, []), KIR2DS3_GENE)
    current = gene_alleles_list(panel[candidate_method].get(sample_id, []), KIR2DS3_GENE)
    discard = gene_alleles_list(panel.get("discard", {}).get(sample_id, []), KIR2DS3_GENE)
    likelihood = gene_alleles_list(
        panel.get("likelihood_top5000", {}).get(sample_id, []),
        KIR2DS3_GENE,
    )
    return truth, current, discard, likelihood, by_label_method[(label, candidate_method)]


def build_outputs(
    summary_rows: list[dict[str, str]],
    label: str,
    sample_id: str,
    candidate_method: str,
    top_detail_rows: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str]]:
    """Build rank rows, variant rows, and sample diagnostics."""
    truth, current, discard, likelihood, summary_row = load_sample_calls(
        summary_rows,
        label,
        sample_id,
        candidate_method,
    )
    context, normalized_cn, config_path = load_typing_context(label, sample_id)
    reads_by_gene: dict[str, list[Any]] = context["likelihood_gene_reads"]
    variants_by_gene: dict[str, list[Any]] = context["gene_variants"]
    raw_gene_reads: dict[str, list[Any]] = context["raw_gene_reads"]
    discard_gene_reads: dict[str, list[Any]] = context["discard_gene_reads"]
    reads = reads_by_gene.get(KIR2DS3_BACKBONE, [])
    variants = variants_by_gene.get(KIR2DS3_BACKBONE, [])
    raw_reads = raw_gene_reads.get(KIR2DS3_BACKBONE, [])
    discard_reads = discard_gene_reads.get(KIR2DS3_BACKBONE, [])
    copy_number = normalized_cn.get(KIR2DS3_BACKBONE, 0)
    top_n = resolve_kir2ds3_top_n(summary_row)
    result, positive, negative, allele_variants = run_likelihood_model(
        reads,
        variants,
        copy_number,
        top_n,
        float(context["ambiguity_neutral_prob"]),
    )
    rank_rows = build_rank_rows(
        label,
        sample_id,
        result,
        truth,
        current,
        discard,
        likelihood,
        reads,
        variants,
        allele_variants,
        positive,
        negative,
        top_detail_rows,
    )
    variant_rows = selected_variant_rows(label, sample_id, truth, current, context)
    positive_variants, negative_variants = count_variant_observations(reads)
    diagnostics = {
        "config": str(config_path),
        "truth": allele_text(truth),
        "current": allele_text(current),
        "discard": allele_text(discard),
        "likelihood": allele_text(likelihood),
        "copy_number": str(copy_number),
        "top_n": str(top_n),
        "raw_reads": str(len(raw_reads)),
        "raw_informative_reads": str(sum(1 for read in raw_reads if has_read_evidence(read))),
        "likelihood_reads": str(len(reads)),
        "likelihood_informative_reads": str(
            sum(1 for read in reads if has_read_evidence(read))
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
    }
    return rank_rows, variant_rows, diagnostics


def _rank_row_by_kind(rows: Sequence[dict[str, str]], label: str) -> dict[str, str] | None:
    """Return the first rank row containing a kind label."""
    return next((row for row in rows if label in row["kind"].split(";")), None)


def render_markdown(
    rank_rows: list[dict[str, str]],
    variant_rows: list[dict[str, str]],
    diagnostics: dict[str, str],
    output_tsv: str,
    variant_output_tsv: str,
    max_variant_rows: int,
) -> str:
    """Render the KIR2DS3 rank-failure Markdown report."""
    truth_row = _rank_row_by_kind(rank_rows, "truth")
    current_row = _rank_row_by_kind(rank_rows, "current")
    selected_summary = summarize_variant_rows(variant_rows)
    lines = [
        "# KIR2DS3 Rank-Failure Audit",
        "",
        "This report focuses on the `synthetic-difficult5x12` sample `03` miss where",
        "`KIR2DS3*011_KIR2DS3*016` is the truth but current likelihood-based typing",
        "selects `KIR2DS3*011_KIR2DS3*0020101`.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/inspect_kir2ds3_rank_failure.py",
        "```",
        "",
        f"Rank TSV: `{output_tsv}`",
        f"Variant TSV: `{variant_output_tsv}`",
        "",
        "## Sample Context",
        "",
        "| field | value |",
        "|---|---|",
    ]
    for field in (
        "config",
        "truth",
        "current",
        "discard",
        "likelihood",
        "copy_number",
        "top_n",
        "likelihood_reads",
        "likelihood_informative_reads",
        "variants",
        "exon_variants",
        "positive_variant_ids",
        "negative_variant_ids",
    ):
        lines.append(f"| {field} | {diagnostics[field]} |")

    lines.extend(
        [
            "",
            "## Key Rank Rows",
            "",
            "| rank | kind | genotype | 3-digit | 5-digit | value | gap | private support | truth-only variants | candidate-only variants | candidate-only unsupported | truth-candidate margin | signal |",
            "|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    key_rows = [
        row
        for row in rank_rows
        if int(row["rank"]) <= 10
        or any(label in row["kind"].split(";") for label in ("truth", "current", "discard", "likelihood"))
    ]
    for row in key_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["rank"],
                    row["kind"],
                    row["genotype"],
                    row["three_digit"],
                    row["five_digit"],
                    row["value"],
                    row["gap"],
                    row["private_support"],
                    row["truth_only_variants"],
                    row["candidate_only_variants"],
                    row["candidate_only_unsupported"],
                    row["truth_minus_candidate_net"],
                    row["signal"],
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Selected-vs-Truth Variant Summary",
            "",
            "| truth-only variants | selected-only variants | truth-only exon | selected-only exon | truth-only unsupported | selected-only unsupported | truth net | selected net | truth-selected margin | ambiguous positive ratio | signal |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
            "| "
            + " | ".join(
                [
                    selected_summary["truth_only_variants"],
                    selected_summary["candidate_only_variants"],
                    selected_summary["truth_only_exon"],
                    selected_summary["candidate_only_exon"],
                    selected_summary["truth_only_unsupported"],
                    selected_summary["candidate_only_unsupported"],
                    selected_summary["truth_only_net"],
                    selected_summary["candidate_only_net"],
                    selected_summary["truth_minus_candidate_net"],
                    selected_summary["ambiguous_positive_ratio"],
                    selected_summary["signal"],
                ]
            )
            + " |",
            "",
            "## Top Selected-vs-Truth Variants",
            "",
            "| side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | selected carriers |",
            "|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in sorted(variant_rows, key=variant_detail_sort_key)[:max_variant_rows]:
        lines.append(
            "| "
            + " | ".join(
                [
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
    if truth_row and current_row:
        lines.extend(
            [
                f"* Truth rank is {truth_row['rank']} with likelihood gap {truth_row['gap']} and private support {truth_row['private_support']}.",
                f"* Current selected rank is {current_row['rank']} with private support {current_row['private_support']}.",
                f"* The selected-vs-truth variant margin is {selected_summary['truth_minus_candidate_net']}; positive values favor truth at variant-evidence level.",
                "* The miss is therefore not a near-window selection issue; `KIR2DS3*016` is too far below rank 1 for the existing top-window gates to rescue.",
                "* Any next experiment should test a very narrow `KIR2DS3*0020101` overcall penalty or pre-likelihood candidate pruning, then verify it does not regress the seed5102 sample where variant evidence favors the non-truth candidate.",
            ]
        )
    else:
        lines.append(
            "* The truth or current genotype was not present in the likelihood result list; inspect the TSV for top-n/copy-number filtering first."
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build the KIR2DS3 rank-failure audit report."""
    args = build_parser().parse_args()
    summary_rows = read_sweep_summary(Path(args.summary_tsv))
    rank_rows, variant_rows, diagnostics = build_outputs(
        summary_rows,
        args.label,
        args.sample_id,
        args.candidate_method,
        args.top_detail_rows,
    )
    write_tsv(Path(args.output_tsv), rank_rows)
    write_variant_tsv(Path(args.variant_output_tsv), variant_rows)
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(
        render_markdown(
            rank_rows,
            variant_rows,
            diagnostics,
            args.output_tsv,
            args.variant_output_tsv,
            args.max_variant_rows,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.variant_output_tsv}")
    print(f"Wrote {args.output_md}")


if __name__ == "__main__":
    main()
