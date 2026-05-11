"""Inspect discriminating variants behind remaining KIR2DL5 failures."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from inspect_kir2dl5_remaining_failures import (  # noqa: E402
    KIR2DL5_BACKBONE,
    KIR2DL5_GENES,
    allele_text,
    classify_kir2dl5_cause,
    format_float,
    load_typing_context,
    merged_gene_alleles,
    read_remaining_rows,
)
from inspect_remaining_functional_errors import (  # noqa: E402
    DEFAULT_METHODS,
    load_panel_calls,
    read_sweep_summary,
)

FIELDNAMES = [
    "panel",
    "sample_id",
    "root_cause",
    "truth_full",
    "candidate_full",
    "side",
    "variant_id",
    "position",
    "variant_type",
    "variant_value",
    "in_exon",
    "carrier_allele_count",
    "truth_carriers",
    "candidate_carriers",
    "positive_weight",
    "negative_weight",
    "net_weight",
    "positive_reads",
    "negative_reads",
    "ambiguous_positive_weight",
    "ambiguous_negative_weight",
    "ambiguous_positive_ratio",
    "ambiguous_negative_ratio",
    "support_class",
    "support_direction",
]


@dataclass(frozen=True)
class VariantEvidence:
    """Weighted read evidence for one variant ID."""

    positive_weight: float
    negative_weight: float
    positive_reads: int
    negative_reads: int
    ambiguous_positive_weight: float
    ambiguous_negative_weight: float

    @property
    def net_weight(self) -> float:
        """Positive minus negative support."""
        return self.positive_weight - self.negative_weight

    @property
    def ambiguous_positive_ratio(self) -> float:
        """Fraction of positive support from ambiguous alignments."""
        if self.positive_weight <= 0.0:
            return 0.0
        return self.ambiguous_positive_weight / self.positive_weight

    @property
    def ambiguous_negative_ratio(self) -> float:
        """Fraction of negative support from ambiguous alignments."""
        if self.negative_weight <= 0.0:
            return 0.0
        return self.ambiguous_negative_weight / self.negative_weight


def build_parser() -> argparse.ArgumentParser:
    """Build the discriminating-variant audit parser."""
    parser = argparse.ArgumentParser(
        description="Inspect KIR2DL5 truth-vs-candidate discriminating variants.",
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
        default="enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware",
        help="Candidate method to inspect.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2dl5_discriminating_variants.tsv",
        help="Machine-readable discriminating-variant TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_kir2dl5_discriminating_variants.md",
        help="Markdown audit report path.",
    )
    parser.add_argument(
        "--max-detail-rows",
        type=int,
        default=18,
        help="Maximum discriminating variants to show in Markdown detail.",
    )
    return parser


def read_pair_name(read: object) -> str:
    """Extract the physical read-pair name from a PairRead-like object."""
    left = getattr(read, "l_sam", "")
    right = getattr(read, "r_sam", "")
    if left:
        return str(left).split("\t", 1)[0]
    if right:
        return str(right).split("\t", 1)[0]
    return ""


def read_weight(read: object) -> float:
    """Return the likelihood weight attached to a read alignment."""
    return float(getattr(read, "weight", 1.0))


def is_ambiguous_read(read: object) -> bool:
    """Return whether a read alignment came from an ambiguous read pair."""
    return float(getattr(read, "ambiguous_weight", 0.0)) > 0.0


def read_positive_variant_ids(read: object) -> set[str]:
    """Return positive variant IDs observed by a read."""
    return set(getattr(read, "lpv", []) + getattr(read, "rpv", []))


def read_negative_variant_ids(read: object) -> set[str]:
    """Return negative variant IDs observed by a read."""
    return set(getattr(read, "lnv", []) + getattr(read, "rnv", []))


def build_allele_variant_index(variants: Sequence[object]) -> dict[str, set[str]]:
    """Map allele name to variant IDs carried by that allele."""
    index: dict[str, set[str]] = defaultdict(set)
    for variant in variants:
        variant_id = str(getattr(variant, "id"))
        for allele in getattr(variant, "allele", []):
            index[str(allele)].add(variant_id)
    return dict(index)


def genotype_variant_ids(
    alleles: Sequence[str],
    allele_variants: dict[str, set[str]],
) -> set[str]:
    """Return the union of variant IDs carried by a genotype."""
    variant_ids: set[str] = set()
    for allele in alleles:
        variant_ids.update(allele_variants.get(allele, set()))
    return variant_ids


def discriminating_variant_ids(
    truth: Sequence[str],
    candidate: Sequence[str],
    allele_variants: dict[str, set[str]],
) -> tuple[set[str], set[str]]:
    """Return truth-only and candidate-only genotype-level variant IDs."""
    truth_ids = genotype_variant_ids(truth, allele_variants)
    candidate_ids = genotype_variant_ids(candidate, allele_variants)
    return truth_ids - candidate_ids, candidate_ids - truth_ids


def compute_variant_evidence(
    variant_id: str, reads: Sequence[object]
) -> VariantEvidence:
    """Compute weighted positive/negative support for one variant."""
    positive_weight = 0.0
    negative_weight = 0.0
    ambiguous_positive_weight = 0.0
    ambiguous_negative_weight = 0.0
    positive_names: set[str] = set()
    negative_names: set[str] = set()
    for read in reads:
        weight = read_weight(read)
        ambiguous = is_ambiguous_read(read)
        name = read_pair_name(read)
        if variant_id in read_positive_variant_ids(read):
            positive_weight += weight
            positive_names.add(name)
            if ambiguous:
                ambiguous_positive_weight += weight
        if variant_id in read_negative_variant_ids(read):
            negative_weight += weight
            negative_names.add(name)
            if ambiguous:
                ambiguous_negative_weight += weight
    return VariantEvidence(
        positive_weight=positive_weight,
        negative_weight=negative_weight,
        positive_reads=len(positive_names),
        negative_reads=len(negative_names),
        ambiguous_positive_weight=ambiguous_positive_weight,
        ambiguous_negative_weight=ambiguous_negative_weight,
    )


def support_class(evidence: VariantEvidence) -> str:
    """Classify one variant using the private-support score thresholds."""
    if evidence.positive_weight == 0.0 and evidence.negative_weight == 0.0:
        return "unobserved"
    if (
        evidence.positive_weight >= 2.0
        and evidence.positive_weight >= evidence.negative_weight
    ):
        return "supported"
    if evidence.negative_weight >= 5.0 and evidence.positive_weight == 0.0:
        return "unsupported"
    if evidence.positive_weight > evidence.negative_weight:
        return "positive_lean"
    if evidence.negative_weight > evidence.positive_weight:
        return "negative_lean"
    return "mixed"


def support_direction(side: str, evidence: VariantEvidence) -> str:
    """Explain whether this variant's evidence supports truth or candidate."""
    if evidence.positive_weight == 0.0 and evidence.negative_weight == 0.0:
        return "no_observed_signal"
    if abs(evidence.net_weight) < 1e-9:
        return "mixed"
    positive_owner = "truth" if side == "truth_only" else "candidate"
    negative_owner = "candidate" if side == "truth_only" else "truth"
    return (
        f"supports_{positive_owner}"
        if evidence.net_weight > 0
        else f"supports_{negative_owner}"
    )


def carriers_of_interest(variant: object, alleles: Sequence[str]) -> list[str]:
    """Return truth/candidate alleles that carry this variant."""
    variant_alleles = set(getattr(variant, "allele", []))
    return sorted(allele for allele in alleles if allele in variant_alleles)


def build_variant_rows_for_case(
    label: str,
    sample_id: str,
    root_cause: str,
    truth: Sequence[str],
    candidate: Sequence[str],
) -> list[dict[str, str]]:
    """Build discriminating-variant rows for one KIR2DL5 failure case."""
    context, _, _ = load_typing_context(label, sample_id)
    reads_by_gene: dict[str, list[Any]] = context["likelihood_gene_reads"]
    variants_by_gene: dict[str, list[Any]] = context["gene_variants"]
    reads = reads_by_gene.get(KIR2DL5_BACKBONE, [])
    variants = variants_by_gene.get(KIR2DL5_BACKBONE, [])
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
        for variant_id in sorted(variant_ids, key=lambda key: variant_map[key]):
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


def remaining_key(row: dict[str, str]) -> tuple[str, str]:
    """Build a sample key from a remaining-error row."""
    return row["panel"], row["sample_id"]


def build_discriminating_variant_rows(
    summary_rows: list[dict[str, str]],
    remaining_rows: list[dict[str, str]],
    candidate_method: str,
) -> list[dict[str, str]]:
    """Build discriminating-variant rows for remaining KIR2DL5 failures."""
    methods: tuple[str, ...] = DEFAULT_METHODS
    if candidate_method not in methods:
        methods = (*methods, candidate_method)
    calls = load_panel_calls(summary_rows, methods)
    rows_by_sample: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in remaining_rows:
        rows_by_sample.setdefault(remaining_key(row), []).append(row)

    variant_rows: list[dict[str, str]] = []
    for label, sample_id in sorted(rows_by_sample):
        panel = calls[label]
        truth = merged_gene_alleles(panel["truth"].get(sample_id, []))
        candidate = merged_gene_alleles(panel[candidate_method].get(sample_id, []))
        discard = merged_gene_alleles(panel.get("discard", {}).get(sample_id, []))
        likelihood = merged_gene_alleles(
            panel.get("likelihood_top5000", {}).get(sample_id, [])
        )
        root_cause = classify_kir2dl5_cause(truth, candidate, discard, likelihood)
        variant_rows.extend(
            build_variant_rows_for_case(label, sample_id, root_cause, truth, candidate)
        )
    return variant_rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write discriminating-variant rows as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _float(row: dict[str, str], field: str) -> float:
    """Parse a float field from an audit row."""
    return float(row[field] or "0")


def summarize_sample_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """Summarize discriminating variants by sample."""
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["panel"], row["sample_id"])].append(row)

    summaries: list[dict[str, str]] = []
    for (panel, sample_id), sample_rows in sorted(grouped.items()):
        truth_rows = [row for row in sample_rows if row["side"] == "truth_only"]
        candidate_rows = [row for row in sample_rows if row["side"] == "candidate_only"]
        truth_net = sum(_float(row, "net_weight") for row in truth_rows)
        candidate_net = sum(_float(row, "net_weight") for row in candidate_rows)
        margin = truth_net - candidate_net
        positive_weight = sum(_float(row, "positive_weight") for row in sample_rows)
        ambiguous_positive = sum(
            _float(row, "ambiguous_positive_weight") for row in sample_rows
        )
        ambiguous_ratio = (
            ambiguous_positive / positive_weight if positive_weight else 0.0
        )
        if margin > 0.0:
            signal = "variant_signal_favors_truth"
        elif margin < 0.0:
            signal = "variant_signal_favors_candidate"
        else:
            signal = "variant_signal_mixed"
        summaries.append(
            {
                "panel": panel,
                "sample_id": sample_id,
                "root_cause": sample_rows[0]["root_cause"],
                "truth_only_variants": str(len(truth_rows)),
                "candidate_only_variants": str(len(candidate_rows)),
                "truth_only_exon": str(
                    sum(1 for row in truth_rows if row["in_exon"] == "yes")
                ),
                "candidate_only_exon": str(
                    sum(1 for row in candidate_rows if row["in_exon"] == "yes")
                ),
                "truth_only_supported": str(
                    sum(1 for row in truth_rows if row["support_class"] == "supported")
                ),
                "candidate_only_supported": str(
                    sum(
                        1
                        for row in candidate_rows
                        if row["support_class"] == "supported"
                    )
                ),
                "truth_only_net": format_float(truth_net),
                "candidate_only_net": format_float(candidate_net),
                "truth_minus_candidate_net": format_float(margin),
                "ambiguous_positive_ratio": format_float(ambiguous_ratio),
                "signal": signal,
            }
        )
    return summaries


def variant_detail_sort_key(row: dict[str, str]) -> tuple[int, float, int]:
    """Sort discriminating variants for Markdown detail."""
    observed_rank = 0 if abs(_float(row, "net_weight")) > 0.0 else 1
    exon_rank = 0 if row["in_exon"] == "yes" else 1
    return observed_rank, -abs(_float(row, "net_weight")), exon_rank


def render_markdown(
    rows: list[dict[str, str]],
    output_tsv: str,
    max_detail_rows: int,
) -> str:
    """Render a concise Markdown report."""
    summaries = summarize_sample_rows(rows)
    direction_counts = Counter(row["support_direction"] for row in rows)
    lines = [
        "# KIR2DL5 Discriminating Variant Audit",
        "",
        "This report expands the remaining KIR2DL5 failures to genotype-level",
        "truth-only and candidate-only variants using the current likelihood",
        "multi-map weights.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/inspect_kir2dl5_discriminating_variants.py",
        "```",
        "",
        f"Full TSV: `{output_tsv}`",
        "",
        "## Sample Summary",
        "",
        "| panel | sample | root cause | truth-only variants | candidate-only variants | truth-only exon | candidate-only exon | truth net | candidate net | truth-candidate margin | ambiguous positive ratio | signal |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["sample_id"],
                    row["root_cause"],
                    row["truth_only_variants"],
                    row["candidate_only_variants"],
                    row["truth_only_exon"],
                    row["candidate_only_exon"],
                    row["truth_only_net"],
                    row["candidate_only_net"],
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
    for row in sorted(rows, key=variant_detail_sort_key)[:max_detail_rows]:
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
    if not rows:
        lines.extend(
            [
                "* No KIR2DL5 discriminating-variant rows are present because the current candidate has no remaining KIR2DL5A/B functional errors in this stress sweep.",
                "* The unsupported-overcall guard should remain constrained to KIR2DL5 and to cases with multiple unsupported selected-only variants.",
            ]
        )
    else:
        lines.extend(
            [
                "* KIR2DL5 failures should be interpreted from the truth-only and candidate-only variant evidence above.",
                "* A positive truth-candidate margin means observed discriminating-variant evidence favors the truth genotype despite the current candidate.",
                "* A negative margin means variant-level evidence itself agrees with the current candidate, so a broader guard may be unsafe.",
            ]
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build the KIR2DL5 discriminating-variant audit report."""
    args = build_parser().parse_args()
    summary_rows = read_sweep_summary(Path(args.summary_tsv))
    remaining_rows = read_remaining_rows(Path(args.remaining_tsv), KIR2DL5_GENES)
    rows = build_discriminating_variant_rows(
        summary_rows,
        remaining_rows,
        args.candidate_method,
    )
    write_tsv(Path(args.output_tsv), rows)
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(
        render_markdown(rows, args.output_tsv, args.max_detail_rows),
        encoding="utf-8",
    )
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.output_md}")


if __name__ == "__main__":
    main()
