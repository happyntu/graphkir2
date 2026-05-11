"""Audit remaining KIR2DL5A/B functional typing failures."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from copy import deepcopy
from pathlib import Path
from typing import Any, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from inspect_remaining_functional_errors import (  # noqa: E402
    DEFAULT_METHODS,
    extract_sample_id,
    gene_name,
    load_panel_calls,
    read_sweep_summary,
)


KIR2DL5_GENES = ("KIR2DL5A", "KIR2DL5B")
KIR2DL5_BACKBONE = "KIR2DL5*BACKBONE"


def build_parser() -> argparse.ArgumentParser:
    """Build the KIR2DL5 failure audit parser."""
    parser = argparse.ArgumentParser(
        description="Audit KIR2DL5A/B remaining functional typing failures.",
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
        default="enhancedgate_kir2dl5guard_geneaware",
        help="Candidate method to audit.",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/functional-stress-sweep/kir2dl5_failure_audit.tsv",
        help="Machine-readable KIR2DL5 audit TSV.",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_kir2dl5_failure_audit.md",
        help="Markdown audit report path.",
    )
    return parser


def split_allele_text(raw: str) -> list[str]:
    """Split a Graph-KIR allele field while preserving duplicate alleles."""
    return [token for token in raw.split("_") if token]


def allele_text(alleles: Sequence[str]) -> str:
    """Render alleles in stable TSV-friendly order."""
    return "_".join(sorted(alleles))


def typing_backbone_alias(gene_or_backbone: str) -> str:
    """Return the merged typing backbone for KIR2DL5A/B."""
    if gene_or_backbone.startswith("KIR2DL5A") or gene_or_backbone.startswith("KIR2DL5B"):
        return KIR2DL5_BACKBONE
    return gene_or_backbone


def merged_gene_alleles(
    alleles: Sequence[str],
    genes: Sequence[str] = KIR2DL5_GENES,
) -> list[str]:
    """Return alleles belonging to the merged KIR2DL5A/B audit group."""
    gene_set = set(genes)
    return [allele for allele in alleles if gene_name(allele) in gene_set]


def gene_copy_counts(
    alleles: Sequence[str],
    genes: Sequence[str] = KIR2DL5_GENES,
) -> Counter[str]:
    """Count alleles by pure gene name for the audit group."""
    counts: Counter[str] = Counter({gene: 0 for gene in genes})
    for allele in alleles:
        gene = gene_name(allele)
        if gene in counts:
            counts[gene] += 1
    return counts


def gene_copy_text(counts: Counter[str]) -> str:
    """Render KIR2DL5A/B/merged copy counts."""
    merged = sum(counts[gene] for gene in KIR2DL5_GENES)
    return f"KIR2DL5A={counts['KIR2DL5A']};KIR2DL5B={counts['KIR2DL5B']};merged={merged}"


def classify_kir2dl5_cause(
    truth: Sequence[str],
    candidate: Sequence[str],
    discard: Sequence[str] | None = None,
    likelihood: Sequence[str] | None = None,
) -> str:
    """Classify remaining KIR2DL5A/B failures at merged-backbone scope."""
    truth_key = sorted(truth)
    candidate_key = sorted(candidate)
    if len(truth_key) != len(candidate_key):
        return "merged_copy_count_mismatch"
    if gene_copy_counts(truth) != gene_copy_counts(candidate):
        return "ab_assignment_mismatch_on_merged_backbone"
    if truth_key == candidate_key:
        return "resolved"
    if discard is not None and candidate_key == sorted(discard):
        return "shared_allele_substitution"
    if likelihood is not None and candidate_key == sorted(likelihood):
        return "unresolved_likelihood_substitution"
    return "allele_substitution"


def parse_gene_top_n_spec(spec: str) -> dict[str, int]:
    """Parse comma-separated `GENE:TOPN` overrides."""
    values: dict[str, int] = {}
    for field in spec.split(","):
        item = field.strip()
        if not item:
            continue
        gene, sep, raw_value = item.partition(":")
        if not sep:
            raise ValueError(f"Invalid gene top-n field: {item}")
        values[gene.strip()] = int(raw_value.strip())
    return values


def resolve_gene_top_n(gene: str, summary_row: dict[str, str]) -> int:
    """Resolve the effective current-candidate top-n for one gene."""
    top_n = int(summary_row.get("top_n") or "600")
    base_top_n = int(summary_row.get("base_top_n") or "0")
    gene_base_top_ns = parse_gene_top_n_spec(summary_row.get("gene_base_top_ns", ""))
    if gene in gene_base_top_ns:
        return min(top_n, gene_base_top_ns[gene])
    if base_top_n > 0:
        return base_top_n
    return top_n


def read_remaining_rows(path: Path, genes: Sequence[str]) -> list[dict[str, str]]:
    """Read remaining-error rows for selected genes."""
    gene_set = set(genes)
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row for row in rows if row.get("gene", "") in gene_set]


def find_candidate_rank(result: Any, alleles: Sequence[str]) -> tuple[str, str, str]:
    """Return 1-based rank, value, and top-value gap for one allele tuple."""
    if not alleles or result.isFail():
        return "", "", ""
    target = tuple(sorted(alleles))
    top_value = float(result.value[0])
    for index, allele_names in enumerate(result.allele_name):
        if tuple(sorted(allele_names)) == target:
            value = float(result.value[index])
            return str(index + 1), format_float(value), format_float(top_value - value)
    return "", "", ""


def format_float(value: float) -> str:
    """Format floats compactly for TSV/Markdown reports."""
    return f"{value:.3f}"


def has_read_evidence(read: object) -> bool:
    """Return whether a read has positive or negative variant observations."""
    return bool(read.lpv + read.lnv + read.rpv + read.rnv)  # type: ignore[attr-defined]


def count_variant_observations(reads: Sequence[object]) -> tuple[int, int]:
    """Count unique positive and negative variant IDs observed by reads."""
    positive: set[str] = set()
    negative: set[str] = set()
    for read in reads:
        positive.update(read.lpv + read.rpv)  # type: ignore[attr-defined]
        negative.update(read.lnv + read.rnv)  # type: ignore[attr-defined]
    return len(positive), len(negative)


def support_score_text(
    alleles: Sequence[str],
    allele_variants: dict[str, set[str]],
    positive: Any,
    negative: Any,
) -> str:
    """Compute private-support score when all requested alleles are model candidates."""
    if not alleles or any(allele not in allele_variants for allele in alleles):
        return ""
    from graphkir2.typing.private_support import private_support_score  # type: ignore[import-not-found]

    return format_float(private_support_score(list(alleles), allele_variants, positive, negative))


def rank_evidence(
    reads: list[Any],
    variants: list[Any],
    copy_number: int,
    top_n: int,
    exon_weight: float,
    ambiguity_likelihood: bool,
    ambiguity_neutral_prob: float,
    min_fraction_ratio: float,
    truth: Sequence[str],
    candidate: Sequence[str],
) -> dict[str, str]:
    """Run one KIR2DL5 typing model and collect rank/support diagnostics."""
    from graphkir.typing_mulit_allele import AlleleTyping
    from graphkir2.typing.private_support import collect_variant_support  # type: ignore[import-not-found]

    if not copy_number:
        return {
            "selected": "",
            "selected_rank": "",
            "truth_rank": "",
            "candidate_rank": "",
            "truth_value_gap": "",
            "candidate_value_gap": "",
            "selected_support": "",
            "truth_support": "",
            "candidate_support": "",
            "reads_after_filter": "0",
            "fail_reason": "cn_zero_or_missing",
        }
    if not reads:
        return {
            "selected": "",
            "selected_rank": "",
            "truth_rank": "",
            "candidate_rank": "",
            "truth_value_gap": "",
            "candidate_value_gap": "",
            "selected_support": "",
            "truth_support": "",
            "candidate_support": "",
            "reads_after_filter": "0",
            "fail_reason": "no_reads_for_backbone",
        }
    if not variants:
        return {
            "selected": "",
            "selected_rank": "",
            "truth_rank": "",
            "candidate_rank": "",
            "truth_value_gap": "",
            "candidate_value_gap": "",
            "selected_support": "",
            "truth_support": "",
            "candidate_support": "",
            "reads_after_filter": "0",
            "fail_reason": "no_variants_for_backbone",
        }

    model = AlleleTyping(
        reads,
        variants,
        force_homo=False,
        top_n=top_n,
        variant_correction=True,
        exon_weight=exon_weight,
        ambiguity_likelihood=ambiguity_likelihood,
        ambiguity_neutral_prob=ambiguity_neutral_prob,
    )
    if not model.probs.shape[0]:
        return {
            "selected": "",
            "selected_rank": "",
            "truth_rank": "",
            "candidate_rank": "",
            "truth_value_gap": "",
            "candidate_value_gap": "",
            "selected_support": "",
            "truth_support": "",
            "candidate_support": "",
            "reads_after_filter": "0",
            "fail_reason": "no_informative_reads_after_filter",
        }

    result = model.typing(copy_number)
    selected = result.selectBest(min_fraction_ratio=min_fraction_ratio)
    selected_rank, _, _ = find_candidate_rank(result, selected)
    truth_rank, _, truth_gap = find_candidate_rank(result, truth)
    candidate_rank, _, candidate_gap = find_candidate_rank(result, candidate)
    positive, negative = collect_variant_support(model.reads)
    allele_variants = {
        allele: {
            str(variant.id)
            for variant in variants
            if allele in variant.allele
        }
        for allele in model.allele_to_id
    }
    return {
        "selected": allele_text(selected),
        "selected_rank": selected_rank,
        "truth_rank": truth_rank,
        "candidate_rank": candidate_rank,
        "truth_value_gap": truth_gap,
        "candidate_value_gap": candidate_gap,
        "selected_support": support_score_text(selected, allele_variants, positive, negative),
        "truth_support": support_score_text(truth, allele_variants, positive, negative),
        "candidate_support": support_score_text(candidate, allele_variants, positive, negative),
        "reads_after_filter": str(model.getReadsNum()),
        "fail_reason": "",
    }


def config_path_for_label(label: str) -> Path:
    """Return the benchmark config whose generated paths should be audited."""
    enhanced = Path("benchmarks/configs") / f"{label}-conditional-kir2ds3-enhancedgate.json"
    if enhanced.exists():
        return enhanced
    return Path("benchmarks/configs") / f"{label}.json"


def load_typing_context(
    label: str,
    sample_id: str,
) -> tuple[dict[str, Any], dict[str, int], Path]:
    """Load reads, variants, CN hints, and the config path for one sample."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir.hisat2 import loadReadsAndVariantsData, removeMultipleMapped
    from graphkir.kir_cn import loadCN
    from graphkir.kir_typing import (
        TypingWithPosNegAllele,
        groupReads,
        groupVariants,
        likelihoodAmbiguousMapped,
    )
    from graphkir2.benchmark import load_benchmark_preset  # type: ignore[import-not-found]
    from graphkir2.core.pipeline import GraphKir2Pipeline  # type: ignore[import-not-found]
    from graphkir2.io.manifest import load_sample_manifest  # type: ignore[import-not-found]

    config_path = config_path_for_label(label)
    preset = load_benchmark_preset(str(config_path))
    plan = GraphKir2Pipeline(preset.to_run_config()).build_typing_plan()
    sample_plan = next(
        sample
        for sample in plan.samples
        if extract_sample_id(sample.sample_id) == sample_id
    )
    manifest = load_sample_manifest(preset.input_csv)
    cn_hints = {sample.sample_id: sample.copy_number_hint for sample in manifest.samples}

    raw_reads_data = loadReadsAndVariantsData(sample_plan.variant_json)
    likelihood_reads_data = likelihoodAmbiguousMapped(deepcopy(raw_reads_data))
    discard_reads_data = removeMultipleMapped(deepcopy(raw_reads_data))

    likelihood_gene_reads = groupReads(likelihood_reads_data["reads"])
    gene_variants = groupVariants(raw_reads_data["variants"])
    dummy_model = TypingWithPosNegAllele.__new__(TypingWithPosNegAllele)
    dummy_model._gene_reads = likelihood_gene_reads
    dummy_model._gene_variants = gene_variants
    raw_cn = loadCN(cn_hints.get(sample_plan.sample_id, "") or sample_plan.cn_tsv)
    normalized_cn = dummy_model._normalizeGeneCopyNumbers(raw_cn)

    context = {
        "sample_plan": sample_plan,
        "raw_gene_reads": groupReads(raw_reads_data["reads"]),
        "likelihood_gene_reads": likelihood_gene_reads,
        "discard_gene_reads": groupReads(discard_reads_data["reads"]),
        "gene_variants": gene_variants,
        "raw_cn": raw_cn,
        "normalized_cn": normalized_cn,
        "select_min_fraction_ratio": sample_plan.select_min_fraction_ratio,
        "ambiguity_neutral_prob": sample_plan.ambiguity_neutral_prob,
    }
    return context, normalized_cn, config_path


def inspect_kir2dl5_sample(
    label: str,
    sample_id: str,
    summary_row: dict[str, str],
    truth: Sequence[str],
    candidate: Sequence[str],
) -> dict[str, str]:
    """Inspect read/variant/rank diagnostics for one remaining KIR2DL5 sample."""
    context, normalized_cn, config_path = load_typing_context(label, sample_id)
    raw_gene_reads: dict[str, list[Any]] = context["raw_gene_reads"]
    likelihood_gene_reads: dict[str, list[Any]] = context["likelihood_gene_reads"]
    discard_gene_reads: dict[str, list[Any]] = context["discard_gene_reads"]
    gene_variants: dict[str, list[Any]] = context["gene_variants"]
    raw_cn: dict[str, int] = context["raw_cn"]

    raw_reads = raw_gene_reads.get(KIR2DL5_BACKBONE, [])
    likelihood_reads = likelihood_gene_reads.get(KIR2DL5_BACKBONE, [])
    discard_reads = discard_gene_reads.get(KIR2DL5_BACKBONE, [])
    variants = gene_variants.get(KIR2DL5_BACKBONE, [])
    copy_number = normalized_cn.get(KIR2DL5_BACKBONE, 0)
    top_n = resolve_gene_top_n("KIR2DL5", summary_row)

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
    return {
        "config": str(config_path),
        "typing_backbone": KIR2DL5_BACKBONE,
        "raw_cn_kir2dl5a": str(raw_cn.get("KIR2DL5A*BACKBONE", 0)),
        "raw_cn_kir2dl5b": str(raw_cn.get("KIR2DL5B*BACKBONE", 0)),
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
        **{
            f"likelihood_{key}": value
            for key, value in likelihood_evidence.items()
        },
        **{
            f"discard_{key}": value
            for key, value in discard_evidence.items()
        },
    }


def remaining_key(row: dict[str, str]) -> tuple[str, str]:
    """Build a sample key from a remaining-error row."""
    return row["panel"], row["sample_id"]


def build_audit_rows(
    summary_rows: list[dict[str, str]],
    remaining_rows: list[dict[str, str]],
    candidate_method: str,
) -> list[dict[str, str]]:
    """Build KIR2DL5A/B audit rows for sample-level remaining failures."""
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
    for (label, sample_id), sample_errors in sorted(rows_by_sample.items()):
        summary_row = by_label_method[(label, candidate_method)]
        panel = calls[label]
        truth = merged_gene_alleles(panel["truth"].get(sample_id, []))
        candidate = merged_gene_alleles(panel[candidate_method].get(sample_id, []))
        discard = merged_gene_alleles(panel.get("discard", {}).get(sample_id, []))
        likelihood = merged_gene_alleles(
            panel.get("likelihood_top5000", {}).get(sample_id, [])
        )
        cause = classify_kir2dl5_cause(truth, candidate, discard, likelihood)
        evidence = inspect_kir2dl5_sample(label, sample_id, summary_row, truth, candidate)
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
                "truth_gene_counts": gene_copy_text(gene_copy_counts(truth)),
                "candidate_gene_counts": gene_copy_text(gene_copy_counts(candidate)),
                "discard_gene_counts": gene_copy_text(gene_copy_counts(discard)),
                "likelihood_gene_counts": gene_copy_text(gene_copy_counts(likelihood)),
                **evidence,
            }
        )
    return audit_rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write audit rows as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else [
        "panel",
        "sample_id",
        "root_cause",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def render_markdown(rows: list[dict[str, str]], output_tsv: str) -> str:
    """Render the KIR2DL5 audit Markdown report."""
    cause_counts = Counter(row["root_cause"] for row in rows)
    lines = [
        "# KIR2DL5 Remaining Failure Audit",
        "",
        "This report audits the remaining KIR2DL5A/B functional errors for the",
        "current functional-stress lead at the merged KIR2DL5 typing-backbone scope.",
        "",
        "Command:",
        "",
        "```bash",
        "python benchmarks/scripts/inspect_kir2dl5_remaining_failures.py",
        "```",
        "",
        f"Full TSV: `{output_tsv}`",
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
            "## Detail",
            "",
            "| panel | sample | root cause | truth | candidate | truth CN | candidate CN | normalized CN | likelihood selected | truth rank | truth gap | truth support | candidate support | discard selected | discard truth rank |",
            "|---|---:|---|---|---|---|---|---:|---|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["panel"],
                    row["sample_id"],
                    row["root_cause"],
                    row["truth_full"],
                    row["candidate_full"],
                    row["truth_gene_counts"],
                    row["candidate_gene_counts"],
                    row["normalized_cn"],
                    row["likelihood_selected"],
                    row["likelihood_truth_rank"] or "not in top-n",
                    row["likelihood_truth_value_gap"] or "",
                    row["likelihood_truth_support"] or "",
                    row["likelihood_candidate_support"] or "",
                    row["discard_selected"],
                    row["discard_truth_rank"] or "not in top-n",
                ]
            )
            + " |"
        )

    lines.extend(["", "## Interpretation", ""])
    if not rows:
        lines.extend(
            [
                "* No KIR2DL5A/B 3-digit or 5-digit remaining errors are present for the current candidate.",
                "* The KIR2DL5-only unsupported-overcall guard resolved the previous A/B placement and KIR2DL5A*001-vs-*012 overcall cases without leaving a KIR2DL5 functional miss in this stress sweep.",
                "* Keep the guard narrow: it should continue to require multiple unsupported selected-only variants before replacing a likelihood-selected KIR2DL5 call.",
            ]
        )
    else:
        lines.extend(
            [
                "* Remaining rows indicate KIR2DL5 still needs case-specific follow-up before promoting broader changes.",
                "* Prefer inspecting discriminating variants for these rows before changing KIR2DS3/KIR2DS5 guards.",
            ]
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Build the KIR2DL5 remaining-failure audit report."""
    args = build_parser().parse_args()
    summary_rows = read_sweep_summary(Path(args.summary_tsv))
    remaining_rows = read_remaining_rows(Path(args.remaining_tsv), KIR2DL5_GENES)
    audit_rows = build_audit_rows(summary_rows, remaining_rows, args.candidate_method)
    write_tsv(Path(args.output_tsv), audit_rows)
    output_md = Path(args.output_md)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(render_markdown(audit_rows, args.output_tsv), encoding="utf-8")
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.output_md}")


if __name__ == "__main__":
    main()
