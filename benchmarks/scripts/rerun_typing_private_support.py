"""Re-run typing with private-variant support reranking ablations."""

from __future__ import annotations

import argparse
import sys
from copy import deepcopy
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    """Build the private-support rerun parser."""
    parser = argparse.ArgumentParser(
        description="Ablate private-variant support reranking.",
    )
    parser.add_argument("--config", required=True, help="Benchmark preset JSON path.")
    parser.add_argument("--output-tsv", required=True, help="Prediction TSV path.")
    parser.add_argument("--top-n", type=int, default=5000, help="Typing top-n.")
    parser.add_argument(
        "--base-top-n",
        type=int,
        default=None,
        help=(
            "Override config lower top-n for genes outside private-support, "
            "fallback, or tie-break target sets."
        ),
    )
    parser.add_argument(
        "--gene-base-top-ns",
        default=None,
        help="Override config gene-specific non-target top-n values, e.g. `KIR2DL1:1000`.",
    )
    parser.add_argument(
        "--neutralize-cross-gene",
        action="store_true",
        help="Set cross-gene ambiguous read target weights to zero before typing.",
    )
    parser.add_argument(
        "--neutralize-cross-gene-groups",
        default=None,
        help="Optional comma-separated groups such as `KIR2DS3/KIR2DS5`.",
    )
    parser.add_argument(
        "--private-support-genes",
        default=None,
        help="Optional comma-separated pure genes to rerank, e.g. `KIR2DS3`.",
    )
    parser.add_argument(
        "--private-support-lambda",
        type=float,
        default=None,
        help="Weight for private-variant support reranking.",
    )
    parser.add_argument(
        "--private-support-window",
        type=float,
        default=None,
        help="Only rerank candidates within this log-likelihood window.",
    )
    parser.add_argument(
        "--private-support-condition-alleles",
        default=None,
        help="Optional comma-separated allele prefixes required before rescue.",
    )
    parser.add_argument(
        "--private-support-cross-gene-ratio",
        type=float,
        default=None,
        help="Optional minimum cross-gene ratio for selected-private positive support.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-genes",
        default=None,
        help="Optional comma-separated genes allowed to use discard-style fallback.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-residual-alleles",
        default=None,
        help="Optional allele prefixes that trigger fallback if still present after rescue.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-introduced-alleles",
        default=None,
        help="Optional allele prefixes that trigger fallback if introduced by rescue.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-introduced-max-ratio",
        type=float,
        default=None,
        help="Optional maximum cross-gene ratio for introduced-allele fallback.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-max-score",
        type=float,
        default=None,
        help="Optional maximum base private-support score allowed for fallback.",
    )
    parser.add_argument(
        "--private-support-discard-fallback-residual-min-ratio",
        type=float,
        default=None,
        help="Optional minimum cross-gene ratio for residual-allele fallback.",
    )
    parser.add_argument(
        "--functional-discard-fallback-genes",
        default=None,
        help="Optional comma-separated genes allowed to use discard evidence for weak functional likelihood calls.",
    )
    parser.add_argument(
        "--functional-discard-fallback-resolution",
        type=int,
        default=None,
        help="Optional allele field length for functional discard fallback comparison.",
    )
    parser.add_argument(
        "--functional-discard-fallback-min-score-delta",
        type=float,
        default=None,
        help="Optional private-support score advantage required for functional discard fallback.",
    )
    parser.add_argument(
        "--functional-discard-fallback-max-score",
        type=float,
        default=None,
        help="Optional maximum selected-call private-support score allowed for functional fallback.",
    )
    parser.add_argument(
        "--functional-discard-fallback-promoted-alleles",
        default=None,
        help="Optional allele prefixes whose likelihood-only promotion should be guarded.",
    )
    parser.add_argument(
        "--functional-discard-fallback-protected-alleles",
        default=None,
        help="Optional discard allele prefixes to protect from configured promotions.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-genes",
        default="",
        help="Optional comma-separated genes allowed to guard unsupported candidate-only overcalls.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-window",
        type=float,
        default=25.0,
        help="Maximum likelihood gap for unsupported-overcall guard alternatives.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-min-unsupported-delta",
        type=int,
        default=2,
        help="Minimum unsupported-variant count improvement required by the overcall guard.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-min-net-delta",
        type=float,
        default=20.0,
        help="Minimum negative-minus-positive support improvement required by the overcall guard.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-negative-threshold",
        type=float,
        default=5.0,
        help="Negative support threshold for unsupported candidate-only variants.",
    )
    parser.add_argument(
        "--unsupported-overcall-guard-max-positive",
        type=float,
        default=1.0,
        help="Maximum positive support allowed for unsupported candidate-only variants.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-genes",
        default="",
        help="Optional genes allowed to guard selected allele-prefix unsupported overcalls.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-alleles",
        default="",
        help="Allele prefixes required in the selected call before targeted unsupported guard runs.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-window",
        type=float,
        default=25.0,
        help="Maximum likelihood gap for targeted unsupported-overcall guard alternatives.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-min-unsupported-delta",
        type=int,
        default=1,
        help="Minimum targeted unsupported-variant improvement required by the guard.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-min-net-delta",
        type=float,
        default=20.0,
        help="Minimum targeted negative-minus-positive support improvement required by the guard.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-preserve-non-target-resolution",
        type=int,
        default=0,
        help="Optional allele field length that targeted alternatives must preserve outside the target prefixes.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-negative-threshold",
        type=float,
        default=5.0,
        help="Negative support threshold for targeted unsupported candidate-only variants.",
    )
    parser.add_argument(
        "--targeted-unsupported-overcall-guard-max-positive",
        type=float,
        default=1.0,
        help="Maximum positive support allowed for targeted unsupported candidate-only variants.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-genes",
        default="",
        help="Optional genes allowed to run a second rank-wide targeted unsupported-overcall guard.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-alleles",
        default="",
        help="Allele prefixes required in the selected call before the rank-wide guard runs.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-window",
        type=float,
        default=400.0,
        help="Maximum likelihood gap for rank-wide unsupported-overcall alternatives.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-min-unsupported-delta",
        type=int,
        default=1,
        help="Minimum targeted unsupported-variant improvement required by the rank-wide guard.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-min-net-delta",
        type=float,
        default=20.0,
        help="Minimum targeted negative-minus-positive support improvement required by the rank-wide guard.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-max-selected-support",
        type=float,
        default=0.0,
        help="Maximum selected private-support score allowed before the rank-wide guard runs.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-preserve-non-target-resolution",
        type=int,
        default=0,
        help="Optional allele field length that rank-wide alternatives must preserve outside the target prefixes.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-negative-threshold",
        type=float,
        default=5.0,
        help="Negative support threshold for rank-wide unsupported candidate-only variants.",
    )
    parser.add_argument(
        "--rankwide-unsupported-overcall-guard-max-positive",
        type=float,
        default=1.0,
        help="Maximum positive support allowed for rank-wide unsupported candidate-only variants.",
    )
    return parser


def main() -> None:
    """Run private-support reranking."""
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
    from graphkir.typing_mulit_allele import AlleleTyping
    from graphkir2.benchmark import load_benchmark_preset
    from graphkir2.core.pipeline import GraphKir2Pipeline
    from graphkir2.io.manifest import load_sample_manifest
    from graphkir2.typing.private_support import (
        choose_targeted_top_n,
        apply_functional_promotion_guard,
        collect_variant_support,
        neutralize_cross_gene_reads,
        parse_gene_groups,
        parse_gene_set,
        parse_gene_top_n_spec,
        parse_name_set,
        private_positive_cross_gene_ratio,
        private_support_score,
        pure_gene,
        select_against_unsupported_candidate_only_variants,
        select_with_highest_suffix_tie_break,
        select_with_private_support,
        should_apply_conditional_private_support,
        should_use_functional_discard_fallback,
        should_use_discard_fallback,
    )

    args = build_parser().parse_args()
    preset = load_benchmark_preset(args.config)
    run_config = preset.to_run_config()
    plan = GraphKir2Pipeline(run_config).build_typing_plan()
    manifest = load_sample_manifest(preset.input_csv)
    cn_hints = {sample.sample_id: sample.copy_number_hint for sample in manifest.samples}
    neutralize_group_spec = (
        args.neutralize_cross_gene_groups
        if args.neutralize_cross_gene_groups is not None
        else run_config.typing.cross_gene_neutralization_groups
    )
    private_support_gene_spec = (
        args.private_support_genes
        if args.private_support_genes is not None
        else run_config.typing.private_support_genes
    )
    private_support_lambda = (
        args.private_support_lambda
        if args.private_support_lambda is not None
        else run_config.typing.private_support_lambda
    )
    private_support_window = (
        args.private_support_window
        if args.private_support_window is not None
        else run_config.typing.private_support_window
    )
    private_support_condition_spec = (
        args.private_support_condition_alleles
        if args.private_support_condition_alleles is not None
        else run_config.typing.private_support_condition_alleles
    )
    private_support_cross_gene_ratio = (
        args.private_support_cross_gene_ratio
        if args.private_support_cross_gene_ratio is not None
        else run_config.typing.private_support_cross_gene_ratio
    )
    fallback_gene_spec = (
        args.private_support_discard_fallback_genes
        if args.private_support_discard_fallback_genes is not None
        else run_config.typing.private_support_discard_fallback_genes
    )
    fallback_residual_spec = (
        args.private_support_discard_fallback_residual_alleles
        if args.private_support_discard_fallback_residual_alleles is not None
        else run_config.typing.private_support_discard_fallback_residual_alleles
    )
    fallback_introduced_spec = (
        args.private_support_discard_fallback_introduced_alleles
        if args.private_support_discard_fallback_introduced_alleles is not None
        else run_config.typing.private_support_discard_fallback_introduced_alleles
    )
    fallback_introduced_max_ratio = (
        args.private_support_discard_fallback_introduced_max_ratio
        if args.private_support_discard_fallback_introduced_max_ratio is not None
        else run_config.typing.private_support_discard_fallback_introduced_max_ratio
    )
    fallback_max_score = (
        args.private_support_discard_fallback_max_score
        if args.private_support_discard_fallback_max_score is not None
        else run_config.typing.private_support_discard_fallback_max_score
    )
    fallback_residual_min_ratio = (
        args.private_support_discard_fallback_residual_min_ratio
        if args.private_support_discard_fallback_residual_min_ratio is not None
        else run_config.typing.private_support_discard_fallback_residual_min_ratio
    )
    functional_fallback_gene_spec = (
        args.functional_discard_fallback_genes
        if args.functional_discard_fallback_genes is not None
        else run_config.typing.functional_discard_fallback_genes
    )
    functional_fallback_resolution = (
        args.functional_discard_fallback_resolution
        if args.functional_discard_fallback_resolution is not None
        else run_config.typing.functional_discard_fallback_resolution
    )
    functional_fallback_min_score_delta = (
        args.functional_discard_fallback_min_score_delta
        if args.functional_discard_fallback_min_score_delta is not None
        else run_config.typing.functional_discard_fallback_min_score_delta
    )
    functional_fallback_max_score = (
        args.functional_discard_fallback_max_score
        if args.functional_discard_fallback_max_score is not None
        else run_config.typing.functional_discard_fallback_max_score
    )
    functional_fallback_promoted_spec = (
        args.functional_discard_fallback_promoted_alleles
        if args.functional_discard_fallback_promoted_alleles is not None
        else run_config.typing.functional_discard_fallback_promoted_alleles
    )
    functional_fallback_protected_spec = (
        args.functional_discard_fallback_protected_alleles
        if args.functional_discard_fallback_protected_alleles is not None
        else run_config.typing.functional_discard_fallback_protected_alleles
    )
    base_top_n = (
        args.base_top_n
        if args.base_top_n is not None
        else run_config.typing.base_top_n
    )
    gene_base_top_n_spec = (
        args.gene_base_top_ns
        if args.gene_base_top_ns is not None
        else run_config.typing.gene_base_top_ns
    )
    neutralize_groups = parse_gene_groups(neutralize_group_spec)
    private_support_genes = parse_gene_set(private_support_gene_spec)
    private_support_condition_alleles = parse_name_set(private_support_condition_spec)
    fallback_genes = parse_gene_set(fallback_gene_spec)
    fallback_residual_alleles = parse_name_set(fallback_residual_spec)
    fallback_introduced_alleles = parse_name_set(fallback_introduced_spec)
    functional_fallback_genes = parse_gene_set(functional_fallback_gene_spec)
    functional_fallback_promoted_alleles = parse_name_set(functional_fallback_promoted_spec)
    functional_fallback_protected_alleles = parse_name_set(functional_fallback_protected_spec)
    unsupported_overcall_guard_genes = parse_gene_set(args.unsupported_overcall_guard_genes)
    targeted_unsupported_overcall_guard_genes = parse_gene_set(
        args.targeted_unsupported_overcall_guard_genes
    )
    targeted_unsupported_overcall_guard_alleles = parse_name_set(
        args.targeted_unsupported_overcall_guard_alleles
    )
    rankwide_unsupported_overcall_guard_genes = parse_gene_set(
        args.rankwide_unsupported_overcall_guard_genes
    )
    rankwide_unsupported_overcall_guard_alleles = parse_name_set(
        args.rankwide_unsupported_overcall_guard_alleles
    )
    gene_base_top_ns = parse_gene_top_n_spec(gene_base_top_n_spec)
    has_conditional_gate = bool(
        private_support_condition_alleles or private_support_cross_gene_ratio > 0.0
    )
    highest_suffix_tie_break_genes = parse_gene_set(
        run_config.typing.highest_suffix_tie_break_genes
    )
    high_top_n_genes = (
        private_support_genes
        | fallback_genes
        | highest_suffix_tie_break_genes
    )

    rows: list[dict[str, str]] = []
    for sample in plan.samples:
        raw_reads_data = loadReadsAndVariantsData(sample.variant_json)
        reads_data = likelihoodAmbiguousMapped(deepcopy(raw_reads_data))
        if (args.neutralize_cross_gene or neutralize_groups) and not has_conditional_gate:
            neutralize_cross_gene_reads(
                reads_data["reads"],
                neutralize_groups,
                target_genes=private_support_genes,
            )
        gene_reads = groupReads(reads_data["reads"])
        gene_variants = groupVariants(reads_data["variants"])
        dummy_model = TypingWithPosNegAllele.__new__(TypingWithPosNegAllele)
        dummy_model._gene_reads = gene_reads
        dummy_model._gene_variants = gene_variants
        cn = dummy_model._normalizeGeneCopyNumbers(
            loadCN(cn_hints.get(sample.sample_id, "") or sample.cn_tsv)
        )

        alleles: list[str] = []
        warnings: list[str] = []
        for gene, copy_number in cn.items():
            if not copy_number:
                continue
            gene_name = pure_gene(gene)
            gene_top_n = choose_targeted_top_n(
                gene_name,
                args.top_n,
                base_top_n if base_top_n > 0 else None,
                high_top_n_genes,
                gene_base_top_ns,
            )
            model = AlleleTyping(
                gene_reads[gene],
                gene_variants[gene],
                top_n=gene_top_n,
                variant_correction=True,
                exon_weight=2.0,
                ambiguity_likelihood=True,
                ambiguity_neutral_prob=0.999,
            )
            result = model.typing(copy_number)
            positive, negative = collect_variant_support(model.reads)
            allele_variants = {
                allele: {
                    str(variant.id)
                    for variant in gene_variants[gene]
                    if allele in variant.allele
                }
                for allele in model.allele_to_id
            }
            support_lambda = private_support_lambda
            selected_without_rescue = result.selectBest(
                min_fraction_ratio=sample.select_min_fraction_ratio
            )
            conditional_cross_gene_ratio = 0.0
            base_private_support = private_support_score(
                selected_without_rescue,
                allele_variants,
                positive,
                negative,
            )
            if private_support_genes and gene_name not in private_support_genes:
                support_lambda = 0.0
            if support_lambda and (
                private_support_condition_alleles
                or private_support_cross_gene_ratio > 0.0
            ):
                if private_support_cross_gene_ratio > 0.0:
                    conditional_cross_gene_ratio = private_positive_cross_gene_ratio(
                        reads_data["reads"],
                        selected_without_rescue,
                        allele_variants,
                        neutralize_groups,
                    )
                should_rescue = should_apply_conditional_private_support(
                    selected_without_rescue,
                    reads_data["reads"],
                    allele_variants,
                    private_support_condition_alleles,
                    neutralize_groups,
                    private_support_cross_gene_ratio,
                )
                if should_rescue and neutralize_groups:
                    rescue_reads_data = deepcopy(reads_data["reads"])
                    neutralize_cross_gene_reads(
                        rescue_reads_data,
                        neutralize_groups,
                        target_genes=private_support_genes,
                    )
                    rescue_gene_reads = groupReads(rescue_reads_data)
                    model = AlleleTyping(
                        rescue_gene_reads[gene],
                        gene_variants[gene],
                        top_n=gene_top_n,
                        variant_correction=True,
                        exon_weight=2.0,
                        ambiguity_likelihood=True,
                        ambiguity_neutral_prob=0.999,
                    )
                    result = model.typing(copy_number)
                    positive, negative = collect_variant_support(model.reads)
                    allele_variants = {
                        allele: {
                            str(variant.id)
                            for variant in gene_variants[gene]
                            if allele in variant.allele
                        }
                        for allele in model.allele_to_id
                    }
                elif not should_rescue:
                    support_lambda = 0.0
            selected = select_with_private_support(
                result,
                allele_variants,
                positive,
                negative,
                support_lambda,
                private_support_window,
                sample.select_min_fraction_ratio,
            )
            if (
                gene_name in fallback_genes
                and should_use_discard_fallback(
                    selected,
                    selected_without_rescue,
                    fallback_residual_alleles,
                    fallback_introduced_alleles,
                    conditional_cross_gene_ratio,
                    fallback_introduced_max_ratio,
                    base_private_support,
                    fallback_max_score,
                    fallback_residual_min_ratio,
                )
            ):
                discard_reads_data = removeMultipleMapped(deepcopy(raw_reads_data))
                discard_gene_reads = groupReads(discard_reads_data["reads"])
                if gene in discard_gene_reads:
                    fallback_model = AlleleTyping(
                        discard_gene_reads[gene],
                        gene_variants[gene],
                        top_n=gene_top_n,
                        variant_correction=True,
                        exon_weight=1.0,
                        ambiguity_likelihood=False,
                    )
                    fallback_result = fallback_model.typing(copy_number)
                    selected = fallback_result.selectBest(
                        min_fraction_ratio=sample.select_min_fraction_ratio
                    )
                    model = fallback_model
            if gene_name in functional_fallback_genes:
                selected_support = private_support_score(
                    selected,
                    allele_variants,
                    positive,
                    negative,
                )
                has_promotion_guard = bool(
                    functional_fallback_promoted_alleles
                    and functional_fallback_protected_alleles
                )
                can_use_full_functional_fallback = (
                    selected_support <= functional_fallback_max_score
                )
                if can_use_full_functional_fallback or has_promotion_guard:
                    discard_reads_data = removeMultipleMapped(deepcopy(raw_reads_data))
                    discard_gene_reads = groupReads(discard_reads_data["reads"])
                    if gene in discard_gene_reads:
                        fallback_model = AlleleTyping(
                            discard_gene_reads[gene],
                            gene_variants[gene],
                            top_n=gene_top_n,
                            variant_correction=True,
                            exon_weight=1.0,
                            ambiguity_likelihood=False,
                        )
                        fallback_result = fallback_model.typing(copy_number)
                        fallback_selected = fallback_result.selectBest(
                            min_fraction_ratio=sample.select_min_fraction_ratio
                        )
                        fallback_positive, fallback_negative = collect_variant_support(
                            fallback_model.reads
                        )
                        fallback_support = private_support_score(
                            fallback_selected,
                            allele_variants,
                            fallback_positive,
                            fallback_negative,
                        )
                        if can_use_full_functional_fallback and should_use_functional_discard_fallback(
                            selected,
                            fallback_selected,
                            selected_support,
                            fallback_support,
                            functional_fallback_resolution,
                            functional_fallback_min_score_delta,
                        ):
                            selected = fallback_selected
                            model = fallback_model
                        else:
                            guarded_selected = apply_functional_promotion_guard(
                                selected,
                                fallback_selected,
                                functional_fallback_promoted_alleles,
                                functional_fallback_protected_alleles,
                                functional_fallback_resolution,
                            )
                            if guarded_selected != selected:
                                selected = guarded_selected
            if gene_name in unsupported_overcall_guard_genes:
                selected = select_against_unsupported_candidate_only_variants(
                    result,
                    selected,
                    allele_variants,
                    positive,
                    negative,
                    sample.select_min_fraction_ratio,
                    args.unsupported_overcall_guard_window,
                    args.unsupported_overcall_guard_min_unsupported_delta,
                    args.unsupported_overcall_guard_min_net_delta,
                    args.unsupported_overcall_guard_negative_threshold,
                    args.unsupported_overcall_guard_max_positive,
                )
            if gene_name in targeted_unsupported_overcall_guard_genes:
                selected = select_against_unsupported_candidate_only_variants(
                    result,
                    selected,
                    allele_variants,
                    positive,
                    negative,
                    sample.select_min_fraction_ratio,
                    args.targeted_unsupported_overcall_guard_window,
                    args.targeted_unsupported_overcall_guard_min_unsupported_delta,
                    args.targeted_unsupported_overcall_guard_min_net_delta,
                    args.targeted_unsupported_overcall_guard_negative_threshold,
                    args.targeted_unsupported_overcall_guard_max_positive,
                    targeted_unsupported_overcall_guard_alleles,
                    args.targeted_unsupported_overcall_guard_preserve_non_target_resolution,
                )
            rankwide_selected_support = 0.0
            if gene_name in rankwide_unsupported_overcall_guard_genes:
                rankwide_selected_support = private_support_score(
                    selected,
                    allele_variants,
                    positive,
                    negative,
                )
            if (
                gene_name in rankwide_unsupported_overcall_guard_genes
                and rankwide_selected_support
                <= args.rankwide_unsupported_overcall_guard_max_selected_support
            ):
                selected = select_against_unsupported_candidate_only_variants(
                    result,
                    selected,
                    allele_variants,
                    positive,
                    negative,
                    sample.select_min_fraction_ratio,
                    args.rankwide_unsupported_overcall_guard_window,
                    args.rankwide_unsupported_overcall_guard_min_unsupported_delta,
                    args.rankwide_unsupported_overcall_guard_min_net_delta,
                    args.rankwide_unsupported_overcall_guard_negative_threshold,
                    args.rankwide_unsupported_overcall_guard_max_positive,
                    rankwide_unsupported_overcall_guard_alleles,
                    args.rankwide_unsupported_overcall_guard_preserve_non_target_resolution,
                )
            if gene_name in highest_suffix_tie_break_genes:
                selected = select_with_highest_suffix_tie_break(
                    result,
                    selected,
                    sample.select_min_fraction_ratio,
                )
            alleles.extend(
                allele if allele != "fail" else f"{gene_name}*"
                for allele in selected
            )
            if model.getReadsNum() < 100:
                warnings.append(gene)
        rows.append(
            {
                "name": sample.output_prefix,
                "alleles": "_".join(alleles),
                "warnings": "_".join(warnings),
            }
        )

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
