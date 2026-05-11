"""CLI entrypoint for the refactored Graph-KIR implementation."""

from __future__ import annotations

import argparse

from .benchmark.runner import BenchmarkLabel
from .config import (
    CopyNumberConfig,
    EngineConfig,
    GraphKir2RuntimeConfig,
    IndexConfig,
    InputConfig,
    MappingConfig,
    ReferenceConfig,
    TypingConfig,
)
from .core.pipeline import GraphKir2Pipeline, GraphKir2RunConfig
from .io.manifest import load_sample_manifest


def build_parser() -> argparse.ArgumentParser:
    """Build the graphkir2 CLI parser."""
    parser = argparse.ArgumentParser(
        prog="graphkir2",
        description="Refactored Graph-KIR pipeline scaffold.",
    )
    parser.add_argument(
        "--print-plan",
        action="store_true",
        help="Print the refactor pipeline plan instead of running work.",
    )
    parser.add_argument(
        "--benchmark-label",
        default="dev",
        help="Label used to distinguish benchmark runs.",
    )
    parser.add_argument(
        "--print-layout",
        action="store_true",
        help="Print the current graphkir2 module layout.",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print the current graphkir2 configuration snapshot.",
    )
    parser.add_argument(
        "--print-mapping-plan",
        action="store_true",
        help="Build and print the mapping plan derived from config + manifest.",
    )
    parser.add_argument(
        "--print-cn-plan",
        action="store_true",
        help="Build and print the copy-number plan derived from mapping + config.",
    )
    parser.add_argument(
        "--print-typing-plan",
        action="store_true",
        help="Build and print the typing plan derived from mapping + CN + config.",
    )
    parser.add_argument(
        "--print-benchmark-summary",
        action="store_true",
        help="Build and print an end-to-end benchmark planning summary.",
    )
    parser.add_argument(
        "--print-benchmark-json",
        action="store_true",
        help="Build and print an end-to-end benchmark summary as JSON.",
    )
    parser.add_argument(
        "--print-compare-summary",
        action="store_true",
        help="Print an old-vs-new comparison summary with the legacy command.",
    )
    parser.add_argument(
        "--print-compare-json",
        action="store_true",
        help="Print an old-vs-new comparison summary as JSON.",
    )
    parser.add_argument(
        "--input-csv",
        help="Path to a Graph-KIR style sample manifest CSV.",
    )
    parser.add_argument(
        "--print-manifest",
        action="store_true",
        help="Parse and print the sample manifest specified by --input-csv.",
    )
    parser.add_argument("--thread", type=int, default=2, help="Number of threads.")
    parser.add_argument(
        "--memory-gb", type=int, default=7, help="Memory budget hint in GB."
    )
    parser.add_argument(
        "--engine",
        default="local",
        choices=["local", "docker", "podman", "singularity"],
        help="Execution backend.",
    )
    parser.add_argument(
        "--ref-genome",
        default="hg19",
        choices=["hg19", "hg38"],
        help="Reference genome selection.",
    )
    parser.add_argument(
        "--index-folder", default="index", help="Index folder for Graph-KIR assets."
    )
    parser.add_argument(
        "--ipd-version", default="2100", help="IPD-KIR version label."
    )
    parser.add_argument(
        "--msa-type", default="ab_2dl1s1", help="MSA layout choice."
    )
    parser.add_argument(
        "--msa-no-exon-only-allele",
        action="store_true",
        help="Disable exon-only alleles when building the index model.",
    )
    parser.add_argument(
        "--multi-map-policy",
        default="discard",
        choices=["discard", "best-only", "weighted", "margin", "likelihood"],
        help="How graphkir2 should treat multi-mapped reads in the mapping stage.",
    )
    parser.add_argument(
        "--allele-strategy",
        default="full",
        choices=["full", "exonfirst"],
        help="Allele typing strategy.",
    )
    parser.add_argument(
        "--allele-exon-weight",
        type=float,
        default=1.0,
        help="Extra likelihood weight applied to exon variants during typing.",
    )
    parser.add_argument(
        "--allele-gene-exon-weights",
        default="",
        help="Optional `GENE:WEIGHT,...` overrides for exon weighting.",
    )
    parser.add_argument(
        "--allele-margin-gap",
        type=float,
        default=5.0,
        help="Top1-top2 score gap threshold below which near-tie alignments share weight.",
    )
    parser.add_argument(
        "--allele-margin-scale",
        type=float,
        default=2.0,
        help="Softmax score scale used when graphkir2 keeps near-tie alignments in margin mode.",
    )
    parser.add_argument(
        "--allele-ambiguity-neutral-prob",
        type=float,
        default=0.999,
        help="Neutral allele probability used for non-target read mass in likelihood mode.",
    )
    parser.add_argument(
        "--allele-select-min-fraction-ratio",
        type=float,
        default=0.5,
        help="Minimum selected allele fraction as a ratio of expected 1/CN fraction.",
    )
    parser.add_argument(
        "--allele-cross-gene-neutralization-groups",
        default="",
        help="Comma-separated slash groups to neutralize in cross-gene ambiguous reads, e.g. `KIR2DS3/KIR2DS5`.",
    )
    parser.add_argument(
        "--allele-private-support-genes",
        default="",
        help="Comma-separated pure genes to rerank with private-variant support.",
    )
    parser.add_argument(
        "--allele-private-support-lambda",
        type=float,
        default=0.0,
        help="Weight for private-variant support reranking.",
    )
    parser.add_argument(
        "--allele-private-support-window",
        type=float,
        default=0.0,
        help="Only rerank candidates within this log-likelihood window.",
    )
    parser.add_argument(
        "--allele-private-support-condition-alleles",
        default="",
        help="Comma-separated allele prefixes required before private-support rescue, e.g. `KIR2DS3*00201`.",
    )
    parser.add_argument(
        "--allele-private-support-cross-gene-ratio",
        type=float,
        default=0.0,
        help="Minimum cross-gene share of selected-private positive support before rescue.",
    )
    parser.add_argument(
        "--allele-private-support-discard-fallback-genes",
        default="",
        help="Comma-separated genes allowed to fall back to discard-style typing after rescue.",
    )
    parser.add_argument(
        "--allele-private-support-discard-fallback-residual-alleles",
        default="",
        help="Allele prefixes that trigger discard fallback if still present after rescue.",
    )
    parser.add_argument(
        "--allele-private-support-discard-fallback-introduced-alleles",
        default="",
        help="Allele prefixes that trigger fallback if introduced by rescue.",
    )
    parser.add_argument(
        "--allele-private-support-discard-fallback-introduced-max-ratio",
        type=float,
        default=0.0,
        help="Maximum cross-gene ratio for introduced-allele discard fallback.",
    )
    parser.add_argument(
        "--allele-highest-suffix-tie-break-genes",
        default="",
        help="Comma-separated genes where exact likelihood ties keep the highest 7-digit suffix within the same 5-digit call.",
    )
    parser.add_argument(
        "--cn-diploid-gene",
        default="",
        help="Diploid gene hint used by CN estimation.",
    )
    parser.add_argument(
        "--cn-cohort",
        action="store_true",
        help="Enable cohort-based CN estimation.",
    )
    parser.add_argument(
        "--cn-3dl3-not-diploid",
        action="store_true",
        help="Do not assume KIR3DL3 is diploid.",
    )
    parser.add_argument(
        "--step-skip-extraction",
        action="store_true",
        help="Skip WGS extraction in the future pipeline.",
    )
    parser.add_argument(
        "--output-cohort-name",
        default="cohort",
        help="Cohort output prefix used when planning merged outputs.",
    )
    return parser


def entrypoint() -> None:
    """Run the graphkir2 scaffold CLI."""
    args = build_parser().parse_args()
    config = GraphKir2RunConfig(
        benchmark=BenchmarkLabel(name=args.benchmark_label),
        runtime=GraphKir2RuntimeConfig(
            threads=args.thread,
            memory_gb=args.memory_gb,
        ),
        engine=EngineConfig(name=args.engine),
        reference=ReferenceConfig(genome=args.ref_genome),
        index=IndexConfig(
            index_folder=args.index_folder,
            ipd_version=args.ipd_version,
            msa_type=args.msa_type,
            use_exon_only_alleles=not args.msa_no_exon_only_allele,
        ),
        inputs=InputConfig(
            manifest_path=args.input_csv or "",
            skip_extraction=args.step_skip_extraction,
            output_folder="",
            output_cohort_name=args.output_cohort_name,
        ),
        mapping=MappingConfig(multi_map_policy=args.multi_map_policy),
        copy_number=CopyNumberConfig(
            diploid_gene=args.cn_diploid_gene,
            cohort_mode=args.cn_cohort,
            assume_3dl3_diploid=not args.cn_3dl3_not_diploid,
        ),
        typing=TypingConfig(
            strategy=args.allele_strategy,
            exon_weight=args.allele_exon_weight,
            gene_exon_weights=args.allele_gene_exon_weights,
            margin_gap=args.allele_margin_gap,
            margin_scale=args.allele_margin_scale,
            ambiguity_neutral_prob=args.allele_ambiguity_neutral_prob,
            select_min_fraction_ratio=args.allele_select_min_fraction_ratio,
            cross_gene_neutralization_groups=args.allele_cross_gene_neutralization_groups,
            private_support_genes=args.allele_private_support_genes,
            private_support_lambda=args.allele_private_support_lambda,
            private_support_window=args.allele_private_support_window,
            private_support_condition_alleles=args.allele_private_support_condition_alleles,
            private_support_cross_gene_ratio=args.allele_private_support_cross_gene_ratio,
            private_support_discard_fallback_genes=args.allele_private_support_discard_fallback_genes,
            private_support_discard_fallback_residual_alleles=args.allele_private_support_discard_fallback_residual_alleles,
            private_support_discard_fallback_introduced_alleles=args.allele_private_support_discard_fallback_introduced_alleles,
            private_support_discard_fallback_introduced_max_ratio=args.allele_private_support_discard_fallback_introduced_max_ratio,
            highest_suffix_tie_break_genes=args.allele_highest_suffix_tie_break_genes,
        ),
        print_plan=args.print_plan,
        print_layout=args.print_layout,
        print_config=args.print_config,
    )
    pipeline = GraphKir2Pipeline(config)
    if args.print_manifest:
        if not args.input_csv:
            raise SystemExit("--print-manifest requires --input-csv")
        manifest = load_sample_manifest(args.input_csv)
        print(manifest.describe())
        return
    if args.print_plan:
        print(pipeline.describe())
        return
    if args.print_layout:
        print(pipeline.describe_layout())
        return
    if args.print_config:
        print(pipeline.describe_config())
        return
    if args.print_mapping_plan:
        print(pipeline.describe_mapping_plan())
        return
    if args.print_cn_plan:
        print(pipeline.describe_copy_number_plan())
        return
    if args.print_typing_plan:
        print(pipeline.describe_typing_plan())
        return
    if args.print_benchmark_summary:
        print(pipeline.describe_benchmark_summary())
        return
    if args.print_benchmark_json:
        print(pipeline.benchmark_summary_json())
        return
    if args.print_compare_summary:
        print(pipeline.describe_comparison_summary())
        return
    if args.print_compare_json:
        print(pipeline.comparison_summary_json())
        return
    print(
        "graphkir2 is scaffolded for refactor work. "
        "Implement stages under src/graphkir2/core and compare them via benchmarks/."
    )
