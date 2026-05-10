"""Re-run typing with private-variant support reranking ablations."""

from __future__ import annotations

import argparse
import sys
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
    return parser


def main() -> None:
    """Run private-support reranking."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir.hisat2 import loadReadsAndVariantsData
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
        collect_variant_support,
        neutralize_cross_gene_reads,
        parse_gene_groups,
        parse_gene_set,
        pure_gene,
        select_with_private_support,
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
    neutralize_groups = parse_gene_groups(neutralize_group_spec)
    private_support_genes = parse_gene_set(private_support_gene_spec)

    rows: list[dict[str, str]] = []
    for sample in plan.samples:
        reads_data = likelihoodAmbiguousMapped(loadReadsAndVariantsData(sample.variant_json))
        if args.neutralize_cross_gene or neutralize_groups:
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
            model = AlleleTyping(
                gene_reads[gene],
                gene_variants[gene],
                top_n=args.top_n,
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
            gene_name = pure_gene(gene)
            support_lambda = private_support_lambda
            if private_support_genes and gene_name not in private_support_genes:
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
