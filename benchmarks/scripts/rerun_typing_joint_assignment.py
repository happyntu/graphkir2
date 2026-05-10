"""Re-run typing with a two-pass cross-gene ambiguous-read assignment ablation."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    """Build the joint-assignment rerun parser."""
    parser = argparse.ArgumentParser(
        description="Ablate two-pass cross-gene ambiguous-read assignment.",
    )
    parser.add_argument("--config", required=True, help="Benchmark preset JSON path.")
    parser.add_argument("--output-tsv", required=True, help="Prediction TSV path.")
    parser.add_argument(
        "--top-n",
        type=int,
        default=5000,
        help="Typing top-n parameter for both passes.",
    )
    parser.add_argument(
        "--assignment-margin",
        type=float,
        default=0.0,
        help="Minimum log10 score margin required for hard assignment.",
    )
    parser.add_argument(
        "--keep-uncertain",
        action="store_true",
        help="Keep original weights when cross-gene assignment margin is too small.",
    )
    return parser


def _read_name(read: object) -> str:
    """Extract a read-pair name from a PairRead-like object."""
    if read.l_sam:  # type: ignore[attr-defined]
        return read.l_sam.split("\t", 1)[0]  # type: ignore[attr-defined]
    if read.r_sam:  # type: ignore[attr-defined]
        return read.r_sam.split("\t", 1)[0]  # type: ignore[attr-defined]
    return ""


def _gene(backbone: str) -> str:
    """Extract the pure gene from a backbone name."""
    return backbone.split("*", 1)[0]


def main() -> None:
    """Run two-pass typing with joint assignment."""
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

    args = build_parser().parse_args()
    preset = load_benchmark_preset(args.config)
    plan = GraphKir2Pipeline(preset.to_run_config()).build_typing_plan()
    manifest = load_sample_manifest(preset.input_csv)
    cn_hints = {sample.sample_id: sample.copy_number_hint for sample in manifest.samples}

    rows: list[dict[str, str]] = []
    for sample in plan.samples:
        reads_data = likelihoodAmbiguousMapped(loadReadsAndVariantsData(sample.variant_json))
        gene_reads = groupReads(reads_data["reads"])
        gene_variants = groupVariants(reads_data["variants"])
        dummy_model = TypingWithPosNegAllele.__new__(TypingWithPosNegAllele)
        dummy_model._gene_reads = gene_reads
        dummy_model._gene_variants = gene_variants
        cn = dummy_model._normalizeGeneCopyNumbers(
            loadCN(cn_hints.get(sample.sample_id, "") or sample.cn_tsv)
        )

        first_models: dict[str, AlleleTyping] = {}
        first_call_ids: dict[str, list[int]] = {}
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
            called = result.selectBest(min_fraction_ratio=0.7)
            first_models[gene] = model
            first_call_ids[gene] = [
                model.allele_to_id[allele]
                for allele in called
                if allele in model.allele_to_id
            ]

        reads_by_name: dict[str, list[Any]] = defaultdict(list)
        for read in cast(list[Any], reads_data["reads"]):
            reads_by_name[_read_name(read)].append(read)

        for read_group in reads_by_name.values():
            genes = {_gene(read.backbone) for read in read_group}  # type: ignore[attr-defined]
            if len(genes) <= 1:
                continue
            scored_reads: list[tuple[float, Any]] = []
            for read in read_group:
                backbone = read.backbone  # type: ignore[attr-defined]
                assignment_model = first_models.get(backbone)
                ids = first_call_ids.get(backbone, [])
                if assignment_model is None or not ids:
                    continue
                prob = assignment_model.read2AlleleProb(read)
                if prob is None:
                    continue
                allele_score = max(float(prob[index]) for index in ids)
                score = float(np.log10(max(allele_score, 1e-300)))
                score += float(np.log10(max(read.weight, 1e-6)))  # type: ignore[attr-defined]
                scored_reads.append((score, read))
            if not scored_reads:
                continue
            scored_reads.sort(key=lambda item: item[0], reverse=True)
            best_score, best_read = scored_reads[0]
            second_score = scored_reads[1][0] if len(scored_reads) > 1 else -np.inf
            if best_score - second_score < args.assignment_margin:
                if args.keep_uncertain:
                    continue
                for _, read in scored_reads:
                    read.weight = 0.0  # type: ignore[attr-defined]
                    read.ambiguous_weight = 1.0  # type: ignore[attr-defined]
                continue
            for _, read in scored_reads:
                read.weight = 0.0  # type: ignore[attr-defined]
                read.ambiguous_weight = 1.0  # type: ignore[attr-defined]
            best_read.weight = 1.0  # type: ignore[attr-defined]
            best_read.ambiguous_weight = 0.0  # type: ignore[attr-defined]

        assigned_gene_reads = groupReads(reads_data["reads"])
        alleles: list[str] = []
        warnings: list[str] = []
        for gene, copy_number in cn.items():
            if not copy_number:
                continue
            model = AlleleTyping(
                assigned_gene_reads[gene],
                gene_variants[gene],
                top_n=args.top_n,
                variant_correction=True,
                exon_weight=2.0,
                ambiguity_likelihood=True,
                ambiguity_neutral_prob=0.999,
            )
            result = model.typing(copy_number)
            pure_gene = gene.split("*", 1)[0]
            alleles.extend(
                allele if allele != "fail" else f"{pure_gene}*"
                for allele in result.selectBest(min_fraction_ratio=0.7)
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
