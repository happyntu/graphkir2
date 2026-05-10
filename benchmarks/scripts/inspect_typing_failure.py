"""Inspect per-gene typing failure signals from existing variant JSON files."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the typing failure inspection parser."""
    parser = argparse.ArgumentParser(
        description="Inspect read/variant/candidate evidence for selected genes.",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Benchmark preset JSON path.",
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Comma-separated pure gene names, e.g. KIR2DL5A,KIR2DL5B.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path to write the inspection TSV.",
    )
    return parser


def _pure_gene(backbone: str) -> str:
    """Convert KIRX*BACKBONE to KIRX."""
    return backbone.split("*", 1)[0]


def _typing_backbone_alias(gene_or_backbone: str) -> str:
    """Return the merged graph backbone used for genes represented jointly."""
    if gene_or_backbone.startswith("KIR2DL5A") or gene_or_backbone.startswith("KIR2DL5B"):
        return "KIR2DL5*BACKBONE"
    return gene_or_backbone


def _has_read_evidence(read: object) -> bool:
    """Return whether a PairRead carries any positive or negative variant evidence."""
    return bool(read.lpv + read.lnv + read.rpv + read.rnv)  # type: ignore[attr-defined]


def _collect_variant_counts(reads: Sequence[object]) -> tuple[int, int]:
    """Count unique positive and negative variant IDs observed in reads."""
    positive: set[str] = set()
    negative: set[str] = set()
    for read in reads:
        positive.update(read.lpv + read.rpv)  # type: ignore[attr-defined]
        negative.update(read.lnv + read.rnv)  # type: ignore[attr-defined]
    return len(positive), len(negative)


def _sample_id_from_name(name: str) -> str:
    """Extract sample id from generated synthetic sample name."""
    sample_name = Path(name).name
    parts = sample_name.split(".")
    if len(parts) >= 2 and parts[1].isdigit():
        return parts[1]
    return sample_name


def main() -> None:
    """Inspect selected genes for each sample in a benchmark preset."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir.hisat2 import loadReadsAndVariantsData
    from graphkir.kir_cn import loadCN
    from graphkir.kir_typing import (
        groupReads,
        groupVariants,
        likelihoodAmbiguousMapped,
    )
    from graphkir.typing_mulit_allele import AlleleTyping
    from graphkir2.benchmark import load_benchmark_preset
    from graphkir2.core.pipeline import GraphKir2Pipeline
    from graphkir2.io.manifest import load_sample_manifest

    args = build_parser().parse_args()
    genes = tuple(gene.strip() for gene in args.genes.split(",") if gene.strip())
    preset = load_benchmark_preset(args.config)
    pipeline = GraphKir2Pipeline(preset.to_run_config())
    typing_plan = pipeline.build_typing_plan()
    manifest = load_sample_manifest(preset.input_csv)
    cn_hints = {sample.sample_id: sample.copy_number_hint for sample in manifest.samples}

    rows: list[dict[str, str]] = []
    for sample_plan in typing_plan.samples:
        sample_id = _sample_id_from_name(sample_plan.sample_id)
        cn_path = cn_hints.get(sample_plan.sample_id, "") or sample_plan.cn_tsv
        cn = loadCN(cn_path)
        reads_data = loadReadsAndVariantsData(sample_plan.variant_json)
        likelihood_reads_data = likelihoodAmbiguousMapped(reads_data)
        raw_reads_by_gene = groupReads(reads_data["reads"])
        likelihood_reads_by_gene = groupReads(likelihood_reads_data["reads"])
        variants_by_gene = groupVariants(reads_data["variants"])

        pure_read_counts = Counter(_pure_gene(read.backbone) for read in reads_data["reads"])
        pure_likelihood_read_counts = Counter(
            _pure_gene(read.backbone) for read in likelihood_reads_data["reads"]
        )
        for gene in genes:
            expected_backbone = gene + "*BACKBONE"
            alias_backbone = _typing_backbone_alias(expected_backbone)
            matching_backbones = [
                backbone
                for backbone in sorted(set(raw_reads_by_gene) | set(variants_by_gene) | set(cn))
                if _pure_gene(backbone) == gene
            ]
            if (
                alias_backbone != expected_backbone
                and (
                    alias_backbone in raw_reads_by_gene
                    or alias_backbone in variants_by_gene
                )
            ):
                matching_backbones.append(alias_backbone)
            if not matching_backbones:
                matching_backbones = [expected_backbone]

            for backbone in matching_backbones:
                normalized_cn = sum(
                    value
                    for cn_backbone, value in cn.items()
                    if _typing_backbone_alias(cn_backbone) == backbone
                )
                raw_reads = raw_reads_by_gene.get(backbone, [])
                likelihood_reads = likelihood_reads_by_gene.get(backbone, [])
                variants = variants_by_gene.get(backbone, [])
                raw_informative = [read for read in raw_reads if _has_read_evidence(read)]
                likelihood_informative = [
                    read for read in likelihood_reads if _has_read_evidence(read)
                ]
                positive_variants, negative_variants = _collect_variant_counts(likelihood_reads)
                candidate_alleles = sorted(
                    {
                        allele
                        for variant in variants
                        for allele in variant.allele
                    }
                )
                after_error_correction = ""
                after_empty_filter = ""
                top_result = ""
                fail_reason = ""
                if variants:
                    model = AlleleTyping(
                        likelihood_reads,
                        variants,
                        top_n=50,
                        variant_correction=True,
                        exon_weight=getattr(preset, "allele_exon_weight", 1.0),
                        ambiguity_likelihood=True,
                        ambiguity_neutral_prob=getattr(
                            preset,
                            "allele_ambiguity_neutral_prob",
                            0.999,
                        ),
                        ambiguity_target_weight_power=getattr(
                            preset,
                            "allele_ambiguity_target_weight_power",
                            1.0,
                        ),
                    )
                    after_error_correction = str(len(model.reads))
                    after_empty_filter = str(len(model.probs))
                    if model.probs.shape[0] and normalized_cn:
                        result = model.typing(normalized_cn)
                        top_result = "_".join(result.selectBest())
                    elif not model.probs.shape[0]:
                        fail_reason = "no_informative_reads_after_filter"
                    elif not normalized_cn:
                        fail_reason = "cn_zero_or_missing"
                elif not raw_reads:
                    fail_reason = "no_reads_for_backbone"
                else:
                    fail_reason = "no_variants_for_backbone"

                rows.append(
                    {
                        "sample_id": sample_id,
                        "gene": gene,
                        "backbone": backbone,
                        "cn": str(cn.get(backbone, "")),
                        "normalized_cn": str(normalized_cn or ""),
                        "alias_backbone": alias_backbone if alias_backbone != expected_backbone else "",
                        "raw_reads": str(len(raw_reads)),
                        "raw_informative_reads": str(len(raw_informative)),
                        "likelihood_reads": str(len(likelihood_reads)),
                        "likelihood_informative_reads": str(len(likelihood_informative)),
                        "pure_gene_raw_reads": str(pure_read_counts[gene]),
                        "pure_gene_likelihood_reads": str(pure_likelihood_read_counts[gene]),
                        "variants": str(len(variants)),
                        "exon_variants": str(sum(1 for variant in variants if variant.in_exon)),
                        "candidate_alleles": str(len(candidate_alleles)),
                        "positive_variant_ids": str(positive_variants),
                        "negative_variant_ids": str(negative_variants),
                        "after_error_correction_reads": after_error_correction,
                        "after_empty_filter_reads": after_empty_filter,
                        "top_result": top_result,
                        "fail_reason": fail_reason,
                    }
                )

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "gene",
        "backbone",
        "cn",
        "normalized_cn",
        "alias_backbone",
        "raw_reads",
        "raw_informative_reads",
        "likelihood_reads",
        "likelihood_informative_reads",
        "pure_gene_raw_reads",
        "pure_gene_likelihood_reads",
        "variants",
        "exon_variants",
        "candidate_alleles",
        "positive_variant_ids",
        "negative_variant_ids",
        "after_error_correction_reads",
        "after_empty_filter_reads",
        "top_result",
        "fail_reason",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
