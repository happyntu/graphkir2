"""Benchmark preset loading helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from ..config import (
    CopyNumberConfig,
    EngineConfig,
    GraphKir2RuntimeConfig,
    IndexConfig,
    InputConfig,
    MappingConfig,
    ReferenceConfig,
    TypingConfig,
)
from ..core.pipeline import GraphKir2RunConfig
from .runner import BenchmarkLabel


@dataclass(frozen=True)
class BenchmarkPreset:
    """Serializable benchmark preset for reproducible compare runs."""

    benchmark_label: str = "dev"
    input_csv: str = ""
    output_folder: str = ""
    output_cohort_name: str = "cohort"
    thread: int = 2
    memory_gb: int = 7
    engine: str = "local"
    ref_genome: str = "hg19"
    index_folder: str = "index"
    ipd_version: str = "2100"
    msa_type: str = "ab_2dl1s1"
    use_exon_only_alleles: bool = True
    multi_map_policy: str = "discard"
    allele_strategy: str = "full"
    allele_exon_weight: float = 1.0
    allele_gene_exon_weights: str = ""
    allele_margin_gap: float = 5.0
    allele_margin_scale: float = 2.0
    allele_ambiguity_neutral_prob: float = 0.999
    allele_select_min_fraction_ratio: float = 0.5
    allele_base_top_n: int = 0
    allele_gene_base_top_ns: str = ""
    allele_cross_gene_neutralization_groups: str = ""
    allele_private_support_genes: str = ""
    allele_private_support_lambda: float = 0.0
    allele_private_support_window: float = 0.0
    allele_private_support_condition_alleles: str = ""
    allele_private_support_cross_gene_ratio: float = 0.0
    allele_private_support_discard_fallback_genes: str = ""
    allele_private_support_discard_fallback_residual_alleles: str = ""
    allele_private_support_discard_fallback_introduced_alleles: str = ""
    allele_private_support_discard_fallback_introduced_max_ratio: float = 0.0
    allele_private_support_discard_fallback_max_score: float = 0.0
    allele_private_support_discard_fallback_residual_min_ratio: float = 0.0
    allele_highest_suffix_tie_break_genes: str = ""
    allele_truth_tsv: str = ""
    cn_truth_tsv: str = ""
    cn_diploid_gene: str = ""
    cn_cohort: bool = False
    assume_3dl3_diploid: bool = True
    step_skip_extraction: bool = False

    def to_run_config(self) -> GraphKir2RunConfig:
        """Convert the preset into a graphkir2 pipeline config."""
        return GraphKir2RunConfig(
            benchmark=BenchmarkLabel(name=self.benchmark_label),
            runtime=GraphKir2RuntimeConfig(
                threads=self.thread,
                memory_gb=self.memory_gb,
            ),
            engine=EngineConfig(name=self.engine),
            reference=ReferenceConfig(genome=self.ref_genome),
            index=IndexConfig(
                index_folder=self.index_folder,
                ipd_version=self.ipd_version,
                msa_type=self.msa_type,
                use_exon_only_alleles=self.use_exon_only_alleles,
            ),
            inputs=InputConfig(
                manifest_path=self.input_csv,
                skip_extraction=self.step_skip_extraction,
                output_folder=self.output_folder,
                output_cohort_name=self.output_cohort_name,
            ),
            mapping=MappingConfig(multi_map_policy=self.multi_map_policy),
            copy_number=CopyNumberConfig(
                diploid_gene=self.cn_diploid_gene,
                cohort_mode=self.cn_cohort,
                assume_3dl3_diploid=self.assume_3dl3_diploid,
            ),
            typing=TypingConfig(
                strategy=self.allele_strategy,
                exon_weight=self.allele_exon_weight,
                gene_exon_weights=self.allele_gene_exon_weights,
                margin_gap=self.allele_margin_gap,
                margin_scale=self.allele_margin_scale,
                ambiguity_neutral_prob=self.allele_ambiguity_neutral_prob,
                select_min_fraction_ratio=self.allele_select_min_fraction_ratio,
                base_top_n=self.allele_base_top_n,
                gene_base_top_ns=self.allele_gene_base_top_ns,
                cross_gene_neutralization_groups=self.allele_cross_gene_neutralization_groups,
                private_support_genes=self.allele_private_support_genes,
                private_support_lambda=self.allele_private_support_lambda,
                private_support_window=self.allele_private_support_window,
                private_support_condition_alleles=self.allele_private_support_condition_alleles,
                private_support_cross_gene_ratio=self.allele_private_support_cross_gene_ratio,
                private_support_discard_fallback_genes=self.allele_private_support_discard_fallback_genes,
                private_support_discard_fallback_residual_alleles=self.allele_private_support_discard_fallback_residual_alleles,
                private_support_discard_fallback_introduced_alleles=self.allele_private_support_discard_fallback_introduced_alleles,
                private_support_discard_fallback_introduced_max_ratio=self.allele_private_support_discard_fallback_introduced_max_ratio,
                private_support_discard_fallback_max_score=self.allele_private_support_discard_fallback_max_score,
                private_support_discard_fallback_residual_min_ratio=self.allele_private_support_discard_fallback_residual_min_ratio,
                highest_suffix_tie_break_genes=self.allele_highest_suffix_tie_break_genes,
            ),
        )


def load_benchmark_preset(path: str) -> BenchmarkPreset:
    """Load a benchmark preset from JSON."""
    preset_path = Path(path)
    data = json.loads(preset_path.read_text(encoding="utf-8"))
    return BenchmarkPreset(**data)
