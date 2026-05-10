"""Runtime configuration primitives for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GraphKir2RuntimeConfig:
    """Resource and execution defaults for graphkir2."""

    threads: int = 2
    memory_gb: int = 7


@dataclass(frozen=True)
class EngineConfig:
    """Execution backend configuration."""

    name: str = "local"


@dataclass(frozen=True)
class ReferenceConfig:
    """Reference-genome selection for WGS extraction and related tasks."""

    genome: str = "hg19"


@dataclass(frozen=True)
class IndexConfig:
    """Index and database configuration for graphkir2."""

    index_folder: str = "index"
    ipd_version: str = "2100"
    msa_type: str = "ab_2dl1s1"
    use_exon_only_alleles: bool = True


@dataclass(frozen=True)
class InputConfig:
    """Input-side configuration."""

    manifest_path: str = ""
    skip_extraction: bool = False
    output_folder: str = ""
    output_cohort_name: str = ""


@dataclass(frozen=True)
class MappingConfig:
    """Mapping-stage configuration."""

    multi_map_policy: str = "discard"


@dataclass(frozen=True)
class TypingConfig:
    """Allele-typing configuration."""

    strategy: str = "full"
    exon_weight: float = 1.0
    gene_exon_weights: str = ""
    margin_gap: float = 5.0
    margin_scale: float = 2.0
    ambiguity_neutral_prob: float = 0.999
    select_min_fraction_ratio: float = 0.5
    cross_gene_neutralization_groups: str = ""
    private_support_genes: str = ""
    private_support_lambda: float = 0.0
    private_support_window: float = 0.0
    highest_suffix_tie_break_genes: str = ""


@dataclass(frozen=True)
class CopyNumberConfig:
    """Copy-number estimation configuration."""

    diploid_gene: str = ""
    cohort_mode: bool = False
    assume_3dl3_diploid: bool = True
