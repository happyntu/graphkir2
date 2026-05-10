"""Shared domain models for the refactored implementation."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SampleInput:
    """Input metadata for one sequencing sample."""

    sample_id: str
    read1: str
    read2: str
    copy_number_hint: str = ""
    output_prefix: str = ""


@dataclass(frozen=True)
class GeneCopyNumber:
    """Copy-number prediction for one gene."""

    gene: str
    copy_number: int
    confidence: float = 0.0


@dataclass(frozen=True)
class AlleleCall:
    """Allele call for one gene or merged locus."""

    gene: str
    alleles: tuple[str, ...]
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class SampleResult:
    """Top-level output bundle for one sample."""

    sample_id: str
    copy_numbers: tuple[GeneCopyNumber, ...] = ()
    allele_calls: tuple[AlleleCall, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class PipelinePlan:
    """Named refactor stages tracked by the new implementation."""

    stages: tuple[str, ...] = field(
        default_factory=lambda: (
            "input-loading",
            "mapping",
            "copy-number",
            "typing",
            "reporting",
        )
    )
