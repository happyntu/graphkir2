"""Benchmark metadata for graphkir2."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass

from ..cn.interface import CopyNumberPlan
from .legacy import LegacyCommandPlan
from ..mapping.interface import MappingPlan
from ..typing.interface import TypingPlan


@dataclass(frozen=True)
class BenchmarkLabel:
    """Logical benchmark label for old-vs-new comparisons."""

    name: str


@dataclass(frozen=True)
class BenchmarkSummary:
    """Summary of the refactored pipeline planning surface for comparisons."""

    label: BenchmarkLabel
    manifest_path: str
    mapping: MappingPlan
    copy_number: CopyNumberPlan
    typing: TypingPlan

    def to_json(self) -> str:
        """Serialize the benchmark summary as formatted JSON."""
        return json.dumps(asdict(self), indent=2)

    def describe(self) -> str:
        """Render a compact textual benchmark summary."""
        return "\n".join(
            [
                "graphkir2 benchmark summary",
                f"label={self.label.name}",
                f"manifest={self.manifest_path}",
                f"mapping_samples={len(self.mapping.samples)}",
                f"multi_map_policy={self.mapping.multi_map_policy}",
                f"cn_mode={self.copy_number.mode}",
                f"typing_strategy={self.typing.strategy}",
                f"typing_exon_weight={self.typing.exon_weight}",
                f"typing_gene_exon_weights={self.typing.gene_exon_weights or '<none>'}",
                f"typing_private_support_genes={self.typing.private_support_genes or '<none>'}",
                f"typing_private_support_condition_alleles={self.typing.private_support_condition_alleles or '<none>'}",
                f"typing_private_support_cross_gene_ratio={self.typing.private_support_cross_gene_ratio}",
                f"typing_private_support_discard_fallback_genes={self.typing.private_support_discard_fallback_genes or '<none>'}",
                f"merged_cn={self.copy_number.cohort.merged_cn_tsv}",
                f"merged_alleles={self.typing.merged_allele_tsv}",
            ]
        )


@dataclass(frozen=True)
class ComparisonSummary:
    """High-level old-vs-new planning summary for benchmark scaffolding."""

    benchmark: BenchmarkSummary
    legacy_command: LegacyCommandPlan

    def to_json(self) -> str:
        """Serialize the comparison summary as formatted JSON."""
        return json.dumps(asdict(self), indent=2)

    def describe(self) -> str:
        """Render a compact textual old-vs-new planning summary."""
        return "\n".join(
            [
                "graphkir benchmark compare",
                f"label={self.benchmark.label.name}",
                f"manifest={self.benchmark.manifest_path}",
                f"legacy_command={self.legacy_command.render_shell()}",
                f"legacy_exact_match={self.legacy_command.exact_config_match}",
                f"mapping_samples={len(self.benchmark.mapping.samples)}",
                f"multi_map_policy={self.benchmark.mapping.multi_map_policy}",
                f"cn_mode={self.benchmark.copy_number.mode}",
                f"typing_strategy={self.benchmark.typing.strategy}",
                f"typing_exon_weight={self.benchmark.typing.exon_weight}",
                f"typing_gene_exon_weights={self.benchmark.typing.gene_exon_weights or '<none>'}",
                f"typing_private_support_genes={self.benchmark.typing.private_support_genes or '<none>'}",
                f"typing_private_support_condition_alleles={self.benchmark.typing.private_support_condition_alleles or '<none>'}",
                f"typing_private_support_cross_gene_ratio={self.benchmark.typing.private_support_cross_gene_ratio}",
                f"typing_private_support_discard_fallback_genes={self.benchmark.typing.private_support_discard_fallback_genes or '<none>'}",
                f"graphkir2_cn={self.benchmark.copy_number.cohort.merged_cn_tsv}",
                f"graphkir2_alleles={self.benchmark.typing.merged_allele_tsv}",
            ]
            + [f"legacy_note={note}" for note in self.legacy_command.notes]
        )
