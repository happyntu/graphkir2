"""Typing-stage interfaces for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass

from ..cn.interface import CopyNumberPlan
from ..config.runtime import TypingConfig
from ..mapping.interface import MappingPlan, sanitize_path_component


def get_common_prefix(left: str, right: str) -> str:
    """Return the common character prefix of two artifact paths."""
    size = min(len(left), len(right))
    for index in range(size):
        if left[index] != right[index]:
            return left[:index]
    return left[:size]


def build_typing_suffix(variant_prefix: str, cn_tsv: str, strategy: str) -> str:
    """Reproduce the legacy allele-typing suffix pattern."""
    common = get_common_prefix(variant_prefix, cn_tsv)
    suffix = ".cn" + sanitize_path_component(cn_tsv[len(common) :]) + "."
    if strategy == "exonfirst":
        return suffix + "exonfirst_1"
    return suffix + strategy


@dataclass(frozen=True)
class SampleTypingPlan:
    """Planned typing artifacts for one sample."""

    sample_id: str
    variant_json: str
    cn_tsv: str
    output_prefix: str
    allele_tsv: str
    possible_tsv: str
    strategy: str
    exon_weight: float
    gene_exon_weights: str
    margin_gap: float
    margin_scale: float
    ambiguity_neutral_prob: float
    select_min_fraction_ratio: float
    cross_gene_neutralization_groups: str
    private_support_genes: str
    private_support_lambda: float
    private_support_window: float


@dataclass(frozen=True)
class TypingPlan:
    """Top-level allele-typing plan."""

    strategy: str
    exon_weight: float
    gene_exon_weights: str
    margin_gap: float
    margin_scale: float
    ambiguity_neutral_prob: float
    select_min_fraction_ratio: float
    cross_gene_neutralization_groups: str
    private_support_genes: str
    private_support_lambda: float
    private_support_window: float
    merged_allele_tsv: str
    samples: tuple[SampleTypingPlan, ...]

    def describe(self) -> str:
        """Render a readable typing-plan summary."""
        lines = [
            "graphkir2 typing plan",
            f"strategy={self.strategy}",
            f"exon_weight={self.exon_weight}",
            f"gene_exon_weights={self.gene_exon_weights or '<none>'}",
            f"margin_gap={self.margin_gap}",
            f"margin_scale={self.margin_scale}",
            f"ambiguity_neutral_prob={self.ambiguity_neutral_prob}",
            f"select_min_fraction_ratio={self.select_min_fraction_ratio}",
            f"cross_gene_neutralization_groups={self.cross_gene_neutralization_groups or '<none>'}",
            f"private_support_genes={self.private_support_genes or '<none>'}",
            f"private_support_lambda={self.private_support_lambda}",
            f"private_support_window={self.private_support_window}",
            f"merged_allele_tsv={self.merged_allele_tsv}",
            f"samples={len(self.samples)}",
        ]
        for sample in self.samples:
            lines.extend(
                [
                    f"- sample={sample.sample_id}",
                    f"  variant_json={sample.variant_json}",
                    f"  cn_tsv={sample.cn_tsv}",
                    f"  output_prefix={sample.output_prefix}",
                    f"  allele_tsv={sample.allele_tsv}",
                    f"  possible_tsv={sample.possible_tsv}",
                    f"  strategy={sample.strategy}",
                    f"  exon_weight={sample.exon_weight}",
                    f"  gene_exon_weights={sample.gene_exon_weights or '<none>'}",
                    f"  margin_gap={sample.margin_gap}",
                    f"  margin_scale={sample.margin_scale}",
                    f"  ambiguity_neutral_prob={sample.ambiguity_neutral_prob}",
                    f"  select_min_fraction_ratio={sample.select_min_fraction_ratio}",
                    f"  cross_gene_neutralization_groups={sample.cross_gene_neutralization_groups or '<none>'}",
                    f"  private_support_genes={sample.private_support_genes or '<none>'}",
                    f"  private_support_lambda={sample.private_support_lambda}",
                    f"  private_support_window={sample.private_support_window}",
                ]
            )
        return "\n".join(lines)


@dataclass(frozen=True)
class AlleleTyper:
    """Marker/planner object for the allele typing stage."""

    stage_name: str

    def plan(
        self,
        mapping: MappingPlan,
        copy_number: CopyNumberPlan,
        config: TypingConfig,
        cohort_name: str,
    ) -> TypingPlan:
        """Build a non-executing typing plan from mapping and CN outputs."""
        if len(mapping.samples) != len(copy_number.samples):
            raise ValueError("Mapping and CN plans must contain the same number of samples")

        sample_plans: list[SampleTypingPlan] = []
        for mapping_sample, cn_sample in zip(mapping.samples, copy_number.samples):
            strategy = config.strategy
            suffix = build_typing_suffix(
                mapping_sample.variant_prefix,
                cn_sample.cn_tsv,
                strategy,
            )
            output_prefix = mapping_sample.variant_prefix + suffix
            sample_plans.append(
                SampleTypingPlan(
                    sample_id=mapping_sample.sample_id,
                    variant_json=mapping_sample.variant_prefix + ".json",
                    cn_tsv=cn_sample.cn_tsv,
                    output_prefix=output_prefix,
                    allele_tsv=output_prefix + ".tsv",
                    possible_tsv=output_prefix + ".possible.tsv",
                    strategy=strategy,
                    exon_weight=config.exon_weight,
                    gene_exon_weights=config.gene_exon_weights,
                    margin_gap=config.margin_gap,
                    margin_scale=config.margin_scale,
                    ambiguity_neutral_prob=config.ambiguity_neutral_prob,
                    select_min_fraction_ratio=config.select_min_fraction_ratio,
                    cross_gene_neutralization_groups=config.cross_gene_neutralization_groups,
                    private_support_genes=config.private_support_genes,
                    private_support_lambda=config.private_support_lambda,
                    private_support_window=config.private_support_window,
                )
            )

        return TypingPlan(
            strategy=config.strategy,
            exon_weight=config.exon_weight,
            gene_exon_weights=config.gene_exon_weights,
            margin_gap=config.margin_gap,
            margin_scale=config.margin_scale,
            ambiguity_neutral_prob=config.ambiguity_neutral_prob,
            select_min_fraction_ratio=config.select_min_fraction_ratio,
            cross_gene_neutralization_groups=config.cross_gene_neutralization_groups,
            private_support_genes=config.private_support_genes,
            private_support_lambda=config.private_support_lambda,
            private_support_window=config.private_support_window,
            merged_allele_tsv=cohort_name + ".allele.tsv",
            samples=tuple(sample_plans),
        )
