"""High-level pipeline scaffold for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass

from ..benchmark.legacy import build_legacy_graphkir_command
from ..benchmark.runner import BenchmarkLabel, BenchmarkSummary, ComparisonSummary
from ..cn.interface import CopyNumberEstimator, CopyNumberPlan
from ..config.runtime import (
    CopyNumberConfig,
    EngineConfig,
    GraphKir2RuntimeConfig,
    IndexConfig,
    InputConfig,
    MappingConfig,
    ReferenceConfig,
    TypingConfig,
)
from ..io.layout import ModuleLayout
from ..io.manifest import SampleManifest, load_sample_manifest
from ..mapping.interface import MappingPlan, ReadMapper
from ..typing.interface import AlleleTyper, TypingPlan


@dataclass(frozen=True)
class GraphKir2RunConfig:
    """Runtime configuration for the refactored pipeline scaffold."""

    benchmark: BenchmarkLabel
    runtime: GraphKir2RuntimeConfig
    engine: EngineConfig
    reference: ReferenceConfig
    index: IndexConfig
    inputs: InputConfig
    mapping: MappingConfig
    copy_number: CopyNumberConfig
    typing: TypingConfig
    print_plan: bool = False
    print_layout: bool = False
    print_config: bool = False


class GraphKir2Pipeline:
    """Minimal pipeline object used to anchor the refactor work."""

    def __init__(self, config: GraphKir2RunConfig):
        self.config = config
        self.mapper = ReadMapper(stage_name="mapping")
        self.copy_number = CopyNumberEstimator(stage_name="copy-number")
        self.typer = AlleleTyper(stage_name="typing")
        self.layout = ModuleLayout.default()

    def describe(self) -> str:
        """Describe the intended refactor phases."""
        return "\n".join(
            [
                "graphkir2 scaffold",
                f"benchmark_label={self.config.benchmark.name}",
                f"reference={self.config.reference.genome}",
                f"index_folder={self.config.index.index_folder}",
                "phases:",
                "1. extract stable domain models from legacy graphkir",
                "2. redesign IO and pipeline stage boundaries",
                "3. reimplement core mapping/CN/typing flow with tests",
                "4. benchmark graphkir2 against legacy graphkir",
                "stages:",
                f"- {self.mapper.stage_name}",
                f"- {self.copy_number.stage_name}",
                f"- {self.typer.stage_name}",
            ]
        )

    def describe_layout(self) -> str:
        """Describe the graphkir2 module layout."""
        lines = [
            "graphkir2 module layout",
            f"threads={self.config.runtime.threads}",
            f"memory_gb={self.config.runtime.memory_gb}",
        ]
        lines.extend(f"- {entry}" for entry in self.layout.entries)
        return "\n".join(lines)

    def describe_config(self) -> str:
        """Describe the current graphkir2 configuration snapshot."""
        return "\n".join(
            [
                "graphkir2 config",
                f"benchmark_label={self.config.benchmark.name}",
                f"threads={self.config.runtime.threads}",
                f"memory_gb={self.config.runtime.memory_gb}",
                f"engine={self.config.engine.name}",
                f"reference={self.config.reference.genome}",
                f"index_folder={self.config.index.index_folder}",
                f"ipd_version={self.config.index.ipd_version}",
                f"msa_type={self.config.index.msa_type}",
                f"use_exon_only_alleles={self.config.index.use_exon_only_alleles}",
                f"manifest_path={self.config.inputs.manifest_path or '<unset>'}",
                f"skip_extraction={self.config.inputs.skip_extraction}",
                f"multi_map_policy={self.config.mapping.multi_map_policy}",
                f"typing_strategy={self.config.typing.strategy}",
                f"typing_exon_weight={self.config.typing.exon_weight}",
                f"typing_gene_exon_weights={self.config.typing.gene_exon_weights or '<none>'}",
                f"typing_margin_gap={self.config.typing.margin_gap}",
                f"typing_margin_scale={self.config.typing.margin_scale}",
                f"typing_ambiguity_neutral_prob={self.config.typing.ambiguity_neutral_prob}",
                f"typing_select_min_fraction_ratio={self.config.typing.select_min_fraction_ratio}",
                f"typing_base_top_n={self.config.typing.base_top_n or '<none>'}",
                f"typing_gene_base_top_ns={self.config.typing.gene_base_top_ns or '<none>'}",
                f"typing_cross_gene_neutralization_groups={self.config.typing.cross_gene_neutralization_groups or '<none>'}",
                f"typing_private_support_genes={self.config.typing.private_support_genes or '<none>'}",
                f"typing_private_support_lambda={self.config.typing.private_support_lambda}",
                f"typing_private_support_window={self.config.typing.private_support_window}",
                f"typing_private_support_condition_alleles={self.config.typing.private_support_condition_alleles or '<none>'}",
                f"typing_private_support_cross_gene_ratio={self.config.typing.private_support_cross_gene_ratio}",
                f"typing_private_support_discard_fallback_genes={self.config.typing.private_support_discard_fallback_genes or '<none>'}",
                f"typing_private_support_discard_fallback_residual_alleles={self.config.typing.private_support_discard_fallback_residual_alleles or '<none>'}",
                f"typing_private_support_discard_fallback_introduced_alleles={self.config.typing.private_support_discard_fallback_introduced_alleles or '<none>'}",
                f"typing_private_support_discard_fallback_introduced_max_ratio={self.config.typing.private_support_discard_fallback_introduced_max_ratio}",
                f"typing_private_support_discard_fallback_max_score={self.config.typing.private_support_discard_fallback_max_score}",
                f"typing_private_support_discard_fallback_residual_min_ratio={self.config.typing.private_support_discard_fallback_residual_min_ratio}",
                f"typing_highest_suffix_tie_break_genes={self.config.typing.highest_suffix_tie_break_genes or '<none>'}",
                f"diploid_gene={self.config.copy_number.diploid_gene or '<none>'}",
                f"cohort_mode={self.config.copy_number.cohort_mode}",
                f"assume_3dl3_diploid={self.config.copy_number.assume_3dl3_diploid}",
            ]
        )

    def load_manifest(self) -> SampleManifest:
        """Load the configured sample manifest."""
        if not self.config.inputs.manifest_path:
            raise ValueError("No manifest configured. Set --input-csv first.")
        return load_sample_manifest(self.config.inputs.manifest_path)

    def build_mapping_plan(self) -> MappingPlan:
        """Build the mapping plan derived from the current configuration."""
        manifest = self.load_manifest()
        return self.mapper.plan(
            samples=manifest.samples,
            runtime=self.config.runtime,
            index=self.config.index,
            reference=self.config.reference,
            inputs=self.config.inputs,
            mapping=self.config.mapping,
        )

    def build_copy_number_plan(self, mapping: MappingPlan | None = None) -> CopyNumberPlan:
        """Build the copy-number plan derived from mapping + CN configuration."""
        mapping_plan = mapping or self.build_mapping_plan()
        cohort_name = self.config.inputs.output_cohort_name or "cohort"
        return self.copy_number.plan(
            mapping=mapping_plan,
            config=self.config.copy_number,
            cohort_name=cohort_name,
        )

    def build_typing_plan(
        self,
        mapping: MappingPlan | None = None,
        copy_number: CopyNumberPlan | None = None,
    ) -> TypingPlan:
        """Build the typing plan derived from mapping + CN + typing config."""
        mapping_plan = mapping or self.build_mapping_plan()
        copy_number_plan = copy_number or self.build_copy_number_plan(mapping_plan)
        cohort_name = self.config.inputs.output_cohort_name or "cohort"
        return self.typer.plan(
            mapping=mapping_plan,
            copy_number=copy_number_plan,
            config=self.config.typing,
            cohort_name=cohort_name,
        )

    def build_benchmark_summary(self) -> BenchmarkSummary:
        """Build an end-to-end benchmark planning summary."""
        mapping = self.build_mapping_plan()
        copy_number = self.build_copy_number_plan(mapping)
        typing = self.build_typing_plan(mapping, copy_number)
        return BenchmarkSummary(
            label=self.config.benchmark,
            manifest_path=self.config.inputs.manifest_path or "<unset>",
            mapping=mapping,
            copy_number=copy_number,
            typing=typing,
        )

    def build_comparison_summary(self) -> ComparisonSummary:
        """Build a benchmark summary paired with the equivalent legacy command."""
        return ComparisonSummary(
            benchmark=self.build_benchmark_summary(),
            legacy_command=build_legacy_graphkir_command(self.config),
        )

    def describe_mapping_plan(self) -> str:
        """Describe the mapping plan derived from the current configuration."""
        return self.build_mapping_plan().describe()

    def describe_copy_number_plan(self) -> str:
        """Describe the CN plan derived from mapping + CN configuration."""
        return self.build_copy_number_plan().describe()

    def describe_typing_plan(self) -> str:
        """Describe the typing plan derived from mapping + CN + typing config."""
        return self.build_typing_plan().describe()

    def describe_benchmark_summary(self) -> str:
        """Describe the end-to-end benchmark summary."""
        return self.build_benchmark_summary().describe()

    def benchmark_summary_json(self) -> str:
        """Serialize the end-to-end benchmark summary as JSON."""
        return self.build_benchmark_summary().to_json()

    def describe_comparison_summary(self) -> str:
        """Describe the old-vs-new planning summary."""
        return self.build_comparison_summary().describe()

    def comparison_summary_json(self) -> str:
        """Serialize the old-vs-new planning summary as JSON."""
        return self.build_comparison_summary().to_json()
