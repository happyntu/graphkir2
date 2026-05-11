"""Legacy graphkir command planning helpers."""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.pipeline import GraphKir2RunConfig


@dataclass(frozen=True)
class LegacyCommandPlan:
    """Equivalent legacy CLI invocation derived from graphkir2 config."""

    program: str
    args: list[str]
    exact_config_match: bool = True
    notes: tuple[str, ...] = ()

    def render_shell(self) -> str:
        """Render the legacy command as a shell-safe string."""
        return shlex.join([self.program, *self.args])


def build_legacy_graphkir_command(config: GraphKir2RunConfig) -> LegacyCommandPlan:
    """Translate graphkir2 scaffold config into an equivalent legacy command."""
    args: list[str] = [
        "--thread",
        str(config.runtime.threads),
        "--engine",
        config.engine.name,
        "--ref-genome",
        config.reference.genome,
        "--index-folder",
        config.index.index_folder,
        "--ipd-version",
        config.index.ipd_version,
        "--msa-type",
        config.index.msa_type,
        "--allele-strategy",
        config.typing.strategy,
    ]
    if config.inputs.manifest_path:
        args.extend(["--input-csv", config.inputs.manifest_path])
    if config.inputs.output_folder:
        args.extend(["--output-folder", config.inputs.output_folder])
    if config.inputs.output_cohort_name:
        args.extend(["--output-cohort-name", config.inputs.output_cohort_name])
    if not config.index.use_exon_only_alleles:
        args.append("--msa-no-exon-only-allele")
    if config.copy_number.diploid_gene:
        args.extend(["--cn-diploid-gene", config.copy_number.diploid_gene])
    if config.copy_number.cohort_mode:
        args.append("--cn-cohort")
    if not config.copy_number.assume_3dl3_diploid:
        args.append("--cn-3dl3-not-diploid")
    if config.inputs.skip_extraction:
        args.append("--step-skip-extraction")
    notes: list[str] = []
    exact_config_match = True
    if config.mapping.multi_map_policy != "discard":
        exact_config_match = False
        notes.append(
            "legacy graphkir only supports discard-style multi-map handling; "
            "this command is the closest baseline, not an exact execution match"
        )
    if config.typing.exon_weight != 1.0:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no exon-weight flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.gene_exon_weights:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no gene-specific exon-weight flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.mapping.multi_map_policy == "margin":
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no margin-based multi-map mode; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.mapping.multi_map_policy == "likelihood":
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no ambiguity-likelihood multi-map mode; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.margin_gap != 5.0 or config.typing.margin_scale != 2.0:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no margin-tuning flags; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.ambiguity_neutral_prob != 0.999:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no ambiguity-likelihood neutral-probability flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.select_min_fraction_ratio != 0.5:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no candidate fraction-threshold flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.cross_gene_neutralization_groups:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no cross-gene ambiguity neutralization flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if (
        config.typing.private_support_genes
        or config.typing.private_support_lambda != 0.0
        or config.typing.private_support_window != 0.0
        or config.typing.private_support_condition_alleles
        or config.typing.private_support_cross_gene_ratio != 0.0
        or config.typing.private_support_discard_fallback_genes
        or config.typing.private_support_discard_fallback_residual_alleles
        or config.typing.private_support_discard_fallback_introduced_alleles
        or config.typing.private_support_discard_fallback_introduced_max_ratio != 0.0
    ):
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no private-support reranking flags; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    if config.typing.highest_suffix_tie_break_genes:
        exact_config_match = False
        notes.append(
            "legacy graphkir CLI has no highest-suffix tie-break flag; "
            "this graphkir2 config represents a refactor-only typing enhancement"
        )
    return LegacyCommandPlan(
        program="graphkir",
        args=args,
        exact_config_match=exact_config_match,
        notes=tuple(notes),
    )
