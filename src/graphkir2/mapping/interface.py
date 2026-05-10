"""Mapping-stage interfaces for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass

from ..config.runtime import (
    GraphKir2RuntimeConfig,
    IndexConfig,
    InputConfig,
    MappingConfig,
    ReferenceConfig,
)
from ..model import SampleInput


def sanitize_path_component(value: str) -> str:
    """Match the legacy suffix sanitization used by graphkir."""
    return value.replace(".", "_").replace("/", "_").replace("\\", "_")


def build_msa_prefix(index: IndexConfig) -> str:
    """Build the legacy-style MSA prefix from index configuration."""
    base = f"{index.index_folder}/kir_{index.ipd_version}"
    if index.use_exon_only_alleles:
        return f"{base}_withexon_{index.msa_type}.leftalign"
    return f"{base}_{index.msa_type}.leftalign"


def build_wgs_index_path(index_folder: str, genome: str) -> str:
    """Return the expected WGS index path used by legacy graphkir."""
    if genome == "hg19":
        return f"{index_folder}/hs37d5.fa.gz"
    if genome == "hg38":
        return f"{index_folder}/hs38noalt.fa.gz"
    raise ValueError(f"Unsupported reference genome: {genome}")


def build_mapping_evidence_suffix(policy: str) -> str:
    """Return the suffix applied to mapping evidence artifacts for a policy."""
    if policy == "discard":
        return ".no_multi"
    if policy == "best-only":
        return ".best"
    if policy == "weighted":
        return ".weighted"
    if policy == "margin":
        return ".margin"
    if policy == "likelihood":
        return ".likelihood"
    raise ValueError(f"Unsupported multi-map policy: {policy}")


@dataclass(frozen=True)
class IndexBuildPlan:
    """Planned reference/index artifacts for mapping."""

    msa_prefix: str
    index_ref: str
    graph_index: str
    wgs_index: str


@dataclass(frozen=True)
class SampleMappingPlan:
    """Planned mapping artifacts for one sample."""

    sample_id: str
    output_prefix: str
    source_read1: str
    source_read2: str
    extraction_required: bool
    extracted_read1: str
    extracted_read2: str
    graph_bam: str
    variant_prefix: str
    evidence_bam: str
    depth_tsv: str


@dataclass(frozen=True)
class MappingPlan:
    """Top-level mapping plan covering index and sample outputs."""

    index: IndexBuildPlan
    multi_map_policy: str
    samples: tuple[SampleMappingPlan, ...]

    def describe(self) -> str:
        """Render a readable mapping-plan summary."""
        lines = [
            "graphkir2 mapping plan",
            f"msa_prefix={self.index.msa_prefix}",
            f"index_ref={self.index.index_ref}",
            f"graph_index={self.index.graph_index}",
            f"wgs_index={self.index.wgs_index}",
            f"multi_map_policy={self.multi_map_policy}",
            f"samples={len(self.samples)}",
        ]
        for sample in self.samples:
            lines.extend(
                [
                    f"- sample={sample.sample_id}",
                    f"  output_prefix={sample.output_prefix}",
                    f"  source_reads={sample.source_read1} | {sample.source_read2}",
                    f"  extraction_required={sample.extraction_required}",
                    f"  extracted_reads={sample.extracted_read1} | {sample.extracted_read2}",
                    f"  graph_bam={sample.graph_bam}",
                    f"  variant_prefix={sample.variant_prefix}",
                    f"  evidence_bam={sample.evidence_bam}",
                    f"  depth_tsv={sample.depth_tsv}",
                ]
            )
        return "\n".join(lines)


@dataclass(frozen=True)
class ReadMapper:
    """Marker object for the mapping stage."""

    stage_name: str

    def plan(
        self,
        samples: tuple[SampleInput, ...],
        runtime: GraphKir2RuntimeConfig,
        index: IndexConfig,
        reference: ReferenceConfig,
        inputs: InputConfig,
        mapping: MappingConfig,
    ) -> MappingPlan:
        """Build a non-executing mapping plan from graphkir2 configuration."""
        del runtime
        msa_prefix = build_msa_prefix(index)
        index_ref = msa_prefix + ".mut01"
        graph_index = index_ref + ".graph"
        wgs_index = build_wgs_index_path(index.index_folder, reference.genome)
        index_plan = IndexBuildPlan(
            msa_prefix=msa_prefix,
            index_ref=index_ref,
            graph_index=graph_index,
            wgs_index=wgs_index,
        )

        graph_suffix = "." + sanitize_path_component(graph_index)
        wgs_suffix = "." + sanitize_path_component(wgs_index)
        evidence_suffix = build_mapping_evidence_suffix(mapping.multi_map_policy)
        sample_plans: list[SampleMappingPlan] = []
        for sample in samples:
            if inputs.skip_extraction:
                graph_input_prefix = sample.output_prefix
                extracted_read1 = sample.read1
                extracted_read2 = sample.read2
            else:
                graph_input_prefix = sample.output_prefix + wgs_suffix + ".extract"
                extracted_read1 = graph_input_prefix + ".read.1.fq.gz"
                extracted_read2 = graph_input_prefix + ".read.2.fq.gz"

            graph_prefix = graph_input_prefix + graph_suffix
            variant_prefix = graph_prefix + ".variant"
            evidence_bam = variant_prefix + evidence_suffix + ".bam"
            depth_tsv = variant_prefix + evidence_suffix + ".depth.tsv"
            sample_plans.append(
                SampleMappingPlan(
                    sample_id=sample.sample_id,
                    output_prefix=sample.output_prefix,
                    source_read1=sample.read1,
                    source_read2=sample.read2,
                    extraction_required=not inputs.skip_extraction,
                    extracted_read1=extracted_read1,
                    extracted_read2=extracted_read2,
                    graph_bam=graph_prefix + ".bam",
                    variant_prefix=variant_prefix,
                    evidence_bam=evidence_bam,
                    depth_tsv=depth_tsv,
                )
            )

        return MappingPlan(
            index=index_plan,
            multi_map_policy=mapping.multi_map_policy,
            samples=tuple(sample_plans),
        )
