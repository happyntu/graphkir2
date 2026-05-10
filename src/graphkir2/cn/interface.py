"""Copy-number stage interfaces for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config.runtime import CopyNumberConfig
from ..mapping.interface import MappingPlan


@dataclass(frozen=True)
class SampleCopyNumberPlan:
    """Planned CN artifacts for one sample."""

    sample_id: str
    depth_tsv: str
    cn_tsv: str
    model_json: str
    diploid_gene: str


@dataclass(frozen=True)
class CohortCopyNumberPlan:
    """Planned cohort-level CN artifacts."""

    cohort_prefix: str
    model_json: str
    merged_cn_tsv: str


@dataclass(frozen=True)
class CopyNumberPlan:
    """Top-level copy-number plan."""

    mode: str
    select_mode: str
    algorithm: str
    base_dev: float
    cohort: CohortCopyNumberPlan
    samples: tuple[SampleCopyNumberPlan, ...]

    def describe(self) -> str:
        """Render a readable CN-plan summary."""
        lines = [
            "graphkir2 copy-number plan",
            f"mode={self.mode}",
            f"select_mode={self.select_mode}",
            f"algorithm={self.algorithm}",
            f"base_dev={self.base_dev}",
            f"cohort_prefix={self.cohort.cohort_prefix}",
            f"cohort_model_json={self.cohort.model_json}",
            f"merged_cn_tsv={self.cohort.merged_cn_tsv}",
            f"samples={len(self.samples)}",
        ]
        for sample in self.samples:
            lines.extend(
                [
                    f"- sample={sample.sample_id}",
                    f"  depth_tsv={sample.depth_tsv}",
                    f"  cn_tsv={sample.cn_tsv}",
                    f"  model_json={sample.model_json}",
                    f"  diploid_gene={sample.diploid_gene or '<none>'}",
                ]
            )
        return "\n".join(lines)


@dataclass(frozen=True)
class CopyNumberEstimator:
    """Marker/planner object for the copy-number stage."""

    stage_name: str

    def plan(
        self,
        mapping: MappingPlan,
        config: CopyNumberConfig,
        cohort_name: str,
        select_mode: str = "p75",
        algorithm: str = "LCND",
        base_dev: float = 0.08,
    ) -> CopyNumberPlan:
        """Build a non-executing CN plan from mapping outputs."""
        merged_cn_tsv = cohort_name + ".cn.tsv"
        sample_plans: list[SampleCopyNumberPlan] = []

        if config.cohort_mode:
            suffix = f".{select_mode}.cohort.{algorithm}"
            cohort = CohortCopyNumberPlan(
                cohort_prefix=cohort_name + suffix,
                model_json=cohort_name + suffix + ".json",
                merged_cn_tsv=merged_cn_tsv,
            )
            for sample in mapping.samples:
                cn_tsv = str(Path(sample.depth_tsv).with_suffix(suffix + ".tsv"))
                sample_plans.append(
                    SampleCopyNumberPlan(
                        sample_id=sample.sample_id,
                        depth_tsv=sample.depth_tsv,
                        cn_tsv=cn_tsv,
                        model_json=cohort.model_json,
                        diploid_gene="",
                    )
                )
            mode = "cohort"
        else:
            cohort = CohortCopyNumberPlan(
                cohort_prefix=cohort_name,
                model_json="",
                merged_cn_tsv=merged_cn_tsv,
            )
            for sample in mapping.samples:
                suffix = f".{select_mode}.{algorithm}"
                cn_prefix = str(Path(sample.depth_tsv).with_suffix(suffix))
                sample_plans.append(
                    SampleCopyNumberPlan(
                        sample_id=sample.sample_id,
                        depth_tsv=sample.depth_tsv,
                        cn_tsv=cn_prefix + ".tsv",
                        model_json=cn_prefix + ".json",
                        diploid_gene=config.diploid_gene,
                    )
                )
            mode = "per-sample"

        return CopyNumberPlan(
            mode=mode,
            select_mode=select_mode,
            algorithm=algorithm,
            base_dev=base_dev,
            cohort=cohort,
            samples=tuple(sample_plans),
        )
