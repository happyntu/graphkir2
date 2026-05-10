from __future__ import annotations

from dataclasses import dataclass

from graphkir2.cn.interface import (
    CohortCopyNumberPlan,
    CopyNumberPlan,
    SampleCopyNumberPlan,
)
from graphkir2.config.runtime import TypingConfig
from graphkir2.mapping.interface import IndexBuildPlan, MappingPlan, SampleMappingPlan
from graphkir2.typing.interface import AlleleTyper
from graphkir2.typing.private_support import (
    collect_variant_support,
    neutralize_cross_gene_reads,
    parse_gene_groups,
    select_with_private_support,
)


@dataclass
class DummyRead:
    backbone: str
    l_sam: str
    r_sam: str
    weight: float
    ambiguous_weight: float
    lpv: list[str]
    rpv: list[str]
    lnv: list[str]
    rnv: list[str]


class DummyResult:
    value = [100.0, 98.0]
    allele_name = [["KIR2DS3*0020101", "KIR2DS3*0011201"], ["KIR2DS3*0010301", "KIR2DS3*0011201"]]

    def isFail(self) -> bool:
        return False

    def selectBest(self, min_fraction_ratio: float = 0.5) -> list[str]:
        del min_fraction_ratio
        return self.allele_name[0]


def test_private_support_neutralizes_configured_cross_gene_group_only() -> None:
    reads = [
        DummyRead("KIR2DS3*BACKBONE", "pair1\tA", "", 0.5, 0.5, [], [], [], []),
        DummyRead("KIR2DS5*BACKBONE", "pair1\tB", "", 0.5, 0.5, [], [], [], []),
        DummyRead("KIR2DS4*BACKBONE", "pair2\tA", "", 0.5, 0.5, [], [], [], []),
        DummyRead("KIR2DS5*BACKBONE", "pair2\tB", "", 0.5, 0.5, [], [], [], []),
    ]

    neutralize_cross_gene_reads(reads, parse_gene_groups("KIR2DS3/KIR2DS5"))

    assert [read.weight for read in reads] == [0.0, 0.0, 0.5, 0.5]
    assert [read.ambiguous_weight for read in reads] == [1.0, 1.0, 0.5, 0.5]


def test_private_support_can_select_lower_likelihood_supported_candidate() -> None:
    reads = [
        DummyRead("KIR2DS3*BACKBONE", "pair1\tA", "", 1.0, 0.0, ["v_truth"], [], ["v_false"], []),
        DummyRead("KIR2DS3*BACKBONE", "pair2\tA", "", 1.0, 0.0, ["v_truth"], [], ["v_false"], []),
        DummyRead("KIR2DS3*BACKBONE", "pair3\tA", "", 1.0, 0.0, [], [], ["v_false"], []),
        DummyRead("KIR2DS3*BACKBONE", "pair4\tA", "", 1.0, 0.0, [], [], ["v_false"], []),
        DummyRead("KIR2DS3*BACKBONE", "pair5\tA", "", 1.0, 0.0, [], [], ["v_false"], []),
    ]
    positive, negative = collect_variant_support(reads)
    allele_variants = {
        "KIR2DS3*0020101": {"v_false"},
        "KIR2DS3*0010301": {"v_truth"},
        "KIR2DS3*0011201": set(),
    }

    selected = select_with_private_support(
        DummyResult(),
        allele_variants,
        positive,
        negative,
        support_lambda=10.0,
        support_window=50.0,
        min_fraction_ratio=0.7,
    )

    assert selected == ["KIR2DS3*0010301", "KIR2DS3*0011201"]


def test_typing_plan_carries_private_support_config() -> None:
    mapping = MappingPlan(
        index=IndexBuildPlan("msa", "ref", "graph", "wgs"),
        multi_map_policy="likelihood",
        samples=(
            SampleMappingPlan(
                sample_id="01",
                output_prefix="out/01",
                source_read1="r1.fq.gz",
                source_read2="r2.fq.gz",
                extraction_required=False,
                extracted_read1="r1.fq.gz",
                extracted_read2="r2.fq.gz",
                graph_bam="out/01.bam",
                variant_prefix="out/01.variant",
                evidence_bam="out/01.variant.likelihood.bam",
                depth_tsv="out/01.variant.likelihood.depth.tsv",
            ),
        ),
    )
    copy_number = CopyNumberPlan(
        mode="per-sample",
        select_mode="p75",
        algorithm="LCND",
        base_dev=0.08,
        cohort=CohortCopyNumberPlan("cohort", "", "cohort.cn.tsv"),
        samples=(
            SampleCopyNumberPlan(
                sample_id="01",
                depth_tsv="out/01.variant.likelihood.depth.tsv",
                cn_tsv="out/01.cn.tsv",
                model_json="out/01.cn.json",
                diploid_gene="",
            ),
        ),
    )
    config = TypingConfig(
        cross_gene_neutralization_groups="KIR2DS3/KIR2DS5",
        private_support_genes="KIR2DS3",
        private_support_lambda=10.0,
        private_support_window=50.0,
    )

    plan = AlleleTyper("typing").plan(mapping, copy_number, config, "cohort")

    assert plan.cross_gene_neutralization_groups == "KIR2DS3/KIR2DS5"
    assert plan.private_support_genes == "KIR2DS3"
    assert plan.private_support_lambda == 10.0
    assert plan.private_support_window == 50.0
    assert plan.samples[0].private_support_genes == "KIR2DS3"
