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
    choose_targeted_top_n,
    collect_variant_support,
    neutralize_cross_gene_reads,
    parse_gene_groups,
    parse_name_set,
    private_positive_cross_gene_ratio,
    select_with_highest_suffix_tie_break,
    select_with_private_support,
    should_apply_conditional_private_support,
    should_use_discard_fallback,
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


class DummyTieResult:
    n = 2
    value = [-10.0, -10.0, -10.0, -11.0]
    allele_name = [
        ["KIR2DS4*0030103", "KIR2DS4*0010102"],
        ["KIR2DS4*0030103", "KIR2DS4*0010101"],
        ["KIR2DS4*0030103", "KIR2DS4*0010109"],
        ["KIR2DS4*0030101", "KIR2DS4*0010109"],
    ]
    fraction = [
        [0.52, 0.48],
        [0.52, 0.48],
        [0.52, 0.48],
        [0.52, 0.48],
    ]

    def isFail(self) -> bool:
        return False


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


def test_choose_targeted_top_n_keeps_high_top_n_for_target_genes() -> None:
    target_genes = frozenset({"KIR2DS3", "KIR2DS4"})

    assert choose_targeted_top_n("KIR2DS3", 5000, 600, target_genes) == 5000
    assert choose_targeted_top_n("KIR3DL3", 5000, 600, target_genes) == 600
    assert choose_targeted_top_n("KIR3DL3", 5000, None, target_genes) == 5000


def test_private_support_can_neutralize_only_target_gene_in_cross_gene_group() -> None:
    reads = [
        DummyRead("KIR2DS3*BACKBONE", "pair1\tA", "", 0.5, 0.5, [], [], [], []),
        DummyRead("KIR2DS5*BACKBONE", "pair1\tB", "", 0.5, 0.5, [], [], [], []),
    ]

    neutralize_cross_gene_reads(
        reads,
        parse_gene_groups("KIR2DS3/KIR2DS5"),
        target_genes=frozenset({"KIR2DS3"}),
    )

    assert [read.weight for read in reads] == [0.0, 0.5]
    assert [read.ambiguous_weight for read in reads] == [1.0, 0.5]


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


def test_conditional_private_support_uses_cross_gene_ratio_gate() -> None:
    reads = [
        DummyRead(
            "KIR2DS3*BACKBONE",
            "pair1\tA",
            "",
            4.0,
            0.0,
            ["v_private"],
            [],
            [],
            [],
        ),
        DummyRead("KIR2DS5*BACKBONE", "pair1\tB", "", 1.0, 0.0, [], [], [], []),
        DummyRead(
            "KIR2DS3*BACKBONE",
            "pair2\tA",
            "",
            1.0,
            0.0,
            ["v_private"],
            [],
            [],
            [],
        ),
    ]
    selected = ["KIR2DS3*0020101", "KIR2DS3*0011201"]
    allele_variants = {
        "KIR2DS3*0020101": {"v_private"},
        "KIR2DS3*0011201": set(),
    }
    groups = parse_gene_groups("KIR2DS3/KIR2DS5")

    ratio = private_positive_cross_gene_ratio(
        reads,
        selected,
        allele_variants,
        groups,
    )

    assert ratio == 0.8
    assert should_apply_conditional_private_support(
        selected,
        reads,
        allele_variants,
        parse_name_set("KIR2DS3*00201"),
        groups,
        0.8,
    )
    assert not should_apply_conditional_private_support(
        selected,
        reads,
        allele_variants,
        parse_name_set("KIR2DS3*00201"),
        groups,
        0.85,
    )
    assert not should_apply_conditional_private_support(
        ["KIR2DS3*0010301", "KIR2DS3*0011201"],
        reads,
        allele_variants,
        parse_name_set("KIR2DS3*00201"),
        groups,
        0.8,
    )


def test_discard_fallback_detects_residual_and_low_ratio_introduced_calls() -> None:
    assert should_use_discard_fallback(
        ["KIR2DS3*0010311", "KIR2DS3*0020101"],
        ["KIR2DS3*0010311", "KIR2DS3*0020101"],
        parse_name_set("KIR2DS3*00201"),
        parse_name_set("KIR2DS3*00103"),
        cross_gene_ratio=0.81,
        introduced_max_cross_gene_ratio=0.90,
        private_support=-30.0,
        max_private_support=-20.0,
        residual_min_cross_gene_ratio=0.70,
    )
    assert should_use_discard_fallback(
        ["KIR2DS3*012", "KIR2DS3*0010301"],
        ["KIR2DS3*012", "KIR2DS3*00108"],
        parse_name_set("KIR2DS3*00201"),
        parse_name_set("KIR2DS3*00103"),
        cross_gene_ratio=0.88,
        introduced_max_cross_gene_ratio=0.90,
        private_support=-30.0,
        max_private_support=-20.0,
        residual_min_cross_gene_ratio=0.70,
    )
    assert not should_use_discard_fallback(
        ["KIR2DS3*0011201", "KIR2DS3*0010301"],
        ["KIR2DS3*0011201", "KIR2DS3*0020101"],
        parse_name_set("KIR2DS3*00201"),
        parse_name_set("KIR2DS3*00103"),
        cross_gene_ratio=0.91,
        introduced_max_cross_gene_ratio=0.90,
        private_support=-30.0,
        max_private_support=-20.0,
        residual_min_cross_gene_ratio=0.70,
    )
    assert not should_use_discard_fallback(
        ["KIR2DS3*0010311", "KIR2DS3*0020101"],
        ["KIR2DS3*0010311", "KIR2DS3*0020101"],
        parse_name_set("KIR2DS3*00201"),
        parse_name_set("KIR2DS3*00103"),
        cross_gene_ratio=0.81,
        introduced_max_cross_gene_ratio=0.90,
        private_support=10.0,
        max_private_support=-20.0,
        residual_min_cross_gene_ratio=0.70,
    )
    assert not should_use_discard_fallback(
        ["KIR2DS3*011", "KIR2DS3*0020101"],
        ["KIR2DS3*011", "KIR2DS3*0020101"],
        parse_name_set("KIR2DS3*00201"),
        parse_name_set("KIR2DS3*00103"),
        cross_gene_ratio=0.61,
        introduced_max_cross_gene_ratio=0.90,
        private_support=-30.0,
        max_private_support=-20.0,
        residual_min_cross_gene_ratio=0.70,
    )


def test_highest_suffix_tie_break_keeps_same_five_digit_call() -> None:
    selected = select_with_highest_suffix_tie_break(
        DummyTieResult(),
        ["KIR2DS4*0030103", "KIR2DS4*0010102"],
        min_fraction_ratio=0.7,
    )

    assert selected == ["KIR2DS4*0030103", "KIR2DS4*0010109"]


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
        base_top_n=600,
        private_support_genes="KIR2DS3",
        private_support_lambda=10.0,
        private_support_window=50.0,
        private_support_condition_alleles="KIR2DS3*00201",
        private_support_cross_gene_ratio=0.8,
        private_support_discard_fallback_genes="KIR2DS3",
        private_support_discard_fallback_residual_alleles="KIR2DS3*00201",
        private_support_discard_fallback_introduced_alleles="KIR2DS3*00103",
        private_support_discard_fallback_introduced_max_ratio=0.9,
        private_support_discard_fallback_max_score=-20.0,
        private_support_discard_fallback_residual_min_ratio=0.7,
        highest_suffix_tie_break_genes="KIR2DS4",
    )

    plan = AlleleTyper("typing").plan(mapping, copy_number, config, "cohort")

    assert plan.cross_gene_neutralization_groups == "KIR2DS3/KIR2DS5"
    assert plan.base_top_n == 600
    assert plan.private_support_genes == "KIR2DS3"
    assert plan.private_support_lambda == 10.0
    assert plan.private_support_window == 50.0
    assert plan.private_support_condition_alleles == "KIR2DS3*00201"
    assert plan.private_support_cross_gene_ratio == 0.8
    assert plan.private_support_discard_fallback_genes == "KIR2DS3"
    assert plan.private_support_discard_fallback_residual_alleles == "KIR2DS3*00201"
    assert plan.private_support_discard_fallback_introduced_alleles == "KIR2DS3*00103"
    assert plan.private_support_discard_fallback_introduced_max_ratio == 0.9
    assert plan.private_support_discard_fallback_max_score == -20.0
    assert plan.private_support_discard_fallback_residual_min_ratio == 0.7
    assert plan.highest_suffix_tie_break_genes == "KIR2DS4"
    assert plan.samples[0].private_support_genes == "KIR2DS3"
    assert plan.samples[0].base_top_n == 600
    assert plan.samples[0].private_support_condition_alleles == "KIR2DS3*00201"
    assert plan.samples[0].private_support_cross_gene_ratio == 0.8
    assert plan.samples[0].private_support_discard_fallback_genes == "KIR2DS3"
    assert plan.samples[0].private_support_discard_fallback_introduced_max_ratio == 0.9
    assert plan.samples[0].private_support_discard_fallback_max_score == -20.0
    assert plan.samples[0].private_support_discard_fallback_residual_min_ratio == 0.7
    assert plan.samples[0].highest_suffix_tie_break_genes == "KIR2DS4"
