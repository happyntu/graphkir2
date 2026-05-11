from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2dl5_remaining_failures import (  # noqa: E402
    classify_kir2dl5_cause,
    find_candidate_rank,
    gene_copy_counts,
    gene_copy_text,
    merged_gene_alleles,
    resolve_gene_top_n,
    typing_backbone_alias,
)


class DummyResult:
    def __init__(self) -> None:
        self.value = [-10.0, -12.5, -13.0]
        self.allele_name = [
            ["KIR2DL5A*028", "KIR2DL5B*01301"],
            ["KIR2DL5A*028", "KIR2DL5A*0010901"],
            ["KIR2DL5A*00107", "KIR2DL5A*01201"],
        ]

    def isFail(self) -> bool:
        return False


def test_typing_backbone_alias_merges_kir2dl5_a_b() -> None:
    assert typing_backbone_alias("KIR2DL5A*BACKBONE") == "KIR2DL5*BACKBONE"
    assert typing_backbone_alias("KIR2DL5B*01301") == "KIR2DL5*BACKBONE"
    assert typing_backbone_alias("KIR2DS3*BACKBONE") == "KIR2DS3*BACKBONE"


def test_merged_gene_alleles_keeps_both_kir2dl5_loci() -> None:
    alleles = [
        "KIR2DL5A*028",
        "KIR2DL5B*01301",
        "KIR2DS3*00201",
    ]

    assert merged_gene_alleles(alleles) == ["KIR2DL5A*028", "KIR2DL5B*01301"]


def test_gene_copy_text_includes_merged_count() -> None:
    counts = gene_copy_counts(["KIR2DL5A*028", "KIR2DL5A*0010901"])

    assert gene_copy_text(counts) == "KIR2DL5A=2;KIR2DL5B=0;merged=2"


def test_classify_kir2dl5_cause_finds_ab_assignment_mismatch() -> None:
    truth = ["KIR2DL5A*028", "KIR2DL5A*0010901"]
    candidate = ["KIR2DL5A*028", "KIR2DL5B*01301"]

    assert (
        classify_kir2dl5_cause(truth, candidate)
        == "ab_assignment_mismatch_on_merged_backbone"
    )


def test_classify_kir2dl5_cause_finds_shared_substitution() -> None:
    truth = ["KIR2DL5A*00107", "KIR2DL5A*00107"]
    candidate = ["KIR2DL5A*00107", "KIR2DL5A*01201"]

    assert (
        classify_kir2dl5_cause(truth, candidate, discard=candidate)
        == "shared_allele_substitution"
    )


def test_classify_kir2dl5_cause_finds_merged_copy_count_mismatch() -> None:
    assert (
        classify_kir2dl5_cause(["KIR2DL5A*00107"], ["KIR2DL5A*00107", "KIR2DL5A*01201"])
        == "merged_copy_count_mismatch"
    )


def test_resolve_gene_top_n_uses_base_top_n_for_non_target_gene() -> None:
    summary_row = {
        "top_n": "5000",
        "base_top_n": "600",
        "gene_base_top_ns": "KIR2DL1:1000",
    }

    assert resolve_gene_top_n("KIR2DL5", summary_row) == 600
    assert resolve_gene_top_n("KIR2DL1", summary_row) == 1000


def test_find_candidate_rank_reports_gap_from_top_value() -> None:
    rank, value, gap = find_candidate_rank(
        DummyResult(),
        ["KIR2DL5A*0010901", "KIR2DL5A*028"],
    )

    assert rank == "2"
    assert value == "-12.500"
    assert gap == "2.500"
