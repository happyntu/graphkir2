from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2ds3_rank_failure import (  # noqa: E402
    contains_allele_prefix,
    functional_key,
    rank_kind,
    select_result_indexes,
)


class DummyResult:
    def __init__(self) -> None:
        self.allele_name = [
            ("KIR2DS3*011", "KIR2DS3*0020101"),
            ("KIR2DS3*011", "KIR2DS3*00108"),
            ("KIR2DS3*011", "KIR2DS3*016"),
        ]

    def isFail(self) -> bool:
        return False


def test_functional_key_preserves_copy_count_at_resolution() -> None:
    assert (
        functional_key(["KIR2DS3*011", "KIR2DS3*0020101"], 5)
        == "KIR2DS3*00201_KIR2DS3*011"
    )


def test_contains_allele_prefix_matches_prefix_only() -> None:
    alleles = ["KIR2DS3*011", "KIR2DS3*0020101"]

    assert contains_allele_prefix(alleles, "KIR2DS3*002")
    assert not contains_allele_prefix(alleles, "KIR2DS3*016")


def test_rank_kind_tags_truth_and_top_rows() -> None:
    truth = ["KIR2DS3*011", "KIR2DS3*016"]
    current = ["KIR2DS3*011", "KIR2DS3*0020101"]

    assert (
        rank_kind(
            current,
            truth,
            current,
            current,
            current,
            rank=1,
            top_detail_rows=10,
        )
        == "top;current;discard;likelihood"
    )
    assert (
        rank_kind(
            truth,
            truth,
            current,
            current,
            current,
            rank=692,
            top_detail_rows=10,
        )
        == "truth"
    )


def test_select_result_indexes_includes_top_and_requested_far_rank() -> None:
    result = DummyResult()

    assert select_result_indexes(
        result,
        [["KIR2DS3*011", "KIR2DS3*016"]],
        top_detail_rows=1,
    ) == [0, 2]
