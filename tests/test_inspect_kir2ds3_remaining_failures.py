from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2ds3_remaining_failures import (  # noqa: E402
    classify_kir2ds3_cause,
    gene_alleles_list,
    gene_copy_text,
    resolve_kir2ds3_top_n,
)
from inspect_kir2ds5_remaining_failures import summarize_variant_rows  # noqa: E402


def test_gene_alleles_list_keeps_only_kir2ds3() -> None:
    alleles = [
        "KIR2DS3*011",
        "KIR2DS5*0270102",
        "KIR2DS3*0020101",
    ]

    assert gene_alleles_list(alleles) == ["KIR2DS3*011", "KIR2DS3*0020101"]


def test_gene_copy_text_counts_kir2ds3_only() -> None:
    assert gene_copy_text(["KIR2DS3*011", "KIR2DS5*0270102"]) == "KIR2DS3=1"


def test_classify_kir2ds3_cause_finds_shared_with_discard() -> None:
    truth = ["KIR2DS3*011", "KIR2DS3*016"]
    candidate = ["KIR2DS3*011", "KIR2DS3*0020101"]

    assert (
        classify_kir2ds3_cause(truth, candidate, discard=candidate, likelihood=candidate)
        == "shared_with_discard"
    )


def test_classify_kir2ds3_cause_finds_all_methods_shifted() -> None:
    truth = ["KIR2DS3*0010306", "KIR2DS3*0011302"]
    candidate = ["KIR2DS3*0010306", "KIR2DS3*0010301"]
    discard = ["KIR2DS3*0011302", "KIR2DS3*00108"]
    likelihood = ["KIR2DS3*0010306", "KIR2DS3*0020101"]

    assert (
        classify_kir2ds3_cause(truth, candidate, discard=discard, likelihood=likelihood)
        == "all_methods_disagree_or_shifted"
    )


def test_classify_kir2ds3_cause_prioritizes_candidate_regression() -> None:
    truth = ["KIR2DS3*011", "KIR2DS3*016"]
    candidate = ["KIR2DS3*011", "KIR2DS3*0020101"]

    assert (
        classify_kir2ds3_cause(truth, candidate, discard=truth, likelihood=candidate)
        == "candidate_regression"
    )


def test_resolve_kir2ds3_top_n_uses_private_support_top_n() -> None:
    summary_row = {
        "top_n": "5000",
        "base_top_n": "600",
        "gene_base_top_ns": "KIR2DL1:1000",
    }

    assert resolve_kir2ds3_top_n(summary_row) == 5000


def test_summarize_variant_rows_can_identify_kir2ds3_truth_margin() -> None:
    rows = [
        {
            "side": "truth_only",
            "in_exon": "yes",
            "support_class": "supported",
            "support_direction": "supports_truth",
            "net_weight": "5.000",
            "positive_weight": "5.000",
            "ambiguous_positive_weight": "2.500",
        },
        {
            "side": "candidate_only",
            "in_exon": "yes",
            "support_class": "unsupported",
            "support_direction": "supports_truth",
            "net_weight": "-7.000",
            "positive_weight": "0.000",
            "ambiguous_positive_weight": "0.000",
        },
    ]

    summary = summarize_variant_rows(rows)

    assert summary["truth_minus_candidate_net"] == "12.000"
    assert summary["truth_only_supported"] == "1"
    assert summary["candidate_only_unsupported"] == "1"
    assert summary["signal"] == "variant_signal_favors_truth"
