from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2ds5_remaining_failures import (  # noqa: E402
    classify_kir2ds5_cause,
    gene_alleles_list,
    gene_copy_text,
    summarize_variant_rows,
)


def test_gene_alleles_list_keeps_only_kir2ds5() -> None:
    alleles = [
        "KIR2DS5*024",
        "KIR2DS3*00201",
        "KIR2DS5*0270102",
    ]

    assert gene_alleles_list(alleles) == ["KIR2DS5*024", "KIR2DS5*0270102"]


def test_gene_copy_text_counts_kir2ds5_only() -> None:
    assert gene_copy_text(["KIR2DS5*024", "KIR2DS3*00201"]) == "KIR2DS5=1"


def test_classify_kir2ds5_cause_finds_shared_with_discard() -> None:
    truth = ["KIR2DS5*024", "KIR2DS5*0020503"]
    candidate = ["KIR2DS5*0270102", "KIR2DS5*0020503"]

    assert (
        classify_kir2ds5_cause(truth, candidate, discard=candidate, likelihood=candidate)
        == "shared_with_discard"
    )


def test_classify_kir2ds5_cause_finds_unresolved_likelihood_pattern() -> None:
    truth = ["KIR2DS5*0210101", "KIR2DS5*0020129"]
    candidate = ["KIR2DS5*0210102", "KIR2DS5*0270102"]
    discard = ["KIR2DS5*034", "KIR2DS5*034"]

    assert (
        classify_kir2ds5_cause(truth, candidate, discard=discard, likelihood=candidate)
        == "unresolved_likelihood_pattern"
    )


def test_classify_kir2ds5_cause_prioritizes_candidate_regression() -> None:
    truth = ["KIR2DS5*024", "KIR2DS5*0020503"]
    candidate = ["KIR2DS5*0270102", "KIR2DS5*0020503"]

    assert (
        classify_kir2ds5_cause(truth, candidate, discard=truth, likelihood=candidate)
        == "candidate_regression"
    )


def test_classify_kir2ds5_cause_finds_copy_count_mismatch() -> None:
    assert (
        classify_kir2ds5_cause(["KIR2DS5*024"], ["KIR2DS5*024", "KIR2DS5*0270102"])
        == "copy_count_mismatch"
    )


def test_summarize_variant_rows_computes_truth_favoring_margin() -> None:
    rows = [
        {
            "side": "truth_only",
            "in_exon": "yes",
            "support_class": "supported",
            "support_direction": "supports_truth",
            "net_weight": "4.000",
            "positive_weight": "5.000",
            "ambiguous_positive_weight": "1.000",
        },
        {
            "side": "candidate_only",
            "in_exon": "no",
            "support_class": "unsupported",
            "support_direction": "supports_truth",
            "net_weight": "-3.000",
            "positive_weight": "0.000",
            "ambiguous_positive_weight": "0.000",
        },
    ]

    summary = summarize_variant_rows(rows)

    assert summary["truth_only_supported"] == "1"
    assert summary["candidate_only_unsupported"] == "1"
    assert summary["truth_minus_candidate_net"] == "7.000"
    assert summary["ambiguous_positive_ratio"] == "0.200"
    assert summary["signal"] == "variant_signal_favors_truth"


def test_summarize_variant_rows_computes_candidate_favoring_margin() -> None:
    rows = [
        {
            "side": "truth_only",
            "in_exon": "yes",
            "support_class": "unsupported",
            "support_direction": "supports_candidate",
            "net_weight": "-4.000",
            "positive_weight": "0.000",
            "ambiguous_positive_weight": "0.000",
        },
        {
            "side": "candidate_only",
            "in_exon": "yes",
            "support_class": "supported",
            "support_direction": "supports_candidate",
            "net_weight": "3.000",
            "positive_weight": "3.000",
            "ambiguous_positive_weight": "1.000",
        },
    ]

    summary = summarize_variant_rows(rows)

    assert summary["truth_minus_candidate_net"] == "-7.000"
    assert summary["supports_candidate_variants"] == "2"
    assert summary["signal"] == "variant_signal_favors_candidate"
