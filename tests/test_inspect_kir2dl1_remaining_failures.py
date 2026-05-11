from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2dl1_remaining_failures import (  # noqa: E402
    classify_kir2dl1_cause,
    gene_alleles_list,
    gene_copy_text,
    resolve_kir2dl1_top_n,
)


def test_gene_alleles_list_keeps_only_kir2dl1() -> None:
    alleles = [
        "KIR2DL1*0030242",
        "KIR2DS1*0020101",
        "KIR2DL1*00303",
    ]

    assert gene_alleles_list(alleles) == ["KIR2DL1*0030242", "KIR2DL1*00303"]


def test_gene_copy_text_counts_kir2dl1_only() -> None:
    assert gene_copy_text(["KIR2DL1*0030242", "KIR2DS1*0020101"]) == "KIR2DL1=1"


def test_classify_kir2dl1_cause_finds_shared_with_discard() -> None:
    truth = ["KIR2DL1*0030242", "KIR2DL1*0030229"]
    candidate = ["KIR2DL1*0030242", "KIR2DL1*00303"]

    assert (
        classify_kir2dl1_cause(truth, candidate, discard=candidate, likelihood=[])
        == "shared_with_discard"
    )


def test_classify_kir2dl1_cause_prioritizes_shared_with_discard() -> None:
    truth = ["KIR2DL1*0030242", "KIR2DL1*0030229"]
    candidate = ["KIR2DL1*0030242", "KIR2DL1*00303"]

    assert (
        classify_kir2dl1_cause(truth, candidate, discard=candidate, likelihood=truth)
        == "shared_with_discard"
    )


def test_resolve_kir2dl1_top_n_uses_gene_specific_limit() -> None:
    summary_row = {
        "top_n": "5000",
        "base_top_n": "600",
        "gene_base_top_ns": "KIR2DL1:1000,KIR2DS3:5000",
    }

    assert resolve_kir2dl1_top_n(summary_row) == 1000
