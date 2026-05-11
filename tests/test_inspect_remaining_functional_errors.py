from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_remaining_functional_errors import (  # noqa: E402
    build_error_rows,
    classify_error,
    extract_sample_id,
    gene_counter,
    normalize_allele,
)


def test_extract_sample_id_matches_generated_prediction_names() -> None:
    raw = "benchmarks/generated/synthetic-x/synthetic-x.05.index.variant.cn.full"

    assert extract_sample_id(raw) == "05"


def test_normalize_allele_strips_graphkir_suffixes() -> None:
    assert normalize_allele("KIR2DL1*0030242$") == "KIR2DL1*0030242"
    assert normalize_allele("KIR2DL1*0030242e01") == "KIR2DL1*0030242"
    assert normalize_allele("") == ""


def test_gene_counter_limits_resolution_and_keeps_copy_counts() -> None:
    alleles = ["KIR2DL1*0030242", "KIR2DL1*0030229", "KIR2DS1*0020101"]

    assert gene_counter(alleles, "KIR2DL1", 5) == Counter({"KIR2DL1*00302": 2})


def test_classify_error_prefers_candidate_regression_when_discard_is_correct() -> None:
    truth = Counter({"KIR2DL1*003": 2})
    candidate = Counter({"KIR2DL1*003": 1, "KIR2DL1*004": 1})

    assert classify_error(
        truth,
        candidate,
        discard=truth,
        likelihood=candidate,
        enhancedgate=candidate,
    ) == "candidate_regression"


def test_build_error_rows_compares_candidate_to_truth_and_methods() -> None:
    calls = {
        "panel": {
            "truth": {"00": ["KIR2DL1*0030242", "KIR2DL1*0030229"]},
            "discard": {"00": ["KIR2DL1*0030242", "KIR2DL1*00303"]},
            "likelihood_top5000": {"00": ["KIR2DL1*0030205", "KIR2DL1*0040110"]},
            "enhancedgate_geneaware": {"00": ["KIR2DL1*0030205", "KIR2DL1*0040110"]},
            "candidate": {"00": ["KIR2DL1*0030242", "KIR2DL1*00303"]},
        }
    }

    rows = build_error_rows(
        calls,
        candidate_method="candidate",
        methods=("discard", "likelihood_top5000", "enhancedgate_geneaware", "candidate"),
    )

    assert rows == [
        {
            "panel": "panel",
            "sample_id": "00",
            "gene": "KIR2DL1",
            "resolution": "5",
            "cause_hint": "shared_with_discard",
            "truth": "KIR2DL1*00302_KIR2DL1*00302",
            "candidate": "KIR2DL1*00302_KIR2DL1*00303",
            "missing": "KIR2DL1*00302",
            "extra": "KIR2DL1*00303",
            "truth_full": "KIR2DL1*0030242_KIR2DL1*0030229",
            "candidate_full": "KIR2DL1*0030242_KIR2DL1*00303",
            "discard": "KIR2DL1*00302_KIR2DL1*00303",
            "discard_full": "KIR2DL1*0030242_KIR2DL1*00303",
            "likelihood_top5000": "KIR2DL1*00302_KIR2DL1*00401",
            "likelihood_top5000_full": "KIR2DL1*0030205_KIR2DL1*0040110",
            "enhancedgate_geneaware": "KIR2DL1*00302_KIR2DL1*00401",
            "enhancedgate_geneaware_full": "KIR2DL1*0030205_KIR2DL1*0040110",
            "candidate": "KIR2DL1*00302_KIR2DL1*00303",
            "candidate_full": "KIR2DL1*0030242_KIR2DL1*00303",
        }
    ]
