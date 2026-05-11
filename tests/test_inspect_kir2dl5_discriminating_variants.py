from dataclasses import dataclass, field
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from inspect_kir2dl5_discriminating_variants import (  # noqa: E402
    VariantEvidence,
    build_allele_variant_index,
    carriers_of_interest,
    compute_variant_evidence,
    discriminating_variant_ids,
    summarize_sample_rows,
    support_class,
    support_direction,
)


@dataclass
class DummyVariant:
    id: str
    allele: list[str]
    pos: int = 10
    typ: str = "single"
    val: str = "A"
    in_exon: bool = True


@dataclass
class DummyRead:
    name: str
    weight: float = 1.0
    ambiguous_weight: float = 0.0
    lpv: list[str] = field(default_factory=list)
    rpv: list[str] = field(default_factory=list)
    lnv: list[str] = field(default_factory=list)
    rnv: list[str] = field(default_factory=list)

    @property
    def l_sam(self) -> str:
        return f"{self.name}\t0\tKIR2DL5"

    @property
    def r_sam(self) -> str:
        return ""


def test_discriminating_variant_ids_compare_truth_and_candidate_genotypes() -> None:
    variants = [
        DummyVariant("v1", ["KIR2DL5A*00107"]),
        DummyVariant("v2", ["KIR2DL5A*01201"]),
        DummyVariant("shared", ["KIR2DL5A*00107", "KIR2DL5A*01201"]),
    ]
    allele_variants = build_allele_variant_index(variants)

    truth_only, candidate_only = discriminating_variant_ids(
        ["KIR2DL5A*00107"],
        ["KIR2DL5A*01201"],
        allele_variants,
    )

    assert truth_only == {"v1"}
    assert candidate_only == {"v2"}


def test_compute_variant_evidence_tracks_ambiguous_weight_and_read_counts() -> None:
    reads = [
        DummyRead("r1", weight=0.25, ambiguous_weight=0.75, lpv=["v1"]),
        DummyRead("r2", weight=1.0, lnv=["v1"]),
        DummyRead("r3", weight=0.5, ambiguous_weight=0.5, rpv=["v1"], rnv=["v1"]),
    ]

    evidence = compute_variant_evidence("v1", reads)

    assert evidence.positive_weight == 0.75
    assert evidence.negative_weight == 1.5
    assert evidence.positive_reads == 2
    assert evidence.negative_reads == 2
    assert evidence.ambiguous_positive_weight == 0.75
    assert evidence.ambiguous_negative_weight == 0.5
    assert evidence.ambiguous_positive_ratio == 1.0


def test_support_class_uses_private_support_thresholds() -> None:
    assert support_class(VariantEvidence(2.0, 1.0, 2, 1, 0.0, 0.0)) == "supported"
    assert support_class(VariantEvidence(0.0, 5.0, 0, 5, 0.0, 0.0)) == "unsupported"
    assert support_class(VariantEvidence(0.0, 0.0, 0, 0, 0.0, 0.0)) == "unobserved"


def test_support_direction_accounts_for_variant_owner() -> None:
    assert (
        support_direction("truth_only", VariantEvidence(3.0, 1.0, 3, 1, 0.0, 0.0))
        == "supports_truth"
    )
    assert (
        support_direction("candidate_only", VariantEvidence(0.0, 2.0, 0, 2, 0.0, 0.0))
        == "supports_truth"
    )


def test_carriers_of_interest_limits_to_truth_and_candidate_alleles() -> None:
    variant = DummyVariant(
        "v1",
        ["KIR2DL5A*00107", "KIR2DL5A*01201", "KIR2DL5A*999"],
    )

    assert carriers_of_interest(
        variant,
        ["KIR2DL5A*01201", "KIR2DL5A*00107"],
    ) == ["KIR2DL5A*00107", "KIR2DL5A*01201"]


def test_summarize_sample_rows_computes_truth_candidate_margin() -> None:
    rows = [
        {
            "panel": "panel",
            "sample_id": "01",
            "root_cause": "shared_allele_substitution",
            "side": "truth_only",
            "in_exon": "yes",
            "support_class": "supported",
            "net_weight": "4.000",
            "positive_weight": "5.000",
            "ambiguous_positive_weight": "1.000",
        },
        {
            "panel": "panel",
            "sample_id": "01",
            "root_cause": "shared_allele_substitution",
            "side": "candidate_only",
            "in_exon": "no",
            "support_class": "supported",
            "net_weight": "1.000",
            "positive_weight": "3.000",
            "ambiguous_positive_weight": "1.000",
        },
    ]

    summary = summarize_sample_rows(rows)

    assert summary[0]["truth_minus_candidate_net"] == "3.000"
    assert summary[0]["ambiguous_positive_ratio"] == "0.250"
    assert summary[0]["signal"] == "variant_signal_favors_truth"
