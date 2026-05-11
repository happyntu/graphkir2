from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from profile_targeted_top_n import (  # noqa: E402
    ProfileVariant,
    bundle_path,
    build_profile_variants,
    build_profile_variants_with_gene_overrides,
    parse_time_metrics,
    prediction_path,
    write_tsv,
)


def test_parse_time_metrics_extracts_runtime_and_rss() -> None:
    stderr = "tool stderr\n__GK_RUNTIME_SECONDS__=12.34\n__GK_MAXRSS_KB__=2048\n"

    runtime, rss = parse_time_metrics(stderr)

    assert runtime == 12.34
    assert rss == 2.0


def test_profile_paths_are_stable(tmp_path: Path) -> None:
    variant = ProfileVariant("targeted_base600", 600)

    assert prediction_path(tmp_path, "synthetic-x", variant) == (
        tmp_path / "synthetic-x.targeted_base600.allele.tsv"
    )
    assert bundle_path(tmp_path, "synthetic-x", variant) == (
        tmp_path / "synthetic-x.targeted_base600.bundle"
    )


def test_build_profile_variants_defaults_and_custom_values() -> None:
    assert build_profile_variants([]) == (
        ProfileVariant("full_top5000", 0),
        ProfileVariant("targeted_base600", 600),
        ProfileVariant("geneaware_base600_kir2dl1_1000", 600, "KIR2DL1:1000"),
    )
    assert build_profile_variants([0, 1000]) == (
        ProfileVariant("full_top5000", 0),
        ProfileVariant("targeted_base1000", 1000),
    )
    assert build_profile_variants_with_gene_overrides([600], ["KIR2DL1:1000"]) == (
        ProfileVariant("geneaware_base600_kir2dl1_1000", 600, "KIR2DL1:1000"),
    )


def test_write_tsv_uses_expected_columns(tmp_path: Path) -> None:
    path = tmp_path / "summary.tsv"
    write_tsv(
        path,
        [
            {
                "label": "synthetic-x",
                "variant": "full_top5000",
                "top_n": "5000",
                "base_top_n": "0",
                "gene_base_top_ns": "",
                "runtime_seconds": "1.000",
                "max_rss_mb": "100.0",
                "three_digit_f1": "1.000000",
                "five_digit_f1": "1.000000",
                "seven_digit_f1": "0.900000",
                "prediction_tsv": "pred.tsv",
                "bundle_dir": "bundle",
            }
        ],
    )

    assert path.read_text(encoding="utf-8").splitlines() == [
        "label\tvariant\ttop_n\tbase_top_n\tgene_base_top_ns\truntime_seconds\tmax_rss_mb\tthree_digit_f1\tfive_digit_f1\tseven_digit_f1\tprediction_tsv\tbundle_dir",
        "synthetic-x\tfull_top5000\t5000\t0\t\t1.000\t100.0\t1.000000\t1.000000\t0.900000\tpred.tsv\tbundle",
    ]
