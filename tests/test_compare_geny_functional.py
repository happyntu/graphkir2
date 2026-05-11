from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from compare_geny_functional import (  # noqa: E402
    normalize_geny_allele,
    parse_geny_output,
    run_comparison,
    write_allele_tsv,
)


def test_parse_geny_output_normalizes_gene_prefixed_and_bare_alleles(
    tmp_path: Path,
) -> None:
    geny_output = tmp_path / "geny.txt"
    geny_output.write_text(
        "\n".join(
            [
                "Sample: HG001",
                "[ilp] KIR2DL1 0030201",
                "[ilp] KIR2DL1 KIR2DL1*0040102",
                "[ilp] KIR2DS4 *0010101",
                "Sample: HG002",
                "[ilp] KIR2DL1 KIR2DL1*0020101+",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    assert normalize_geny_allele("KIR2DL1", "0030201") == "KIR2DL1*0030201"
    assert parse_geny_output(geny_output) == {
        "HG001": [
            "KIR2DL1*0030201",
            "KIR2DL1*0040102",
            "KIR2DS4*0010101",
        ],
        "HG002": ["KIR2DL1*0020101"],
    }


def test_write_allele_tsv_uses_graphkir_compatible_format(tmp_path: Path) -> None:
    output_tsv = tmp_path / "geny.allele.tsv"

    write_allele_tsv({"HG002": ["KIR2DL1*001"], "HG001": []}, output_tsv)

    assert output_tsv.read_text(encoding="utf-8").splitlines() == [
        "id\talleles",
        "HG001\t",
        "HG002\tKIR2DL1*001",
    ]


def test_run_comparison_evaluates_geny_and_optional_graphkir(
    tmp_path: Path,
) -> None:
    truth_tsv = tmp_path / "truth.tsv"
    truth_tsv.write_text(
        "\n".join(
            [
                "id\talleles",
                "HG001\tKIR2DL1*0030201_KIR2DL1*0040101_KIR2DS4*0010101",
                "HG002\tKIR2DL1*0010101",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    geny_output = tmp_path / "geny.txt"
    geny_output.write_text(
        "\n".join(
            [
                "Sample: HG001",
                "[ilp] KIR2DL1 KIR2DL1*0030201",
                "[ilp] KIR2DL1 KIR2DL1*0040102",
                "[ilp] KIR2DS4 KIR2DS4*0010101",
                "Sample: HG002",
                "[ilp] KIR2DL1 KIR2DL1*0020101",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    graphkir_tsv = tmp_path / "graphkir.tsv"
    graphkir_tsv.write_text(
        "\n".join(
            [
                "id\talleles",
                "HG001\tKIR2DL1*0030201_KIR2DL1*0040101_KIR2DS4*0010101",
                "HG002\tKIR2DL1*0010101",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "out"

    run_comparison(
        truth_tsv,
        geny_output,
        graphkir_tsv,
        output_dir,
        output_dir / "geny.normalized.tsv",
    )

    summary = (output_dir / "summary.tsv").read_text(encoding="utf-8")
    assert "Geny\tevaluated" in summary
    assert "graphkir2\tevaluated" in summary
    assert "\t0.750000\t0.750000\t0.750000\t" in summary
    assert "\t0.500000\t0.500000\t0.500000\t" in summary
    assert "\t1.000000\t1.000000\t1.000000\t" in summary
    assert (output_dir / "per_gene.tsv").exists()
    assert (output_dir / "notes.md").exists()


def test_run_comparison_writes_missing_input_report_without_geny(
    tmp_path: Path,
) -> None:
    truth_tsv = tmp_path / "truth.tsv"
    truth_tsv.write_text("id\talleles\nHG001\tKIR2DL1*001\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    run_comparison(
        truth_tsv,
        tmp_path / "missing_geny.txt",
        None,
        output_dir,
        output_dir / "geny.normalized.tsv",
    )

    missing = (output_dir / "missing_inputs.tsv").read_text(encoding="utf-8")
    assert "geny_output" in missing
    assert "provide raw Geny output" in missing
    assert (output_dir / "summary.tsv").read_text(encoding="utf-8").count("\n") == 1
    assert not (output_dir / "geny.normalized.tsv").exists()
