from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from prepare_hprc_real_mini import (  # noqa: E402
    HprcSample,
    build_enhancedgate_preset,
    build_preset,
    read_hprc_samples,
    read_truth_sample_ids,
    select_available_samples,
    write_manifest,
)


def test_select_available_samples_requires_truth_and_fastq(tmp_path: Path) -> None:
    fastq_root = tmp_path / "fastq"
    fastq_root.mkdir()
    (fastq_root / "SRR1_1.fastq.gz").write_text("r1", encoding="utf-8")
    (fastq_root / "SRR1_2.fastq.gz").write_text("r2", encoding="utf-8")
    samples = [
        HprcSample("HG001", "SRR1"),
        HprcSample("HG002", "SRR2"),
        HprcSample("HG003", "SRR3"),
    ]

    available, missing = select_available_samples(
        samples,
        frozenset({"HG001", "HG002"}),
        (fastq_root,),
        max_samples=4,
    )

    assert [sample.sample.sample_id for sample in available] == ["HG001"]
    assert [sample.sample_id for sample in missing] == ["HG002"]


def test_read_samples_and_truth_tables(tmp_path: Path) -> None:
    cohort = tmp_path / "hprc.csv"
    cohort.write_text(
        "id,sample_id,t1k\nHG001,SRR1,0\nHG002,,1\nHG003,SRR3,0\n",
        encoding="utf-8",
    )
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "id\tname\talleles\nHG001\tHG001\tKIR2DL1*001\nHG003\tHG003\tKIR2DL1*002\n",
        encoding="utf-8",
    )

    assert read_hprc_samples(cohort) == [
        HprcSample("HG001", "SRR1"),
        HprcSample("HG003", "SRR3"),
    ]
    assert read_truth_sample_ids(truth) == frozenset({"HG001", "HG003"})


def test_write_manifest_and_build_preset(tmp_path: Path) -> None:
    read1 = tmp_path / "SRR1_1.fastq.gz"
    read2 = tmp_path / "SRR1_2.fastq.gz"
    read1.write_text("r1", encoding="utf-8")
    read2.write_text("r2", encoding="utf-8")
    sample = select_available_samples(
        [HprcSample("HG001", "SRR1")],
        frozenset({"HG001"}),
        (tmp_path,),
        max_samples=1,
    )[0][0]
    manifest = tmp_path / "manifest.csv"

    write_manifest([sample], manifest, output_dir=Path("benchmarks/results/hprc-real-mini"))

    assert manifest.read_text(encoding="utf-8").splitlines() == [
        "name,r1,r2,cnfile",
        f"benchmarks/results/hprc-real-mini/HG001,{read1},{read2},",
    ]
    preset = build_preset(
        manifest,
        Path("benchmarks/results/hprc-real-mini/cohort"),
        Path("truth.tsv"),
        skip_extraction=True,
    )
    assert preset["benchmark_label"] == "hprc-real-mini"
    assert preset["input_csv"] == str(manifest)
    assert preset["step_skip_extraction"] is True
    assert "allele_base_top_n" not in preset

    enhanced = build_enhancedgate_preset(preset)
    assert enhanced["benchmark_label"] == "hprc-real-mini-enhancedgate"
    assert enhanced["allele_base_top_n"] == 600
    assert enhanced["allele_gene_base_top_ns"] == "KIR2DL1:1000"
    assert enhanced["allele_private_support_genes"] == "KIR2DS3"
