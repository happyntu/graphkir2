from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from prepare_hprc_download_plan import (  # noqa: E402
    build_download_rows,
    render_download_script,
    select_truth_matched_samples,
    write_plan_tsv,
)
from prepare_hprc_real_mini import HprcSample  # noqa: E402


def test_select_truth_matched_samples_preserves_order_and_limit() -> None:
    samples = [
        HprcSample("HG001", "SRR1"),
        HprcSample("HG002", "SRR2"),
        HprcSample("HG003", "SRR3"),
    ]

    assert select_truth_matched_samples(samples, frozenset({"HG003", "HG001"}), 1) == [
        HprcSample("HG001", "SRR1")
    ]
    assert select_truth_matched_samples(samples, frozenset({"HG003", "HG001"}), 0) == [
        HprcSample("HG001", "SRR1"),
        HprcSample("HG003", "SRR3"),
    ]


def test_build_download_rows_uses_prepare_hprc_compatible_names(
    tmp_path: Path,
) -> None:
    rows = build_download_rows(
        [HprcSample("HG002", "SRR14724532")],
        fastq_root=tmp_path / "fastq",
        sra_cache=tmp_path / "sra",
        tmp_dir=tmp_path / "tmp",
        threads=6,
    )

    assert len(rows) == 1
    row = rows[0]
    assert row.read1 == tmp_path / "fastq" / "SRR14724532_1.fastq.gz"
    assert row.read2 == tmp_path / "fastq" / "SRR14724532_2.fastq.gz"
    assert row.sra_path == tmp_path / "sra" / "SRR14724532" / "SRR14724532.sra"
    assert "prefetch SRR14724532" in row.download_command
    assert "fasterq-dump" in row.download_command
    assert "--threads 6" in row.download_command
    assert "SRR14724532_1.fastq.gz" in row.download_command
    assert not row.read1_exists
    assert not row.read2_exists


def test_write_plan_tsv_and_render_script(tmp_path: Path) -> None:
    rows = build_download_rows(
        [HprcSample("HG002", "SRR14724532")],
        fastq_root=tmp_path / "fastq",
        sra_cache=tmp_path / "sra",
        tmp_dir=tmp_path / "tmp",
        threads=4,
    )
    plan_tsv = tmp_path / "plan.tsv"

    write_plan_tsv(rows, plan_tsv)
    plan = plan_tsv.read_text(encoding="utf-8")
    assert "id\tsample_id\tread1\tread2" in plan
    assert "HG002\tSRR14724532" in plan
    assert "false\tfalse" in plan

    script = render_download_script(rows, tmp_path / "fastq")
    assert "set -euo pipefail" in script
    assert "prepare_hprc_real_mini.py" in script
    assert "# HG002 / SRR14724532" in script
