"""Prepare a manual HPRC FASTQ download plan without downloading data."""

from __future__ import annotations

import argparse
import csv
import shlex
from dataclasses import dataclass
from pathlib import Path

from prepare_hprc_real_mini import (
    HprcSample,
    read_hprc_samples,
    read_truth_sample_ids,
)


DEFAULT_COHORT_CSV = "data/cohorts/hprc.csv"
DEFAULT_TRUTH_TSV = "data/groundtruth/hprc_summary_v1_2_e.tsv"
DEFAULT_OUTPUT_DIR = "benchmarks/generated/hprc-download-plan"
DEFAULT_FASTQ_ROOT = "data/hprc_fastq"
DEFAULT_SRA_CACHE = "benchmarks/generated/hprc-download/sra"
DEFAULT_TMP_DIR = "benchmarks/generated/hprc-download/tmp"


@dataclass(frozen=True)
class HprcDownloadRow:
    """One planned HPRC FASTQ download row."""

    sample_id: str
    accession: str
    read1: Path
    read2: Path
    sra_path: Path
    temp_read1: Path
    temp_read2: Path
    download_command: str
    read1_exists: bool
    read2_exists: bool


def select_truth_matched_samples(
    samples: list[HprcSample],
    truth_ids: frozenset[str],
    max_samples: int,
) -> list[HprcSample]:
    """Select cohort samples with truth rows, preserving cohort order."""
    selected = [sample for sample in samples if sample.sample_id in truth_ids]
    if max_samples > 0:
        return selected[:max_samples]
    return selected


def build_download_command(
    accession: str,
    read1: Path,
    read2: Path,
    sra_path: Path,
    temp_read1: Path,
    temp_read2: Path,
    threads: int,
) -> str:
    """Build a manual SRA Toolkit command for one paired FASTQ download."""
    sample_tmp = temp_read1.parent
    cache_root = sra_path.parent.parent
    parts = [
        f"mkdir -p {shlex.quote(str(read1.parent))} {shlex.quote(str(sample_tmp))}",
        f"prefetch {shlex.quote(accession)} --output-directory {shlex.quote(str(cache_root))}",
        (
            f"fasterq-dump {shlex.quote(str(sra_path))} --split-files "
            f"--threads {threads} --outdir {shlex.quote(str(sample_tmp))}"
        ),
        f"gzip -c {shlex.quote(str(temp_read1))} > {shlex.quote(str(read1))}",
        f"gzip -c {shlex.quote(str(temp_read2))} > {shlex.quote(str(read2))}",
        f"rm -f {shlex.quote(str(temp_read1))} {shlex.quote(str(temp_read2))}",
    ]
    return " && ".join(parts)


def build_download_rows(
    samples: list[HprcSample],
    fastq_root: Path,
    sra_cache: Path,
    tmp_dir: Path,
    threads: int,
) -> list[HprcDownloadRow]:
    """Build planned output paths and commands for selected samples."""
    rows: list[HprcDownloadRow] = []
    for sample in samples:
        accession = sample.accession
        read1 = fastq_root / f"{accession}_1.fastq.gz"
        read2 = fastq_root / f"{accession}_2.fastq.gz"
        sra_path = sra_cache / accession / f"{accession}.sra"
        sample_tmp = tmp_dir / accession
        temp_read1 = sample_tmp / f"{accession}_1.fastq"
        temp_read2 = sample_tmp / f"{accession}_2.fastq"
        rows.append(
            HprcDownloadRow(
                sample_id=sample.sample_id,
                accession=accession,
                read1=read1,
                read2=read2,
                sra_path=sra_path,
                temp_read1=temp_read1,
                temp_read2=temp_read2,
                download_command=build_download_command(
                    accession,
                    read1,
                    read2,
                    sra_path,
                    temp_read1,
                    temp_read2,
                    threads,
                ),
                read1_exists=read1.exists(),
                read2_exists=read2.exists(),
            )
        )
    return rows


def write_plan_tsv(rows: list[HprcDownloadRow], path: Path) -> None:
    """Write the download plan as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "id",
                "sample_id",
                "read1",
                "read2",
                "read1_exists",
                "read2_exists",
                "sra_path",
                "download_command",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "id": row.sample_id,
                    "sample_id": row.accession,
                    "read1": str(row.read1),
                    "read2": str(row.read2),
                    "read1_exists": str(row.read1_exists).lower(),
                    "read2_exists": str(row.read2_exists).lower(),
                    "sra_path": str(row.sra_path),
                    "download_command": row.download_command,
                }
            )


def render_download_script(
    rows: list[HprcDownloadRow],
    fastq_root: Path,
) -> str:
    """Render a review-before-running bash download script."""
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Generated helper only. Review disk space before running.",
        "# Requires SRA Toolkit commands: prefetch and fasterq-dump.",
        "# After downloading, validate availability with:",
        (
            "# python benchmarks/scripts/prepare_hprc_real_mini.py "
            f"--fastq-root {shlex.quote(str(fastq_root))} --recursive"
        ),
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"# {row.sample_id} / {row.accession}",
                row.download_command,
                "",
            ]
        )
    return "\n".join(lines)


def write_download_script(
    rows: list[HprcDownloadRow],
    path: Path,
    fastq_root: Path,
) -> None:
    """Write a bash script containing manual download commands."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_download_script(rows, fastq_root), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Create a manual HPRC FASTQ download plan from the cohort and truth "
            "tables. This does not download data."
        )
    )
    parser.add_argument("--cohort-csv", default=DEFAULT_COHORT_CSV)
    parser.add_argument("--truth-tsv", default=DEFAULT_TRUTH_TSV)
    parser.add_argument("--fastq-root", default=DEFAULT_FASTQ_ROOT)
    parser.add_argument("--sra-cache", default=DEFAULT_SRA_CACHE)
    parser.add_argument("--tmp-dir", default=DEFAULT_TMP_DIR)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-samples", type=int, default=4)
    parser.add_argument("--threads", type=int, default=8)
    return parser


def main() -> None:
    """Write the HPRC download plan files."""
    args = build_parser().parse_args()
    cohort_csv = Path(args.cohort_csv)
    truth_tsv = Path(args.truth_tsv)
    fastq_root = Path(args.fastq_root)
    output_dir = Path(args.output_dir)
    samples = select_truth_matched_samples(
        read_hprc_samples(cohort_csv),
        read_truth_sample_ids(truth_tsv),
        args.max_samples,
    )
    rows = build_download_rows(
        samples,
        fastq_root=fastq_root,
        sra_cache=Path(args.sra_cache),
        tmp_dir=Path(args.tmp_dir),
        threads=args.threads,
    )
    plan_tsv = output_dir / "hprc_fastq_download_plan.tsv"
    script_path = output_dir / "download_hprc_fastq.sh"
    write_plan_tsv(rows, plan_tsv)
    write_download_script(rows, script_path, fastq_root)
    print(f"Wrote HPRC FASTQ download plan to {plan_tsv}")
    print(f"Wrote manual download script to {script_path}")
    print(f"Planned samples: {len(rows)}")
    print("This helper did not download data.")


if __name__ == "__main__":
    main()
