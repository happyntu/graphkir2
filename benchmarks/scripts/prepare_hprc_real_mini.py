"""Prepare a small real HPRC benchmark manifest when local FASTQs are available."""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class HprcSample:
    """One HPRC sample row with a sequencing accession."""

    sample_id: str
    accession: str


@dataclass(frozen=True)
class AvailableSample:
    """One sample with a discovered FASTQ pair."""

    sample: HprcSample
    read1: Path
    read2: Path


def default_fastq_roots() -> tuple[Path, ...]:
    """Return default roots to search for local HPRC FASTQ files."""
    roots: list[str] = []
    env_roots = os.environ.get("GRAPHKIR_HPRC_FASTQ_ROOTS") or os.environ.get(
        "HPRC_FASTQ_ROOT"
    )
    if env_roots:
        roots.extend(field for field in env_roots.split(os.pathsep) if field)
    roots.extend(["data_real", "data/hprc_fastq", "benchmarks/data/hprc_fastq"])
    return tuple(dict.fromkeys(Path(root) for root in roots))


def read_hprc_samples(cohort_csv: Path) -> list[HprcSample]:
    """Read HPRC samples with non-empty short-read accessions."""
    samples: list[HprcSample] = []
    with cohort_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (row.get("id") or "").strip()
            accession = (row.get("sample_id") or "").strip()
            if sample_id and accession:
                samples.append(HprcSample(sample_id=sample_id, accession=accession))
    return samples


def read_truth_sample_ids(truth_tsv: Path) -> frozenset[str]:
    """Read sample IDs present in the HPRC truth table."""
    ids: set[str] = set()
    with truth_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for field in ("id", "name", "sample_id"):
                value = (row.get(field) or "").strip()
                if value:
                    ids.add(value)
    return frozenset(ids)


def candidate_fastq_pairs(
    root: Path, sample: HprcSample
) -> tuple[tuple[Path, Path], ...]:
    """Generate likely FASTQ pair paths for common HPRC naming conventions."""
    accession = sample.accession
    sample_id = sample.sample_id
    pairs = [
        (f"{accession}_1.fastq.gz", f"{accession}_2.fastq.gz"),
        (f"{accession}_1.fq.gz", f"{accession}_2.fq.gz"),
        (f"{accession}_R1.fastq.gz", f"{accession}_R2.fastq.gz"),
        (f"{accession}_R1.fq.gz", f"{accession}_R2.fq.gz"),
        (f"{accession}.read.1.fq.gz", f"{accession}.read.2.fq.gz"),
        (f"hprc.{sample_id}.read.1.fq.gz", f"hprc.{sample_id}.read.2.fq.gz"),
        (f"{sample_id}.read.1.fq.gz", f"{sample_id}.read.2.fq.gz"),
        (f"{sample_id}_1.fastq.gz", f"{sample_id}_2.fastq.gz"),
        (f"{sample_id}_R1.fastq.gz", f"{sample_id}_R2.fastq.gz"),
    ]
    return tuple((root / read1, root / read2) for read1, read2 in pairs)


def _find_named_file(root: Path, filename: str) -> Path | None:
    """Find a file by exact name under a root."""
    for path in root.rglob(filename):
        if path.is_file():
            return path
    return None


def find_fastq_pair(
    roots: tuple[Path, ...],
    sample: HprcSample,
    recursive: bool = False,
) -> tuple[Path, Path] | None:
    """Find a local FASTQ pair for one HPRC sample."""
    for root in roots:
        for read1, read2 in candidate_fastq_pairs(root, sample):
            if read1.exists() and read2.exists():
                return read1, read2
            if not recursive or not root.exists():
                continue
            found_read1 = _find_named_file(root, read1.name)
            found_read2 = _find_named_file(root, read2.name)
            if found_read1 and found_read2:
                return found_read1, found_read2
    return None


def select_available_samples(
    samples: list[HprcSample],
    truth_ids: frozenset[str],
    roots: tuple[Path, ...],
    max_samples: int,
    recursive: bool = False,
) -> tuple[list[AvailableSample], list[HprcSample]]:
    """Select samples that have truth and local paired FASTQs."""
    available: list[AvailableSample] = []
    missing: list[HprcSample] = []
    for sample in samples:
        if sample.sample_id not in truth_ids:
            continue
        pair = find_fastq_pair(roots, sample, recursive=recursive)
        if pair is None:
            missing.append(sample)
            continue
        available.append(AvailableSample(sample=sample, read1=pair[0], read2=pair[1]))
        if len(available) >= max_samples:
            break
    return available, missing


def write_manifest(
    samples: list[AvailableSample], path: Path, output_dir: Path
) -> None:
    """Write a Graph-KIR input manifest for available samples."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "r1", "r2", "cnfile"])
        writer.writeheader()
        for sample in samples:
            writer.writerow(
                {
                    "name": str(output_dir / sample.sample.sample_id),
                    "r1": str(sample.read1),
                    "r2": str(sample.read2),
                    "cnfile": "",
                }
            )


def write_missing_report(samples: list[HprcSample], path: Path) -> None:
    """Write samples that have truth but no discovered local FASTQ pair."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "sample_id"])
        writer.writeheader()
        for sample in samples:
            writer.writerow({"id": sample.sample_id, "sample_id": sample.accession})


def build_preset(
    manifest_path: Path,
    output_cohort_name: Path,
    truth_tsv: Path,
    skip_extraction: bool,
) -> dict[str, object]:
    """Build a mini HPRC discard-baseline preset."""
    return {
        "benchmark_label": "hprc-real-mini",
        "input_csv": str(manifest_path),
        "output_cohort_name": str(output_cohort_name),
        "thread": 4,
        "memory_gb": 14,
        "engine": "local",
        "ref_genome": "hg19",
        "index_folder": "index",
        "ipd_version": "2100",
        "msa_type": "ab_2dl1s1",
        "use_exon_only_alleles": True,
        "multi_map_policy": "discard",
        "allele_strategy": "full",
        "allele_truth_tsv": str(truth_tsv),
        "cn_diploid_gene": "",
        "cn_cohort": False,
        "assume_3dl3_diploid": True,
        "step_skip_extraction": skip_extraction,
    }


def build_enhancedgate_preset(base_preset: dict[str, object]) -> dict[str, object]:
    """Build a targeted enhancedgate typing preset from the baseline preset."""
    preset = dict(base_preset)
    preset.update(
        {
            "benchmark_label": "hprc-real-mini-enhancedgate",
            "allele_exon_weight": 2.0,
            "allele_ambiguity_neutral_prob": 0.999,
            "allele_select_min_fraction_ratio": 0.7,
            "allele_base_top_n": 600,
            "allele_gene_base_top_ns": "KIR2DL1:1000",
            "allele_cross_gene_neutralization_groups": "KIR2DS3/KIR2DS5",
            "allele_private_support_genes": "KIR2DS3",
            "allele_private_support_lambda": 10.0,
            "allele_private_support_window": 50.0,
            "allele_private_support_condition_alleles": "KIR2DS3*00201",
            "allele_private_support_cross_gene_ratio": 0.8,
            "allele_private_support_discard_fallback_genes": "KIR2DS3",
            "allele_private_support_discard_fallback_residual_alleles": "KIR2DS3*00201",
            "allele_private_support_discard_fallback_introduced_alleles": "KIR2DS3*00103",
            "allele_private_support_discard_fallback_introduced_max_ratio": 0.885,
            "allele_private_support_discard_fallback_max_score": -20.0,
            "allele_private_support_discard_fallback_residual_min_ratio": 0.7,
            "allele_functional_discard_fallback_genes": "KIR2DL1,KIR2DS5,KIR2DS3",
            "allele_functional_discard_fallback_resolution": 3,
            "allele_functional_discard_fallback_max_score": -100.0,
            "allele_functional_discard_fallback_min_score_delta": 20.0,
            "allele_functional_discard_fallback_promoted_alleles": "KIR2DS5*027,KIR2DS3*00109",
            "allele_functional_discard_fallback_protected_alleles": "KIR2DS5*002,KIR2DS3*00103",
            "allele_highest_suffix_tie_break_genes": "KIR2DS4",
        }
    )
    return preset


def current_lead_guard_args() -> tuple[str, ...]:
    """Return CLI-only guard flags for the current synthetic functional lead."""
    return (
        "--unsupported-overcall-guard-genes",
        "KIR2DL5",
        "--unsupported-overcall-guard-window",
        "25.0",
        "--unsupported-overcall-guard-min-unsupported-delta",
        "2",
        "--unsupported-overcall-guard-min-net-delta",
        "20.0",
        "--targeted-unsupported-overcall-guard-genes",
        "KIR2DS5",
        "--targeted-unsupported-overcall-guard-alleles",
        "KIR2DS5*027,KIR2DS5*010",
        "--targeted-unsupported-overcall-guard-window",
        "25.0",
        "--targeted-unsupported-overcall-guard-min-unsupported-delta",
        "1",
        "--targeted-unsupported-overcall-guard-min-net-delta",
        "20.0",
        "--targeted-unsupported-overcall-guard-preserve-non-target-resolution",
        "5",
        "--rankwide-unsupported-overcall-guard-genes",
        "KIR2DS3",
        "--rankwide-unsupported-overcall-guard-alleles",
        "KIR2DS3*0020101",
        "--rankwide-unsupported-overcall-guard-window",
        "400.0",
        "--rankwide-unsupported-overcall-guard-min-unsupported-delta",
        "1",
        "--rankwide-unsupported-overcall-guard-min-net-delta",
        "20.0",
        "--rankwide-unsupported-overcall-guard-max-selected-support",
        "-100.0",
        "--rankwide-unsupported-overcall-guard-preserve-non-target-resolution",
        "3",
        "--discard-unsupported-overcall-guard-genes",
        "KIR2DL1",
        "--discard-unsupported-overcall-guard-alleles",
        "KIR2DL1*00303",
        "--discard-unsupported-overcall-guard-window",
        "400.0",
        "--discard-unsupported-overcall-guard-min-unsupported-delta",
        "5",
        "--discard-unsupported-overcall-guard-min-net-delta",
        "100.0",
        "--discard-unsupported-overcall-guard-preserve-non-target-resolution",
        "3",
        "--discard-unsupported-overcall-guard-preserve-selected-resolution",
        "3",
    )


def current_lead_typing_command(
    enhanced_config_path: Path,
    output_tsv: Path,
) -> list[str]:
    """Build the recommended current-lead rerun command."""
    return [
        "python",
        "benchmarks/scripts/rerun_typing_private_support.py",
        "--config",
        str(enhanced_config_path),
        "--top-n",
        "5000",
        *current_lead_guard_args(),
        "--output-tsv",
        str(output_tsv),
    ]


def write_preset(preset: dict[str, object], path: Path) -> None:
    """Write a benchmark preset JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(preset, indent=2) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Prepare a small HPRC real-data benchmark when FASTQs exist locally.",
    )
    parser.add_argument("--cohort-csv", default="data/cohorts/hprc.csv")
    parser.add_argument(
        "--truth-tsv", default="data/groundtruth/hprc_summary_v1_2_e.tsv"
    )
    parser.add_argument("--fastq-root", action="append", default=[])
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--max-samples", type=int, default=4)
    parser.add_argument("--skip-extraction", action="store_true")
    parser.add_argument(
        "--manifest-out",
        default="benchmarks/generated/hprc-real-mini/manifest.csv",
    )
    parser.add_argument(
        "--config-out",
        default="benchmarks/generated/hprc-real-mini/hprc_real_mini.json",
    )
    parser.add_argument(
        "--enhanced-config-out",
        default="benchmarks/generated/hprc-real-mini/hprc_real_mini_enhancedgate.json",
    )
    parser.add_argument(
        "--report-out",
        default="benchmarks/results/hprc-real-mini-availability/missing_fastq.tsv",
    )
    return parser


def main() -> None:
    """Prepare the HPRC real mini benchmark files if local data is available."""
    args = build_parser().parse_args()
    cohort_csv = Path(args.cohort_csv)
    truth_tsv = Path(args.truth_tsv)
    roots = (
        tuple(Path(root) for root in args.fastq_root)
        if args.fastq_root
        else default_fastq_roots()
    )
    samples = read_hprc_samples(cohort_csv)
    truth_ids = read_truth_sample_ids(truth_tsv)
    available, missing = select_available_samples(
        samples,
        truth_ids,
        roots,
        max_samples=args.max_samples,
        recursive=args.recursive,
    )
    write_missing_report(missing, Path(args.report_out))
    if not available:
        print("No eligible HPRC samples with paired local FASTQs were found.")
        print(f"Searched roots: {', '.join(str(root) for root in roots)}")
        print(f"Wrote missing report to {args.report_out}")
        return

    manifest_path = Path(args.manifest_out)
    output_dir = Path("benchmarks/results/hprc-real-mini")
    write_manifest(available, manifest_path, output_dir=output_dir)
    config_path = Path(args.config_out)
    preset = build_preset(
        manifest_path,
        output_dir / "cohort",
        truth_tsv,
        skip_extraction=args.skip_extraction,
    )
    write_preset(preset, config_path)
    enhanced_config_path = Path(args.enhanced_config_out)
    write_preset(build_enhancedgate_preset(preset), enhanced_config_path)
    print(f"Wrote manifest to {manifest_path}")
    print(f"Wrote preset to {config_path}")
    print(f"Wrote enhancedgate preset to {enhanced_config_path}")
    print(f"Wrote missing report to {args.report_out}")
    print("Run legacy baseline:")
    print(
        "python benchmarks/scripts/run_compare.py "
        f"--config {config_path} --execute-legacy --evaluate "
        "--output-dir benchmarks/results/hprc-real-mini-discard-eval"
    )
    print("Run enhancedgate typing from the same intermediates:")
    print(
        " ".join(
            current_lead_typing_command(
                enhanced_config_path,
                Path("benchmarks/results/hprc-real-mini/enhancedgate.allele.tsv"),
            )
        )
    )
    print("Evaluate enhancedgate:")
    print(
        "python benchmarks/scripts/run_compare.py "
        f"--config {enhanced_config_path} --evaluate "
        "--pred-tsv benchmarks/results/hprc-real-mini/enhancedgate.allele.tsv "
        "--multi-map-policy-override likelihood "
        "--output-dir benchmarks/results/hprc-real-mini-enhancedgate-eval"
    )


if __name__ == "__main__":
    main()
