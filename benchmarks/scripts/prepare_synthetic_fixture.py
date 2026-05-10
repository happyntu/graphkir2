"""Generate a small synthetic benchmark fixture for graphkir2."""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path
from typing import Any, cast

from Bio import SeqIO
from pyhlamsa import KIRmsa


ART_URL = (
    "https://www.niehs.nih.gov/research/resources/assets/docs/"
    "artbinmountrainier2016.06.05linux64.tgz"
)


def build_parser() -> argparse.ArgumentParser:
    """Build the synthetic fixture generator parser."""
    parser = argparse.ArgumentParser(
        description="Generate a small synthetic Graph-KIR benchmark fixture.",
    )
    parser.add_argument("--label", default="synthetic-smoke", help="Fixture label.")
    parser.add_argument(
        "--output-root",
        default="benchmarks/generated",
        help="Output root for generated synthetic data.",
    )
    parser.add_argument(
        "--config-out",
        help="Optional preset JSON output path. Defaults to benchmarks/configs/<label>.json",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=2,
        help="Number of synthetic samples to generate.",
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=10,
        help="Read depth passed to art_illumina.",
    )
    parser.add_argument(
        "--allele-seed",
        type=int,
        default=2022,
        help="Random seed for allele selection.",
    )
    parser.add_argument(
        "--read-seed",
        type=int,
        default=1031,
        help="Base random seed for read simulation.",
    )
    parser.add_argument(
        "--msa-type",
        default="ab_2dl1s1",
        choices=["merge", "split", "ab", "ab_2dl1s1"],
        help="MSA type to record in the generated preset.",
    )
    parser.add_argument(
        "--msa-no-exon-only-allele",
        action="store_true",
        help="Record a no-exon-only-alleles policy in the generated preset.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing fixture files.",
    )
    parser.add_argument(
        "--gene-whitelist",
        help=(
            "Optional comma-separated gene panel. When set, synthetic alleles are "
            "drawn directly from these genes instead of full haplotypes."
        ),
    )
    return parser


def ensure_art_illumina(repo_root: Path) -> Path:
    """Download and unpack art_illumina when it is missing."""
    art_path = repo_root / "art_illumina"
    if art_path.exists():
        return art_path

    archive_path = repo_root / "artbinmountrainier2016.06.05linux64.tgz"
    try:
        subprocess.run(
            ["wget", "-O", str(archive_path), ART_URL],
            check=True,
            cwd=repo_root,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        request = urllib.request.Request(
            ART_URL,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(request) as response, archive_path.open("wb") as handle:
            handle.write(response.read())
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(repo_root)
    extracted = repo_root / "art_bin_MountRainier" / "art_illumina"
    extracted.rename(art_path)
    shutil.rmtree(extracted.parent)
    archive_path.unlink()
    art_path.chmod(art_path.stat().st_mode | 0o111)
    return art_path


def read_haplotypes(repo_root: Path) -> list[dict[str, str]]:
    """Read haplotype definitions from the repository data folder."""
    import pandas as pd

    haplo_path = repo_root / "data" / "reference" / "KIR_gene_haplotypes.csv"
    df = pd.read_csv(haplo_path)
    rows = cast(list[dict[str, Any]], df.to_dict(orient="records"))
    return [{str(key): str(value) for key, value in row.items()} for row in rows]


def load_full_length_gene_sequences(version: str = "2100") -> dict[str, list[Any]]:
    """Load full-length KIR sequences grouped by gene."""
    kir = KIRmsa(filetype=["gen"], version=version)
    grouped: dict[str, list[Any]] = {}
    for _, msa in kir.items():
        for record in msa.to_records(gap=False):
            sequence = record.seq
            if sequence is None or len(sequence) <= 4200:
                continue
            allele_id = str(record.id)
            grouped.setdefault(allele_id[:7], []).append(record)
    return grouped


def random_two_haplotypes(
    haplotypes: list[dict[str, str]],
) -> tuple[list[str], dict[str, int]]:
    """Randomly select two haplotypes and compute gene copy counts."""
    selected = random.choices(haplotypes, k=2)
    haplo_ids = [str(row["hapID"]) for row in selected]
    gene_count: dict[str, int] = {}
    for row in selected:
        for key, value in row.items():
            if not str(key).startswith("KIR"):
                continue
            gene_count[key] = gene_count.get(key, 0) + int(value)
    return haplo_ids, gene_count


def random_select_alleles(
    haplotypes: list[dict[str, str]],
    gene_sequences: dict[str, list[Any]],
) -> tuple[list[str], list[Any]]:
    """Randomly select allele sequences according to haplotype copy counts."""
    from copy import deepcopy

    haplo_ids, gene_count = random_two_haplotypes(haplotypes)
    alleles: list[Any] = []
    for gene, count in gene_count.items():
        if gene not in gene_sequences or count <= 0:
            continue
        alleles.extend(random.choices(gene_sequences[gene], k=count))
    allele_name_count: dict[str, int] = {}
    copied = [deepcopy(record) for record in alleles]
    for record in copied:
        allele_name_count[record.id] = allele_name_count.get(record.id, 0) + 1
        record.id = f"{record.id}-{allele_name_count[record.id]}"
        record.description = ""
    return haplo_ids, copied


def select_alleles_for_gene_panel(
    genes: list[str],
    gene_sequences: dict[str, list[Any]],
) -> tuple[list[str], list[Any]]:
    """Select two alleles per requested gene for a simpler functional typing panel."""
    from copy import deepcopy

    copied: list[Any] = []
    allele_name_count: dict[str, int] = {}
    for gene in genes:
        if gene not in gene_sequences:
            raise ValueError(f"Gene {gene} is not available in the loaded KIR sequences")
        selected = random.choices(gene_sequences[gene], k=2)
        copied.extend(deepcopy(record) for record in selected)
    for record in copied:
        allele_name_count[record.id] = allele_name_count.get(record.id, 0) + 1
        record.id = f"{record.id}-{allele_name_count[record.id]}"
        record.description = ""
    return ["manual_panel"], copied


def rewrite_sam(file_in: Path, file_out: Path) -> None:
    """Normalize allele suffixes in ART SAM output."""
    import re

    seen: set[str] = set()
    with file_in.open("r", encoding="utf-8") as handle_in, file_out.open(
        "w", encoding="utf-8"
    ) as handle_out:
        for line in handle_in:
            line = re.sub(r"KIR(.*?)\*(\w+?)\-\d\t", r"KIR\1*\2\t", line)
            if line in seen:
                continue
            seen.add(line)
            handle_out.write(line)


def generate_fastq(
    art_path: Path,
    input_fasta: Path,
    output_prefix: Path,
    depth: int,
    seed: int,
) -> None:
    """Generate paired-end reads with ART."""
    cmd = [
        str(art_path),
        "-ss",
        "HS25",
        "-i",
        str(input_fasta),
        "-l",
        "150",
        "-f",
        str(depth),
        "-m",
        "400",
        "-s",
        "10",
        "-sam",
        "-na",
        "-rs",
        str(seed),
        "-o",
        str(output_prefix) + ".read.",
    ]
    subprocess.run(cmd, check=True, cwd=art_path.parent)


def write_truth_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write synthetic truth in the same TSV shape used by the evaluator."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["id", "name", "haplos", "alleles"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a Graph-KIR style sample manifest CSV."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "r1", "r2", "cnfile"])
        writer.writeheader()
        writer.writerows(rows)


def allele_list_to_cn_rows(alleles: list[str]) -> list[dict[str, str]]:
    """Convert a truth allele list into a simple gene-CN TSV row set."""
    counts: dict[str, int] = {}
    for allele in alleles:
        gene = allele.split("*", 1)[0] + "*BACKBONE"
        counts[gene] = counts.get(gene, 0) + 1
    return [
        {"gene": gene, "cn": str(count)}
        for gene, count in sorted(counts.items())
    ]


def write_cn_truth(path: Path, alleles: list[str]) -> None:
    """Write a per-sample truth CN TSV in Graph-KIR-compatible shape."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "cn"], delimiter="\t")
        writer.writeheader()
        writer.writerows(allele_list_to_cn_rows(alleles))


def write_preset(
    path: Path,
    label: str,
    manifest_path: str,
    truth_path: str,
    index_folder: str,
    msa_type: str,
    use_exon_only_alleles: bool,
) -> None:
    """Write a benchmark preset for the generated fixture."""
    preset = {
        "benchmark_label": label,
        "input_csv": manifest_path,
        "output_cohort_name": f"benchmarks/results/{label}/cohort",
        "thread": 2,
        "memory_gb": 7,
        "engine": "local",
        "ref_genome": "hg19",
        "index_folder": index_folder,
        "ipd_version": "2100",
        "msa_type": msa_type,
        "use_exon_only_alleles": use_exon_only_alleles,
        "multi_map_policy": "discard",
        "allele_strategy": "full",
        "allele_exon_weight": 1.0,
        "allele_truth_tsv": truth_path,
        "cn_diploid_gene": "",
        "cn_cohort": False,
        "assume_3dl3_diploid": True,
        "step_skip_extraction": True,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(preset, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    """Generate a small synthetic benchmark fixture."""
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    output_root = repo_root / args.output_root
    dataset_dir = output_root / args.label
    dataset_dir.mkdir(parents=True, exist_ok=True)

    config_out = (
        repo_root / args.config_out
        if args.config_out
        else repo_root / "benchmarks" / "configs" / f"{args.label}.json"
    )
    truth_path = dataset_dir / f"{args.label}_truth.tsv"
    manifest_path = dataset_dir / "manifest.csv"

    art_path = ensure_art_illumina(repo_root)
    haplotypes = read_haplotypes(repo_root)
    gene_sequences = load_full_length_gene_sequences()
    random.seed(args.allele_seed)
    gene_panel = [
        gene.strip()
        for gene in (args.gene_whitelist or "").split(",")
        if gene.strip()
    ]

    truth_rows: list[dict[str, str]] = []
    manifest_rows: list[dict[str, str]] = []

    for sample_index in range(args.samples):
        sample_id = f"{sample_index:02d}"
        sample_prefix = dataset_dir / f"{args.label}.{sample_id}"
        fasta_path = Path(str(sample_prefix) + ".fa")
        if args.force or not fasta_path.exists():
            if gene_panel:
                haplo_ids, alleles = select_alleles_for_gene_panel(gene_panel, gene_sequences)
            else:
                haplo_ids, alleles = random_select_alleles(haplotypes, gene_sequences)
            SeqIO.write(alleles, fasta_path, "fasta")
            truth_alleles = [record.id.split("-")[0] for record in alleles]
            cn_truth_path = Path(str(sample_prefix) + ".truth.cn.tsv")
            write_cn_truth(cn_truth_path, truth_alleles)
            truth_rows.append(
                {
                    "id": sample_id,
                    "name": str(sample_prefix.relative_to(repo_root)).replace("\\", "/"),
                    "haplos": "_".join(haplo_ids),
                    "alleles": "_".join(truth_alleles),
                }
            )
        else:
            raise SystemExit(
                f"{fasta_path} already exists. Re-run with --force or choose another label."
            )

        read_seed = args.read_seed + sample_index
        if args.force or not (sample_prefix.parent / f"{sample_prefix.name}.read.1.fq").exists():
            generate_fastq(
                art_path=art_path,
                input_fasta=fasta_path,
                output_prefix=sample_prefix,
                depth=args.depth,
                seed=read_seed,
            )
            sam_in = sample_prefix.parent / f"{sample_prefix.name}.read..sam"
            sam_out = Path(str(sample_prefix) + ".sam")
            rewrite_sam(sam_in, sam_out)
        manifest_rows.append(
            {
                "name": str(sample_prefix.relative_to(repo_root)).replace("\\", "/"),
                "r1": str(
                    (sample_prefix.parent / f"{sample_prefix.name}.read.1.fq").relative_to(repo_root)
                ).replace("\\", "/"),
                "r2": str(
                    (sample_prefix.parent / f"{sample_prefix.name}.read.2.fq").relative_to(repo_root)
                ).replace("\\", "/"),
                "cnfile": str(
                    Path(str(sample_prefix) + ".truth.cn.tsv").relative_to(repo_root)
                ).replace("\\", "/"),
            }
        )

    write_truth_tsv(truth_path, truth_rows)
    write_manifest(manifest_path, manifest_rows)
    write_preset(
        config_out,
        label=args.label,
        manifest_path=str(manifest_path.relative_to(repo_root)).replace("\\", "/"),
        truth_path=str(truth_path.relative_to(repo_root)).replace("\\", "/"),
        index_folder=str((dataset_dir / f"index_{args.msa_type}").relative_to(repo_root)).replace("\\", "/"),
        msa_type=args.msa_type,
        use_exon_only_alleles=not args.msa_no_exon_only_allele,
    )
    print(f"Synthetic fixture ready: {dataset_dir}")
    print(f"Manifest: {manifest_path.relative_to(repo_root)}")
    print(f"Truth:    {truth_path.relative_to(repo_root)}")
    print(f"Preset:   {config_out.relative_to(repo_root)}")


if __name__ == "__main__":
    main()
