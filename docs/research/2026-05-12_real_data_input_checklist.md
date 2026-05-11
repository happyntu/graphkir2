# Real-Data Input Checklist

This checklist defines the minimum local inputs needed before `graphkir2` can make
a real-data claim against `Geny`. It should be used after the synthetic candidate
has stabilized and before interpreting any HPRC or Geny comparison output.

Current synthetic lead:

```text
enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware
```

The real-data gate is blocked until both Geny output and HPRC paired FASTQs are
available locally. Do not tune the current KIR2DS3 residual synthetic row as a
substitute for this gate.

## Required Inputs

| Input | Default path | Required for | Current action |
|---|---|---|---|
| HPRC truth TSV | `data/groundtruth/hprc_summary_v1_2_e.tsv` | Geny and HPRC evaluation | Already expected in repo data |
| Raw Geny output | `data/geny_hprc44.txt` | Direct Geny functional comparison | Provide locally |
| HPRC paired FASTQs | `data_real`, `data/hprc_fastq`, `benchmarks/data/hprc_fastq`, or `--fastq-root` | HPRC real-mini execution | Provide locally |
| Optional graphkir2 allele TSV | e.g. `benchmarks/results/hprc-real-mini/enhancedgate.allele.tsv` | Side-by-side Geny vs graphkir2 comparison | Generated after HPRC mini run |

Large raw data and generated benchmark outputs are intentionally not committed.
Keep raw FASTQs and Geny outputs local or in external storage.

## Accepted HPRC FASTQ Names

`benchmarks/scripts/prepare_hprc_real_mini.py` searches each root for common
paired-read names using either the sequencing accession from `data/cohorts/hprc.csv`
or the HPRC sample ID:

```text
<accession>_1.fastq.gz / <accession>_2.fastq.gz
<accession>_1.fq.gz / <accession>_2.fq.gz
<accession>_R1.fastq.gz / <accession>_R2.fastq.gz
<accession>_R1.fq.gz / <accession>_R2.fq.gz
<accession>.read.1.fq.gz / <accession>.read.2.fq.gz
hprc.<sample>.read.1.fq.gz / hprc.<sample>.read.2.fq.gz
<sample>.read.1.fq.gz / <sample>.read.2.fq.gz
<sample>_1.fastq.gz / <sample>_2.fastq.gz
<sample>_R1.fastq.gz / <sample>_R2.fastq.gz
```

If files are nested below a root, add `--recursive`. Prefer `--fastq-root` for
explicit one-off runs. Use `GRAPHKIR_HPRC_FASTQ_ROOTS` or `HPRC_FASTQ_ROOT` only
when the same root should be reused repeatedly.

## Download Plan Helper

If HPRC FASTQs are not available locally, generate a reviewable download plan
instead of starting a large download automatically:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/prepare_hprc_download_plan.py"
```

This writes:

```text
benchmarks/generated/hprc-download-plan/hprc_fastq_download_plan.tsv
benchmarks/generated/hprc-download-plan/download_hprc_fastq.sh
```

The helper reads `data/cohorts/hprc.csv`, keeps only samples present in the HPRC
truth TSV, and defaults to the first four samples for a mini sanity run. It
prints SRA Toolkit `prefetch` / `fasterq-dump` commands but does not execute
them. The generated `prefetch` commands use `--max-size 100G` because some HPRC
accessions are larger than the SRA Toolkit 20G default. Increase scope only
after checking disk and runtime budget:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/prepare_hprc_download_plan.py --max-samples 0 --threads 8"
```

## Readiness Checks

Run from Windows PowerShell through the verified WSL environment:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/compare_geny_functional.py"
```

Expected blocked-state output:

```text
benchmarks/results/geny-functional-comparison/missing_inputs.tsv
```

The comparison is ready only when `missing_inputs.tsv` has `exists=true` for both
`truth_tsv` and `geny_output`.

Check HPRC FASTQ availability:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/prepare_hprc_real_mini.py"
```

Expected blocked-state output:

```text
benchmarks/results/hprc-real-mini-availability/missing_fastq.tsv
```

If FASTQs live elsewhere, pass explicit roots:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/prepare_hprc_real_mini.py --fastq-root /mnt/d/path/to/hprc_fastq --recursive"
```

The HPRC mini gate is ready only when the script writes:

```text
benchmarks/generated/hprc-real-mini/manifest.csv
benchmarks/generated/hprc-real-mini/hprc_real_mini.json
benchmarks/generated/hprc-real-mini/hprc_real_mini_enhancedgate.json
```

## Execution Order After Inputs Exist

1. Run `prepare_hprc_real_mini.py` to create the manifest and presets.
2. Run the printed legacy baseline command with `run_compare.py --execute-legacy --evaluate`.
3. Run the printed current-lead `rerun_typing_private_support.py` command to generate `enhancedgate.allele.tsv`.
4. Run the printed enhancedgate evaluation command.
5. Run Geny functional comparison with the graphkir2 allele TSV:

```powershell
wsl -d Ubuntu-24.04 bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate graphkir_env && cd /mnt/d/works/KIR_graph && python benchmarks/scripts/compare_geny_functional.py --geny-output data/geny_hprc44.txt --graphkir-tsv benchmarks/results/hprc-real-mini/enhancedgate.allele.tsv"
```

## Interpretation Rules

Do not claim `graphkir2` beats `Geny` unless:

* `summary.tsv` contains evaluated rows for both `Geny` and `graphkir2`
* `matched_samples` is non-zero and sample IDs match the HPRC truth table
* the claim is stated separately for `3-digit`, `5-digit`, and `7-digit`
* per-gene failures are inspected for KIR2DL1, KIR2DS3, KIR2DS4, KIR2DS5, and KIR2DL5
* runtime and memory are recorded separately from accuracy

If `graphkir2` improves 3/5-digit but regresses CN, 7-digit, or runtime, keep the
result as a candidate finding rather than a final benchmark claim.
