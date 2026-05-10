# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Graph-KIR** is a bioinformatics tool for KIR (Killer Immunoglobulin-like Receptor) gene typing from short-read FASTQ sequencing data. It uses graph-based alignment (via HISAT2) to determine KIR allele types and estimate copy numbers. A companion pipeline `kirpipe` integrates multiple third-party KIR typing tools.

## Development Setup

### WSL (Ubuntu 24.04) — 推薦環境

已驗證可用的環境（2026-05-10）：

| 元件 | 版本 | 路徑 |
|------|------|------|
| Miniconda | latest | `~/miniconda3` |
| conda env | `graphkir_env` | `~/miniconda3/envs/graphkir_env` |
| Python | 3.10.20 | env 內 |
| HISAT2 | 2.2.1 | bioconda |
| samtools | 1.22.1 | bioconda |
| BWA | 0.7.19 | bioconda |
| MUSCLE | 5.3 | bioconda |

```bash
# 啟動 WSL 並進入專案（Windows 路徑掛載於 /mnt/d/）
wsl -d Ubuntu-24.04

# 啟動 conda 環境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate graphkir_env

# 專案路徑
cd /mnt/d/works/KIR_graph
```

**從頭建立環境（若需重建）：**
```bash
# 接受 Anaconda ToS（首次必要）
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda create -n graphkir_env python=3.10 -y
conda activate graphkir_env
conda install -c bioconda -c conda-forge muscle=5.3 hisat2=2.2.1 samtools=1.22.1 bwa=0.7.19 -y
pip install -e .
```

**Docker alternative:**
```bash
docker build -t linnil1/graphkir .
docker run -it --rm -v "$PWD":/data linnil1/graphkir graphkir --help
```

## Commands

**Run graphkir (single sample):**
```bash
graphkir --thread 2 --r1 examples/test00.read1.fq.gz --r2 examples/test00.read2.fq.gz \
    --index-folder example_index --output-folder example_data \
    --output-cohort-name example_data/cohort
```

**Run graphkir (batch via CSV):**
```bash
graphkir --thread 2 --input-csv examples/cohort.csv \
    --index-folder example_index --output-cohort-name example_data/cohort
```

**Type checking:**
```bash
mypy graphkir/ kir/
```

**Format code:**
```bash
black graphkir/ kir/
```

**Build documentation:**
```bash
mkdocs serve   # local preview
mkdocs build   # static site
```

**Tests:** There is no automated test suite. Validation scripts live in `research/`:
- `research/kg_3dx1_test.py` — basic functional test
- `research/kg_main.py` — simulated data validation (100 samples)
- `research/kg_real.py` — real data validation (HPRC cohort)

## Repository Layout

Top-level directories are organized by responsibility:

- `graphkir/` — main installable package and CLI implementation
- `kir/` — wrapper pipeline for Graph-KIR and third-party KIR tools
- `examples/` — example FASTQ inputs and expected outputs for quick local runs
- `data/` — static inputs used by research scripts
  - `data/cohorts/` — cohort manifests such as `hprc.csv`
  - `data/reference/` — static helper inputs such as `KIR_gene_haplotypes.csv` and `hs38.decoy`
  - `data/groundtruth/` — evaluation truth tables such as `hprc_summary_v1_2_e.tsv`
- `research/` — paper workflow, benchmarking, exploratory scripts, and SLURM templates
- `docs/` — MkDocs content and manuscript artifacts

When changing paths in research scripts, prefer keeping static datasets under `data/` and executable logic under `research/`.

## Architecture

### Pipeline stages (`graphkir/`)

The pipeline flows through these sequential stages:

1. **MSA Building** (`kir_msa.py`, `kir_msa_leftalign.py`) — reads the IPD-KIR allele database, constructs and left-aligns Multiple Sequence Alignments. Supports MSA types: `merge`, `split`, `ab`, `ab_2dl1s1`.

2. **Index Creation** (`msa2hisat.py`) — converts MSA to HISAT2 graph index format, tracking variant positions.

3. **Read Mapping** (`hisat2.py`, `wgs.py`) — aligns FASTQ reads to the graph index via HISAT2; `wgs.py` handles extraction from whole-genome sequencing data (hg19/hg38 KIR locus).

4. **Copy Number Estimation** (`kir_cn.py`, `cn_model.py`) — statistical depth-based CN prediction. Supports diploid assumption or full integer CN estimation with cohort-level normalization.

5. **Allele Typing** (`kir_typing.py`, `typing_mulit_allele.py`, `typing_em.py`) — resolves alleles using EM algorithm. Two strategies: standard and `exonfirst` (exon-focused for higher sensitivity). Novel allele discovery in `novel_discover.py`.

6. **Visualization** (`plot.py`) — CN distribution and read mapping statistics using Plotly/Dash.

### Key utilities

- `utils.py` — centralized logger, global resource config (threads/memory), shell command execution
- `external_tools.py` — abstraction layer for external bioinformatics tools with pluggable execution engines (local, Docker, Podman, Singularity)
- `samtools_utils.py` — SAM/BAM file operations

### KIRpipe wrappers (`kir/`)

`kir/kir_pipe.py` defines the base pipeline class. Tool-specific adapters (`graphkir.py`, `ping.py`, `t1k.py`, `sakauekir.py`, `kpi.py`) wrap third-party KIR typing tools into a unified interface.

### Entry points

- `graphkir.main:entrypoint` → CLI command `graphkir`
- `kir.main:main` → CLI command `kirpipe`

## Data Formats

| Format | Purpose |
|--------|---------|
| FASTQ (.fq.gz) | Input sequencing reads |
| CSV | Sample batch input lists |
| BAM/SAM | Internal alignment files |
| `.cn.tsv` | Copy number output per sample |
| `.allele.tsv` | Allele typing output per sample |
| JSON | Complex intermediate data structures |

## Configuration

- `pyproject.toml` — build system, dependencies, mypy and black settings
- `mkdocs.yml` — documentation site configuration
- `Dockerfile` — production container (Python 3.14-slim base)

## Type Checking Notes

mypy is configured in `pyproject.toml`. The following are excluded from type checking:
- Modified Sakaue KIR code
- Scripts in `research/`

Use TypedDicts for complex return types and dataclasses for data structures — this is the existing pattern throughout the codebase.

## External Tool Abstraction

All external bioinformatics tools (HISAT2, samtools, BWA, MUSCLE) are invoked through `external_tools.py`, which supports multiple execution backends. When adding support for a new tool, follow this pattern rather than calling subprocess directly.
