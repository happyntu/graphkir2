# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Graph-KIR** is a bioinformatics tool for KIR (Killer Immunoglobulin-like Receptor) gene typing from short-read FASTQ sequencing data. It uses graph-based alignment (via HISAT2) to determine KIR allele types and estimate copy numbers. A companion pipeline `kirpipe` integrates multiple third-party KIR typing tools.

## GraphKir2 Refactor Objective

The current `graphkir2` refactor is not only a code cleanup. It is a method-improvement
project with a specific benchmark target:

* the primary goal is to make `graphkir2` stronger than `Geny` on `3-digit` and
  `5-digit` functional typing
* the refactor should preserve or improve Graph-KIR's current strengths in:
  * copy number estimation
  * `7-digit` full-resolution typing
  * runtime efficiency

When tradeoffs appear, do not optimize only for `7-digit` or only for architectural
cleanliness. Prefer changes that improve practical short-read functional typing while
keeping the current CN and speed advantages intact.

## Source of Truth for Refactor Priorities

When deciding what to implement next for `graphkir2`, use this priority:

1. benchmark goal: improve `3-digit/5-digit` functional typing against `Geny`
2. preserve Graph-KIR strengths in CN and `7-digit`
3. keep legacy `graphkir/` runnable as the comparison baseline
4. keep manifests, benchmark outputs, and research notes reproducible

Do not treat the current paper framing alone as the final optimization target. The
working development target is stronger functional typing performance, not only defending
the existing `7-digit` and CN story.

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
- `src/graphkir2/` — refactor area for the next-generation implementation
- `examples/` — example FASTQ inputs and expected outputs for quick local runs
- `data/` — static inputs used by research scripts
  - `data/cohorts/` — cohort manifests such as `hprc.csv`
  - `data/reference/` — static helper inputs such as `KIR_gene_haplotypes.csv` and `hs38.decoy`
  - `data/groundtruth/` — evaluation truth tables such as `hprc_summary_v1_2_e.tsv`
- `research/` — paper workflow, benchmarking, exploratory scripts, and SLURM templates
- `benchmarks/` — old-vs-new comparison harness, configs, and result summaries
- `docs/` — MkDocs content and manuscript artifacts

When changing paths in research scripts, prefer keeping static datasets under `data/` and executable logic under `research/`.
When refactoring the implementation, keep the current `graphkir/` package stable as the baseline and put new code under `src/graphkir2/`.

For `graphkir2`, benchmark-oriented design matters as much as implementation:

* architecture and stage contracts belong in `docs/architecture.md`
* benchmark scope and metrics belong in `docs/benchmark_plan.md`
* staged delivery belongs in `docs/roadmap.md`
* accepted durable design decisions belong in `docs/adr/`
* exploratory tuning and failure analysis belong in `docs/research/`

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
- `graphkir.main2:entrypoint` → CLI command `graphkir2`

## Benchmark Priorities

For the legacy paper, Graph-KIR already has a strong story in copy number and `7-digit`
resolution. For the refactor, the benchmark priority is different:

* first priority: beat `Geny` on `3-digit` and `5-digit` functional typing
* second priority: do not regress on CN
* third priority: do not regress on `7-digit`
* fourth priority: preserve practical speed advantages

This means method work should be evaluated not only by aggregate accuracy, but also by:

* `3-digit`
* `5-digit`
* `7-digit`
* copy number
* per-gene failure modes
* runtime and memory

Useful likely targets for improvement include:

* multi-mapped read handling
* allele likelihood redesign
* structural / novel case refinement
* confidence-aware downgrade behavior

Current synthetic functional-target surface:

* `multi_map_policy = likelihood`
* `allele_exon_weight = 2.0`
* `allele_select_min_fraction_ratio = 0.7`
* `top_n = 5000` for typing sweeps
* targeted `KIR2DS3` private-support reranking
* directional `KIR2DS3/KIR2DS5` cross-gene ambiguity neutralization applied to
  the `KIR2DS3` target evidence, not to `KIR2DS5` typing evidence

This surface is represented by
`benchmarks/configs/synthetic-difficult5-functional-target-kir2ds3-private.json`.
It is the current lead for the `3-digit/5-digit` objective. The directional
neutralization variant reaches `synthetic-difficult5` `1.0 / 1.0 / 0.95`,
leaving `KIR2DS4` 7-digit resolution as the main remaining difficult5 tradeoff.

## Synthetic-First Workflow

For `graphkir2`, prefer a synthetic-first benchmark and tuning loop:

1. use synthetic data to find runtime bottlenecks
2. use synthetic data to tune `3-digit/5-digit` functional typing behavior
3. use synthetic difficult cases for multi-map / structural / novel ablation
4. only then run HPRC or other real-data sanity checks

This rule exists because:

* synthetic data is cheaper for repeated profiling and parameter sweeps
* synthetic truth is explicit for `3-digit`, `5-digit`, and `7-digit`
* real-data runs should confirm candidate improvements, not serve as the first search loop

When deciding between two next steps, prefer:

* synthetic profiling or ablation first

over:

* immediate broad real-data benchmarking

unless the task is explicitly about final validation or paper-ready real-data comparison.

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
