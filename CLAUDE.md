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
* `allele_base_top_n = 0` for reduced synthetic panels, which keeps all genes at
  the high `top_n`
* `allele_base_top_n = 600` plus `allele_gene_base_top_ns = KIR2DL1:1000` for
  full-gene smoke or real-data sanity reruns, so only target genes keep the high
  `top_n` and known-sensitive non-target genes retain enough search depth
* targeted `KIR2DS3` private-support reranking
* directional `KIR2DS3/KIR2DS5` cross-gene ambiguity neutralization applied to
  the `KIR2DS3` target evidence, not to `KIR2DS5` typing evidence
* `KIR2DS4` same-5-digit exact-likelihood tie-break that keeps the highest
  full allele suffix

This surface is represented by
`benchmarks/configs/synthetic-difficult5-functional-target-kir2ds3-private.json`.
It is the current lead for the original 4-sample `synthetic-difficult5` panel:
`1.0 / 1.0 / 1.0` without regressing `synthetic-functional8` or
`synthetic-functional8x6`.

Do not promote it as the global default yet. The larger
`synthetic-difficult5x12` overfit check shows that the likelihood/private-support
surface still regresses `KIR2DS3` versus the discard baseline:

* discard baseline: `0.9667 / 0.9583 / 0.925`
* balanced likelihood top5000: `0.925 / 0.925 / 0.9083`
* functional target: `0.9417 / 0.9417 / 0.9167`

`graphkir2` now exposes an experimental conditional rescue gate:

* `allele_private_support_condition_alleles = KIR2DS3*00201`
* `allele_private_support_cross_gene_ratio = 0.8`

This first formalized conditional gate starts from the balanced likelihood call,
checks sample-level `KIR2DS3/KIR2DS5` cross-gene support, then reruns directional
neutralized private-support rescue only when the gate passes. It is represented
by `benchmarks/configs/synthetic-difficult5-conditional-kir2ds3-ratio080.json`
and
`benchmarks/configs/synthetic-difficult5x12-conditional-kir2ds3-ratio080.json`.

The formalized gate keeps `synthetic-difficult5` at `1.0 / 1.0 / 1.0` and
improves `synthetic-difficult5x12` to `0.9583 / 0.9583 / 0.9333`, but still does
not recover discard's `KIR2DS3` robustness. Treat this as the current targeted
experimental candidate, not as the final default.

`graphkir2` now also exposes a formalized second-gate fallback that recomputes
the target gene with `removeMultipleMapped` discard-style evidence inside the
same rerun pipeline. It falls back for residual `KIR2DS3*00201` or low-ratio
introduced-`00103` cases and improves `synthetic-difficult5x12` to
`0.975 / 0.975 / 0.950` while keeping:

* `synthetic-difficult5`: `1.0 / 1.0 / 1.0`
* `synthetic-functional8`: `1.0 / 1.0 / 0.90625`
* `synthetic-functional8x6`: `0.9896 / 0.9896 / 0.90625`

The second gate is represented by
`benchmarks/configs/*-conditional-kir2ds3-secondgate.json`. Keep it as the
current targeted synthetic candidate, not the global default, until it is tested
on broader synthetic seeds and real-data sanity panels.

Broader synthetic seed validation on `synthetic-difficult5x12-seed5101`,
`seed5102`, and `seed5103` showed the original secondgate was not robust enough:

* discard mean: `0.9583 / 0.9528 / 0.9000`
* ratio080 mean: `0.9667 / 0.9611 / 0.9361`
* secondgate mean: `0.9639 / 0.9583 / 0.9361`

However, `KIR2DS3` functional accuracy regresses versus discard:

* discard `KIR2DS3`: `0.9583 / 0.9306 / 0.9306`
* ratio080 `KIR2DS3`: `0.9306 / 0.9028 / 0.8750`
* secondgate `KIR2DS3`: `0.9167 / 0.8889 / 0.8750`

Do not move this method to real-data sanity yet. Next method work must recover
`KIR2DS3` 3/5-digit robustness across seeds while preserving the whole-panel
`7-digit` gain.

The current improved candidate is `enhancedgate`, represented by
`benchmarks/configs/*-conditional-kir2ds3-enhancedgate.json`. It adds two
protections on top of secondgate:

* fallback requires base private-support score `<= -20.0`
* residual `KIR2DS3*00201` fallback requires cross-support ratio `>= 0.70`,
  unless the base call already contains the introduced-risk `KIR2DS3*00103`

Across the original `synthetic-difficult5x12` plus `seed5101/5102/5103`, the
current means are:

* discard mean: `0.9604 / 0.9542 / 0.9062`
* ratio080 mean: `0.9646 / 0.9604 / 0.9354`
* secondgate mean: `0.9667 / 0.9625 / 0.9396`
* enhancedgate mean: `0.9792 / 0.9750 / 0.9500`

`KIR2DS3` means over the same four panels are:

* discard `KIR2DS3`: `0.9583 / 0.9271 / 0.9271`
* enhancedgate `KIR2DS3`: `0.9896 / 0.9688 / 0.9375`

Treat enhancedgate as the current aggregate synthetic lead, not as a final
default. A broader functional stress sweep with
`benchmarks/scripts/run_functional_stress_sweep.py` over `synthetic-functional8`,
`synthetic-functional8x6`, `synthetic-difficult5`, and the four difficult5x12
seed panels showed:

* discard mean: `0.9693 / 0.9643 / 0.8854`
* likelihood-only top5000 mean: `0.9437 / 0.9426 / 0.8946`
* enhancedgate gene-aware top-n mean: `0.9866 / 0.9842 / 0.9446`
* enhancedgate + `KIR2DL1` functional fallback gene-aware mean:
  `0.9881 / 0.9842 / 0.9461`
* enhancedgate + `KIR2DL1` fallback + `KIR2DS5` promotion guard gene-aware mean:
  `0.9905 / 0.9866 / 0.9473`
* enhancedgate functional guard gene-aware mean:
  `0.9905 / 0.9878 / 0.9485`

The `KIR2DL1` fallback candidate enables:

* `allele_functional_discard_fallback_genes = KIR2DL1`
* `allele_functional_discard_fallback_resolution = 3`
* `allele_functional_discard_fallback_max_score = -100.0`
* `allele_functional_discard_fallback_min_score_delta = 20.0`

This fixes the strict `synthetic-functional8x6` `KIR2DL1` 3-digit regression
versus discard (`0.9167 -> 1.0000`) while preserving the `KIR2DS3`, `KIR2DS4`,
and `KIR2DS5` gains on difficult panels. KIR2DL1 still has one 5-digit miss
matching discard's remaining error, so the fallback is a functional regression
guard rather than a full 5-digit KIR2DL1 solution. Do not promote it as a global
default until it has real-data sanity coverage.

The current synthetic lead is
`enhancedgate_kir2dl5_kir2ds5unsupported_geneaware`. It keeps the `KIR2DL1`
functional fallback, narrow promotion guards for `KIR2DS5` and `KIR2DS3`, the
KIR2DL5-only unsupported-overcall guard, and adds a targeted KIR2DS5 selected
allele unsupported-overcall guard:

* `allele_functional_discard_fallback_genes = KIR2DL1,KIR2DS5,KIR2DS3`
* `allele_functional_discard_fallback_promoted_alleles = KIR2DS5*027,KIR2DS3*00109`
* `allele_functional_discard_fallback_protected_alleles = KIR2DS5*002,KIR2DS3*00103`
* `unsupported_overcall_guard_genes = KIR2DL5`
* `unsupported_overcall_guard_window = 25.0`
* `unsupported_overcall_guard_min_unsupported_delta = 2`
* `unsupported_overcall_guard_min_net_delta = 20.0`
* `targeted_unsupported_overcall_guard_genes = KIR2DS5`
* `targeted_unsupported_overcall_guard_alleles = KIR2DS5*027,KIR2DS5*010`
* `targeted_unsupported_overcall_guard_min_unsupported_delta = 1`
* `targeted_unsupported_overcall_guard_min_net_delta = 20.0`
* `targeted_unsupported_overcall_guard_preserve_non_target_resolution = 5`

The promotion guard fires only when likelihood introduces a configured promoted
allele and discard has extra protected copies. When possible, it reuses an
already selected protected full allele suffix instead of blindly copying
discard's 7-digit suffix. The KIR2DL5 guard is separate: it compares the
likelihood-selected KIR2DL5 genotype against nearby top-n alternatives and only
switches when the selected genotype has at least two more unsupported
selected-only variants and at least `20.0` more negative-minus-positive support.
Do not loosen this to one unsupported variant; that created a KIR2DL5B
candidate regression during tuning. The KIR2DS5 targeted guard is separate: it
only scores selected candidate-only variants carried by configured KIR2DS5
target prefixes, and alternatives must preserve the non-target selected
5-digit class. This prevents the KIR2DS5*027 fix from shifting to a different
wrong genotype such as KIR2DS5*010 while still allowing the truth-like
non-target allele to remain fixed.

In the functional stress sweep, the previous functional guard fixed the seed5101
KIR2DS5 candidate regressions and the seed5103 KIR2DS3 5-digit suballele
regression. The KIR2DL5 overcall guard then fixes the remaining KIR2DL5
functional misses, and the KIR2DS5 targeted unsupported guard removes the
remaining KIR2DS5 3/5-digit rows:

* `synthetic-difficult5x12`: `0.9750 / 0.9750 / 0.9500` -> `0.9833 / 0.9833 / 0.9500`
* `synthetic-difficult5x12-seed5101`: `0.9917 / 0.9917 / 0.9750` -> `1.0000 / 1.0000 / 0.9833`
* `synthetic-difficult5x12-seed5103`: `0.9917 / 0.9917 / 0.9500` -> `1.0000 / 1.0000 / 0.9583`
* KIR2DS5 targeted guard on `synthetic-difficult5x12`: `0.9833 / 0.9833 / 0.9500` -> `0.9917 / 0.9917 / 0.9583`
* KIR2DS5 targeted guard on `synthetic-difficult5x12-seed5102`: `0.9750 / 0.9667 / 0.9417` -> `1.0000 / 0.9917 / 0.9750`

The current 7-panel stress mean is:

* `enhancedgate_functionalguard_geneaware`: `0.9905 / 0.9878 / 0.9485`
* `enhancedgate_kir2dl5guard_geneaware`: `0.9940 / 0.9914 / 0.9509`
* `enhancedgate_kir2dl5_kir2ds5unsupported_geneaware`: `0.9988 / 0.9961 / 0.9568`

Remaining functional-error triage is documented in
`docs/research/2026-05-11_remaining_functional_error_triage.md` and generated
by `benchmarks/scripts/inspect_remaining_functional_errors.py`. For the current
candidate, there are no remaining `candidate_regression` rows. Remaining
3/5-digit functional rows are now limited to KIR2DL1 and KIR2DS3:

* KIR2DL1: one shared 5-digit suballele miss matching discard
* KIR2DS3: one shared 3-digit/5-digit miss plus one 5-digit all-methods-shifted row

KIR2DS5 and KIR2DL5 have no remaining 3-digit or 5-digit functional rows in the
current stress sweep. Do not broaden the KIR2DS3/KIR2DS5 ambiguity handling
without a new sample-level audit showing a regression-safe target.

KIR2DL5-specific remaining-error audit is documented in
`docs/research/2026-05-11_kir2dl5_failure_audit.md` and generated by
`benchmarks/scripts/inspect_kir2dl5_remaining_failures.py`. It now shows no
remaining KIR2DL5A/B 3-digit or 5-digit functional errors for the current
candidate.

The follow-up KIR2DL5 discriminating-variant audit is documented in
`docs/research/2026-05-11_kir2dl5_discriminating_variants.md` and generated by
`benchmarks/scripts/inspect_kir2dl5_discriminating_variants.py`. It now has no
remaining rows because the current candidate has no remaining KIR2DL5 functional
errors. KIR2DL1 3-digit is fixed; the remaining KIR2DL1 issue is one shared
5-digit suballele miss.

KIR2DS5-specific remaining-error audit is documented in
`docs/research/2026-05-11_kir2ds5_failure_audit.md` and generated by
`benchmarks/scripts/inspect_kir2ds5_remaining_failures.py`. For the current
candidate, it now has no remaining KIR2DS5 3-digit or 5-digit rows. The prior
audit found truth-favoring margins in four KIR2DS5 samples, with recurrent
unsupported selected `KIR2DS5*0270102` evidence at exon variant `hv1797`; during
tuning, a residual `KIR2DS5*010` alternative was also strongly unsupported.
Keep the targeted guard constrained to selected `KIR2DS5*027/KIR2DS5*010` and
the 5-digit non-target preservation rule unless broader panels show a specific
candidate regression.

KIR2DS3-specific remaining-error audit is documented in
`docs/research/2026-05-11_kir2ds3_failure_audit.md` and generated by
`benchmarks/scripts/inspect_kir2ds3_remaining_failures.py`. It finds two
remaining KIR2DS3 samples for the current candidate:

* `synthetic-difficult5x12` sample `03`: `shared_with_discard`,
  truth `KIR2DS3*011_KIR2DS3*016`, candidate/discard/likelihood
  `KIR2DS3*011_KIR2DS3*0020101`. Variant evidence favors truth strongly:
  candidate-only `KIR2DS3*0020101` has five unsupported variants and a
  truth-candidate net margin of `964.0`. However, the truth genotype is far
  down the likelihood ranking (`rank 692`, gap `383.156`), so a near-window
  unsupported-overcall guard is unlikely to recover it by itself.
* `synthetic-difficult5x12-seed5102` sample `00`: `all_methods_disagree_or_shifted`,
  truth `KIR2DS3*0010306_KIR2DS3*0011302`, candidate
  `KIR2DS3*0010306_KIR2DS3*0010301`, discard
  `KIR2DS3*00108_KIR2DS3*0011302`, likelihood
  `KIR2DS3*0010306_KIR2DS3*0020101`. Variant evidence favors the candidate
  rather than truth because the truth-only `KIR2DS3*0011302` exon variant is
  unsupported.

Do not add a broad KIR2DS3 rollback or a generic KIR2DS3*002 guard from these
two rows alone. The next KIR2DS3 method work should first inspect why
`KIR2DS3*016` is ranked so far below the selected `KIR2DS3*0020101` in sample
`03`, while leaving the seed5102 sample `00` case unchanged unless new evidence
shows a safe way to rescue `00113`.

Operational note: `benchmarks/configs/hprc_real_sanity.json` currently uses
`examples/cohort.csv`, so it is an examples-format smoke run rather than a valid
HPRC accuracy benchmark. Use it to verify the full legacy/rerun plumbing only.
For enhancedgate smoke on that full gene panel, use the committed
`benchmarks/configs/hprc_real_sanity_enhancedgate.json`; applying top5000 to
every gene can exceed a 15GB WSL memory limit. The committed enhancedgate sanity
config includes gene-aware top-n, the `KIR2DL1` functional fallback, and the
`KIR2DS5/KIR2DS3` promotion guards. Synthetic profiling showed plain
`base_top_n = 600` can regress
`synthetic-functional8x6` KIR2DL1, while the gene-aware setting
`base_top_n = 600` plus `KIR2DL1:1000` recovered the observed aggregate
3/5/7-digit F1 and improved runtime.

To prepare a real HPRC mini sanity run, use
`benchmarks/scripts/prepare_hprc_real_mini.py`. It checks
`data/cohorts/hprc.csv` against `data/groundtruth/hprc_summary_v1_2_e.tsv`,
discovers paired FASTQs under `GRAPHKIR_HPRC_FASTQ_ROOTS`,
`HPRC_FASTQ_ROOT`, or the default local roots, and writes ignored generated
manifest/config files only when real sample data exists. It writes a baseline
config and an enhancedgate config with `allele_base_top_n = 600`,
`allele_gene_base_top_ns = KIR2DL1:1000`, the `KIR2DL1` functional fallback, and
the `KIR2DS5/KIR2DS3` promotion guards.
As of
2026-05-11, this workspace does not contain local HPRC KIR FASTQs or `data_real`
intermediates.

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
