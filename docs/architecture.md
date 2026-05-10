# GraphKir2 Architecture

## Goal

`graphkir2` keeps the current `graphkir/` implementation as a stable baseline while
incrementally rebuilding the pipeline around clearer stage boundaries, explicit
artifacts, and benchmarkable failure modes.

The main optimization targets are:

* more accurate handling of multi-mapped reads instead of hard discard
* better recovery on structural/novel KIR cases that currently depress real-data
  `7-digit` performance
* preservation of the current Graph-KIR speed advantage during refactor work

## Pipeline

```text
Input manifest
  ->
input / extraction planning
  ->
reference / index planning
  ->
mapping
  ->
variant evidence + depth
  ->
copy number
  ->
typing
  ->
novel / structural refinement
  ->
merged cohort outputs
  ->
benchmark summary
```

## Source Modules

| Module | Responsibility |
|---|---|
| `config/` | runtime, index, typing, CN, and input configuration objects |
| `io/` | manifest parsing, output layout rules, stable path naming |
| `mapping/` | graph mapping contracts, multi-map policies, mapping artifacts |
| `cn/` | depth-derived copy-number estimation and cohort merge behavior |
| `typing/` | allele typing plans and confidence-aware output contracts |
| `benchmark/` | preset loading, legacy command parity, result artifact schema |
| `core/` | pipeline composition and stage orchestration |

Expected future modules:

| Module | Planned responsibility |
|---|---|
| `extract/` | WGS-to-KIR extraction and direct-FASTQ input normalization |
| `variants/` | per-read variant evidence, depth, and filtering policies |
| `novel/` | second-pass structural and novel-allele refinement |
| `cli/` | optional separation of command handling from library surface |

## Stage Contracts

### 1. Manifest and Inputs

Input should be normalized into typed sample records before pipeline logic begins.

Current artifact:

* `SampleManifest`
* `SampleInput`

### 2. Mapping

Mapping must define stable artifact names independently of the legacy implementation.
It is the first stage where a major method improvement is expected.

Current artifacts:

* graph BAM path
* variant prefix
* depth TSV path

Planned method work:

* `discard` vs `weighted` multi-map handling
* paired-end consistency rescue
* mapping-quality-aware read evidence scoring

### 3. Copy Number

CN remains a first-class output because it is already a practical strength of the
legacy method and a major part of the Graph-KIR paper framing.

Current artifacts:

* per-sample CN TSV
* optional CN model JSON
* merged cohort CN TSV

Planned method work:

* cleaner cohort-mode normalization interface
* explicit diploid-prior policies
* per-gene confidence reporting

### 4. Typing

Typing should consume explicit variant and CN artifacts instead of implicit filename
conventions hidden across the codebase.

Current artifacts:

* allele TSV
* possible-typing TSV
* merged cohort allele TSV

Planned method work:

* confidence-aware `7-digit -> 5-digit -> warning` downgrade path
* more explicit likelihood inputs
* ablation-friendly strategy switches

### 5. Novel and Structural Refinement

Real-data performance is currently limited more by structural variation and novel
polymorphisms than by simple coverage. This stage should become a deliberate second
pass instead of a side utility.

Planned artifacts:

* structural warning report
* novel-like evidence JSON
* refined allele calls

## Artifact Layout

Benchmark output should be reproducible from manifests and written outside `docs/`.

```text
benchmarks/
  configs/
  results/
  plots/
```

Current canonical compare artifact schema:

* `metadata`
* `comparison`
* `metrics`
* `collector`

## Design Principles

* Keep `graphkir/` and `kir/` runnable as the baseline until parity is demonstrated.
* Separate planning, execution, and evaluation surfaces.
* Make every major improvement benchmarkable through manifests and result artifacts.
* Prefer typed contracts over implicit path conventions.
