# GraphKir2 Roadmap

## Phase 0 — Refactor Scaffolding and Benchmark Surface

* keep legacy `graphkir/` and `kir/` stable
* scaffold `src/graphkir2/` modules
* formalize manifests, config objects, and planning outputs
* create canonical benchmark compare artifact schema
* establish `docs/architecture.md`, `docs/benchmark_plan.md`, and ADR/research layout

Status:

* in progress

## Phase 1 — Mapping Evidence Refactor

* turn mapping into an explicit execution stage, not just a naming plan
* separate graph alignment, filtering, and depth artifacts
* add benchmark hooks for mapping-stage resource measurements
* implement ablation-ready multi-map policy surface:
  * `discard`
  * `best-only`
  * `weighted`

Primary success criteria:

* old-vs-new mapping parity on smoke and simulated inputs
* no runtime collapse relative to legacy Graph-KIR

## Phase 2 — Copy Number and Typing Contracts

* make CN execution independent from legacy filename coupling
* expose per-sample vs cohort CN cleanly
* formalize typing inputs, outputs, and confidence states
* add benchmark summaries for `3-digit/5-digit/7-digit` separately
* focus method tuning on functional typing instead of only defending `7-digit`

Primary success criteria:

* measurable improvement toward `3-digit/5-digit > Geny`
* preserve CN strength
* keep typing output merge behavior stable

## Phase 3 — Structural and Novel Second Pass

* promote structural/novel handling into a deliberate stage
* detect long deletions, fusion-like patterns, and novel-like mismatches
* attach warnings and downgrade behavior instead of forcing overconfident calls

Primary success criteria:

* improved difficult-case recall
* clearer failure-mode reports on real data

## Phase 4 — Benchmark-Driven Method Tuning

Synthetic-first execution order:

1. synthetic profiling and parameter sweeps
2. synthetic difficult-case ablations
3. real-data sanity on HPRC-like cohorts
4. broader real-data confirmation only after a candidate method is stable

* run focused ablations from `benchmarks/configs/`
* collect runtime and memory metrics
* compare baseline vs refactored outputs on simulated and HPRC data
* treat HPRC acquisition as a separate provenance gate:
  * verify `fastq_md5` against at least two independent official download paths when possible
  * classify full-size but mismatched files as `source-md5 inconsistency`, not generic download failure
  * only promote real-data runs after the FASTQ provenance is settled
* produce per-gene and failure-mode summaries

Primary success criteria:

* measurable gain in either:
  * `3-digit/5-digit` functional typing
  * `7-digit` performance
  * difficult-case robustness
  * runtime/resource efficiency

## Phase 5 — Promotion Decision

* decide whether `graphkir2` replaces `graphkir`
* freeze a release benchmark suite
* document migration path for users and paper workflow scripts

Promotion rule:

* do not replace the legacy CLI until benchmark evidence shows acceptable parity or
  meaningful improvement on the chosen target metrics
