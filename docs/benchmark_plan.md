# GraphKir2 Benchmark Plan

## Benchmark Goals

The benchmark program for `graphkir2` must answer four separate questions:

1. Does the refactor preserve legacy behavior where parity is expected?
2. Does a method change improve `3-digit/5-digit` functional typing against `Geny`?
3. Which failure modes improve: multi-mapped reads, structural events, or novel cases?
4. Does the change keep or improve runtime and memory behavior?

Primary development target:

* `graphkir2` should beat `Geny` on `3-digit` and `5-digit` functional typing

Secondary guardrails:

* do not regress on copy number
* do not regress on `7-digit`
* preserve practical runtime advantages where possible

## Benchmark Layers

| Layer | Purpose |
|---|---|
| smoke | verify manifests, naming, and CLI behavior |
| planning | verify old-vs-new stage parity before execution code exists |
| ablation | isolate one method change at a time |
| real-data sanity | confirm no obvious regression on HPRC-like inputs |
| release comparison | summarize baseline vs refactored outputs and resources |

## Dataset Progression

| Step | Dataset | Purpose |
|---|---|---|
| 1 | `examples/cohort.csv` | smoke and planning validation |
| 2 | simulated short-read cohort | fast profiling, CN/typing regression checks, and parameter tuning |
| 3 | simulated difficult subset | multi-map / structural / novel ablation |
| 4 | HPRC short-read cohort | real-data sanity after simulated tuning stabilizes |
| 5 | expanded real-data cohort | release-level comparison |

Large benchmark inputs should remain under `data/` or external manifests, not committed
directly into `benchmarks/results/`.

## Synthetic-First Strategy

`graphkir2` should follow a synthetic-first tuning loop similar in spirit to the current
`hybridseed` workflow:

1. use synthetic data to find runtime bottlenecks and unstable method behavior
2. use synthetic data to sweep parameters and compare ablations cheaply
3. freeze a small number of promising settings
4. only then run HPRC or other real-data validation

Why this order is preferred:

* synthetic data is cheaper for repeated profiling
* truth is explicit, so `3-digit/5-digit/7-digit` changes are easier to interpret
* runtime regressions can be caught before spending time on large real cohorts
* real data is then reserved for sanity and generalization checks, not blind search

Real-data benchmarks should not be the first place where method tuning happens.

## Benchmark Matrix

### Current preset categories

* `smoke_examples`
* `simulated_cn_ablation`
* `simulated_typing_ablation`
* `simulated_runtime_tuning`
* `simulated_difficult_cases`
* `hprc_real_sanity`
* `multimap_rescue_ablation`
* `structural_novel_ablation`
* `thread_scaling`
* `index_reuse`

### Method comparisons to support

* legacy `graphkir` vs `graphkir2`
* `graphkir2` vs `Geny` on `3-digit` and `5-digit`
* `discard multi-map` vs `weighted multi-map`
* `full` vs `exonfirst`
* per-sample CN vs cohort CN
* extraction on vs extraction skipped

## Metrics

### Functional metrics

* `3_digit_f1`
* `5_digit_f1`
* `copy_number_f1`
* `allele_f1`
* `typing_call_rate`
* `7_digit_call_rate`
* `downgraded_call_rate`

### Failure-mode metrics

* `multimap_affected_samples`
* `structural_warning_rate`
* `novel_like_rate`
* `fusion_case_recall`
* `segmental_deletion_recall`

### Resource metrics

* `runtime_seconds`
* `peak_rss_mb`
* `index_reuse_hit_rate`
* `samples_per_hour`

### Synthetic-first profiling metrics

* stage-level runtime split
* per-sample runtime variance
* mapping artifact size
* allele candidate-set size
* depth-processing overhead

## Result Layout

Per benchmark label, use:

```text
benchmarks/results/<label>/
  compare.json
  summary.tsv
  per_gene.tsv
  failure_modes.tsv
  notes.md
```

During the current planning-only phase, a single JSON artifact is acceptable:

* `<label>.compare.json`

## Interpretation Rules

* Do not report a `graphkir2` planning artifact as an executed benchmark.
* Do not collapse CN and allele typing into a single headline score.
* Keep `3-digit/5-digit` and `7-digit` claims separate.
* For refactor decisions, prefer changes that improve `3-digit/5-digit` even if they do
  not improve `7-digit`, as long as CN and runtime remain acceptable.
* Track speed claims separately from accuracy claims.
* Treat multi-map rescue and structural refinement as distinct interventions.
* Prefer synthetic-data findings for early tuning decisions, then confirm them on real
  data before treating them as durable improvements.
