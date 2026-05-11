# 2026-05-11 Synthetic Difficult5x12 Overfit Check

## Goal

Check whether the functional-target rules tuned on the 4-sample
`synthetic-difficult5` panel generalize to a larger synthetic difficult panel.

## Setup

Generated `synthetic-difficult5x12` with:

* genes: `KIR2DL4,KIR2DL5,KIR2DS3,KIR2DS4,KIR2DS5`
* samples: `12`
* depth: `20x`
* `msa_type = ab`
* exon-only alleles disabled
* allele seed: `4242`
* read seed: `22000`

Config files:

* `benchmarks/configs/synthetic-difficult5x12.json`
* `benchmarks/configs/synthetic-difficult5x12-likelihood-exon2-minfrac07.json`
* `benchmarks/configs/synthetic-difficult5x12-functional-target-kir2ds3-private.json`

## Cohort Results

| method | 3-digit | 5-digit | 7-digit |
|---|---:|---:|---:|
| discard baseline | 0.9667 | 0.9583 | 0.9250 |
| balanced likelihood top5000 | 0.9250 | 0.9250 | 0.9083 |
| functional target | 0.9417 | 0.9417 | 0.9167 |

## Per-Gene Findings

The main regression is `KIR2DS3`.

| method | KIR2DS3 3-digit | KIR2DS3 5-digit | KIR2DS3 7-digit |
|---|---:|---:|---:|
| discard baseline | 0.9583 | 0.9167 | 0.9167 |
| balanced likelihood top5000 | 0.7083 | 0.7083 | 0.7083 |
| functional target | 0.7917 | 0.7917 | 0.7500 |

Other genes under functional target:

* `KIR2DL4 = 1.0 / 1.0 / 1.0`
* `KIR2DL5A = 0.9744 / 0.9744 / 0.9744`
* `KIR2DL5B = 0.8889 / 0.8889 / 0.8889`
* `KIR2DS4 = 1.0 / 1.0 / 0.9167`
* `KIR2DS5 = 0.9583 / 0.9583 / 0.9583`

## KIR2DS3 Error Pattern

Balanced likelihood frequently pulls one `KIR2DS3` allele toward
`KIR2DS3*00201`.

Examples:

* sample `00`: truth `00103 + 00103`, balanced predicts `00103 + 00201`
* sample `01`: truth `00103 + 00103`, balanced predicts `00103 + 00201`
* sample `02`: truth `00103 + 00103`, balanced predicts `00103 + 00201`
* sample `07`: truth `017 + 015`, balanced predicts `017 + 00201`
* sample `09`: truth `00103 + 016`, balanced predicts `00103 + 00201`

Functional target fixes some but not all of these and introduces additional
misses:

* sample `05`: truth `010 + 00201`, functional predicts `010 + 010`
* sample `06`: truth `00201 + 00201`, functional predicts `00201 + 00109`
* sample `11`: truth `012 + 019`, functional predicts `012 + 00103`

## Interpretation

The 4-sample `synthetic-difficult5` success was real for that fixture, but it is
not enough evidence to promote the current functional target as a global default.

The current functional-target rules are too broad:

* likelihood handling improves some cross-gene ambiguity cases
* `KIR2DS3` private-support reranking rescues the original sample `02`
* but on larger synthetic coverage, the same surface still over-selects or
  over-corrects `KIR2DS3` near-neighbor alleles

## Decision

Keep the current functional-target config as an experimental targeted candidate,
not the default.

Next method work should make `KIR2DS3` private-support and
`KIR2DS3/KIR2DS5` directional neutralization conditional on explicit
contamination evidence. A global per-gene rule is not robust enough.

## Conditional Rescue Ablation

Tested and then formalized a first contamination-gated rescue rule:

* start from balanced likelihood selection
* only enable `KIR2DS3` private-support rescue when the selected allele set has
  a `KIR2DS3*00201` allele
* require that the selected allele set's private-variant positive support is
  mostly carried by `KIR2DS3/KIR2DS5` cross-gene ambiguous reads
* sweep cross-support ratio thresholds: `0.80`, `0.85`, `0.90`

The `0.80` gate is now reproducible through:

* `benchmarks/configs/synthetic-difficult5-conditional-kir2ds3-ratio080.json`
* `benchmarks/configs/synthetic-difficult5x12-conditional-kir2ds3-ratio080.json`

Implementation detail: the gate is evaluated on sample-level reads so the
`KIR2DS3/KIR2DS5` physical-read ambiguity is visible, then the target gene is
rerun with directional neutralization only if the gate passes.

Results:

| method | 3-digit | 5-digit | 7-digit | KIR2DS3 5-digit |
|---|---:|---:|---:|---:|
| discard baseline | 0.9667 | 0.9583 | 0.9250 | 0.9167 |
| balanced likelihood top5000 | 0.9250 | 0.9250 | 0.9083 | 0.7083 |
| functional target | 0.9417 | 0.9417 | 0.9167 | 0.7917 |
| conditional rescue, ratio >= 0.80 | 0.9583 | 0.9583 | 0.9333 | 0.8750 |
| conditional rescue, ratio >= 0.85 | 0.9500 | 0.9500 | 0.9250 | 0.8333 |
| conditional rescue, ratio >= 0.90 | 0.9250 | 0.9250 | 0.9083 | 0.7083 |

Interpretation:

* the ratio-gated rescue is better than the broad functional-target rule
* threshold `0.80` matches discard at `5-digit` and exceeds discard at `7-digit`
* it still trails discard at `3-digit` because `KIR2DS3` remains weaker
* a single cross-support ratio is not sufficient to decide all `KIR2DS3`
  rescue cases

Decision:

Do not promote the ratio gate as a final default yet. Keep it as the current
targeted experimental candidate and add a second criterion that avoids changing
cases where discard already makes the correct `KIR2DS3` functional call.

## Second-Gate Fallback Ablation

The remaining `synthetic-difficult5x12` `KIR2DS3` errors after ratio `0.80`
are:

* sample `01`: discard correct, conditional still keeps `KIR2DS3*00201`
* sample `03`: discard also wrong, so this is not solved by fallback
* sample `11`: discard correct, conditional rescue introduces
  `KIR2DS3*00103`

Tested and then formalized a second gate that falls back to discard-style
`KIR2DS3` evidence only when the first ratio gate already fired and one of
these residual-risk patterns is present:

* conditional call still contains `KIR2DS3*00201`
* conditional rescue introduces `KIR2DS3*00103` when the base likelihood call
  did not have `00103`, and the cross-support ratio is `< 0.90`

The formalized implementation recomputes the target gene inside
`rerun_typing_private_support.py` with `removeMultipleMapped` evidence rather
than reading an external discard prediction table.

Results:

| panel | 3-digit | 5-digit | 7-digit | fallback samples |
|---|---:|---:|---:|---|
| synthetic-difficult5 | 1.0000 | 1.0000 | 1.0000 | `01` |
| synthetic-difficult5x12 | 0.9750 | 0.9750 | 0.9500 | `01,08,11` |
| synthetic-functional8 | 1.0000 | 1.0000 | 0.9063 | none |
| synthetic-functional8x6 | 0.9896 | 0.9896 | 0.9063 | none |

Interpretation:

* this is the best current synthetic result for the difficult panels
* it beats the discard baseline on `synthetic-difficult5x12`
* it does not regress the existing functional panels tested here
* it is now reproducible through the
  `benchmarks/configs/*-conditional-kir2ds3-secondgate.json` presets
* it should remain a targeted synthetic candidate, not the global default,
  until broader seeds and real-data sanity panels are checked
