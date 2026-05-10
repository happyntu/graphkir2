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
