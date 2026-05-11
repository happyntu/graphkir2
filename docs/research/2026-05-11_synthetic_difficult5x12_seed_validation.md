# 2026-05-11 Synthetic Difficult5x12 Seed Validation

## Goal

Check whether the conditional `KIR2DS3` ratio gate and second-gate fallback
generalize beyond the original `synthetic-difficult5x12` fixture.

## Setup

Generated three additional 12-sample difficult panels:

| panel | allele seed | read seed | samples | depth | genes |
|---|---:|---:|---:|---:|---|
| `synthetic-difficult5x12-seed5101` | 5101 | 25101 | 12 | 20x | `KIR2DL4,KIR2DL5,KIR2DS3,KIR2DS4,KIR2DS5` |
| `synthetic-difficult5x12-seed5102` | 5102 | 25102 | 12 | 20x | `KIR2DL4,KIR2DL5,KIR2DS3,KIR2DS4,KIR2DS5` |
| `synthetic-difficult5x12-seed5103` | 5103 | 25103 | 12 | 20x | `KIR2DL4,KIR2DL5,KIR2DS3,KIR2DS4,KIR2DS5` |

All panels used `msa_type = ab` with exon-only alleles disabled.

Compared:

* discard baseline
* conditional `KIR2DS3*00201` ratio gate at `0.80`
* formalized second-gate fallback

## Overall Results

| panel | method | 3-digit | 5-digit | 7-digit |
|---|---|---:|---:|---:|
| seed5101 | discard | 0.9667 | 0.9583 | 0.9250 |
| seed5101 | ratio080 | 0.9667 | 0.9667 | 0.9583 |
| seed5101 | secondgate | 0.9667 | 0.9667 | 0.9583 |
| seed5102 | discard | 0.9333 | 0.9250 | 0.8667 |
| seed5102 | ratio080 | 0.9583 | 0.9500 | 0.9333 |
| seed5102 | secondgate | 0.9583 | 0.9500 | 0.9333 |
| seed5103 | discard | 0.9750 | 0.9750 | 0.9083 |
| seed5103 | ratio080 | 0.9750 | 0.9667 | 0.9167 |
| seed5103 | secondgate | 0.9667 | 0.9583 | 0.9167 |

Mean across the three new seeds:

| method | 3-digit | 5-digit | 7-digit |
|---|---:|---:|---:|
| discard | 0.9583 | 0.9528 | 0.9000 |
| ratio080 | 0.9667 | 0.9611 | 0.9361 |
| secondgate | 0.9639 | 0.9583 | 0.9361 |

## KIR2DS3-Specific Results

| panel | method | KIR2DS3 3-digit | KIR2DS3 5-digit | KIR2DS3 7-digit |
|---|---|---:|---:|---:|
| seed5101 | discard | 0.9583 | 0.9167 | 0.9167 |
| seed5101 | ratio080 | 0.9583 | 0.9583 | 0.9583 |
| seed5101 | secondgate | 0.9583 | 0.9583 | 0.9583 |
| seed5102 | discard | 0.9583 | 0.9167 | 0.9167 |
| seed5102 | ratio080 | 0.9167 | 0.8750 | 0.8333 |
| seed5102 | secondgate | 0.9167 | 0.8750 | 0.8333 |
| seed5103 | discard | 0.9583 | 0.9583 | 0.9583 |
| seed5103 | ratio080 | 0.9167 | 0.8750 | 0.8333 |
| seed5103 | secondgate | 0.8750 | 0.8333 | 0.8333 |

Mean across the three new seeds:

| method | KIR2DS3 3-digit | KIR2DS3 5-digit | KIR2DS3 7-digit |
|---|---:|---:|---:|
| discard | 0.9583 | 0.9306 | 0.9306 |
| ratio080 | 0.9306 | 0.9028 | 0.8750 |
| secondgate | 0.9167 | 0.8889 | 0.8750 |

## Interpretation

The ratio gate and secondgate improve mean whole-panel `7-digit` accuracy, but
they do not satisfy the functional typing goal because `KIR2DS3` 3/5-digit
performance is worse than discard on the broader seeds.

`seed5103` is the clearest warning case:

* discard has the best `3-digit` and `5-digit` overall scores
* ratio080 improves only `7-digit` but regresses `5-digit`
* secondgate regresses both `3-digit` and `5-digit` compared with ratio080 and
  discard

## Decision

Do not promote secondgate to a broader candidate and do not move to real-data
sanity validation yet.

Next method work should target `KIR2DS3` robustness across seeds. The immediate
goal is to preserve the whole-panel `7-digit` gain without sacrificing
`KIR2DS3` 3/5-digit functional calls versus discard.

## Enhanced Fallback Follow-Up

A refined `enhancedgate` was tested after the secondgate failure. It keeps the
conditional ratio rescue but tightens discard-style fallback:

* fallback requires base private-support score `<= -20.0`
* residual `KIR2DS3*00201` fallback requires cross-support ratio `>= 0.70`,
  unless the base call already contains `KIR2DS3*00103`
* introduced `KIR2DS3*00103` fallback uses cross-support ratio `< 0.885`

Primary panel checks:

| panel | 3-digit | 5-digit | 7-digit | KIR2DS3 3-digit | KIR2DS3 5-digit | KIR2DS3 7-digit |
|---|---:|---:|---:|---:|---:|---:|
| synthetic-difficult5 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | 0.9750 | 0.9750 | 0.9500 | 0.9583 | 0.9583 | 0.9167 |
| synthetic-functional8 | 1.0000 | 1.0000 | 0.9063 | n/a | n/a | n/a |
| synthetic-functional8x6 | 0.9896 | 0.9896 | 0.9063 | n/a | n/a | n/a |

Across the original `synthetic-difficult5x12` plus the three new seeds:

| method | 3-digit | 5-digit | 7-digit | KIR2DS3 3-digit | KIR2DS3 5-digit | KIR2DS3 7-digit |
|---|---:|---:|---:|---:|---:|---:|
| discard | 0.9604 | 0.9542 | 0.9062 | 0.9583 | 0.9271 | 0.9271 |
| ratio080 | 0.9646 | 0.9604 | 0.9354 | 0.9167 | 0.8958 | 0.8646 |
| secondgate | 0.9667 | 0.9625 | 0.9396 | 0.9271 | 0.9063 | 0.8854 |
| enhancedgate | 0.9792 | 0.9750 | 0.9500 | 0.9896 | 0.9688 | 0.9375 |

Decision: `enhancedgate` is the current synthetic lead. It is strong enough for
a small real-data sanity check, but not broad benchmarking or default promotion.
