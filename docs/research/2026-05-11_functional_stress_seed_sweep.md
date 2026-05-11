# Functional Stress Seed Sweep

This synthetic-first sweep re-runs only the typing and evaluation stages
on existing generated fixtures. It is intended to expose 3-digit and
5-digit functional typing weaknesses before spending time on real data.

Command:

```bash
python benchmarks/scripts/run_functional_stress_sweep.py
```

Summary TSV: `benchmarks/results/functional-stress-sweep/summary.tsv`
Per-gene TSV: `benchmarks/results/functional-stress-sweep/per_gene_summary.tsv`

## Aggregate Results

| panel | method | runtime s | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| synthetic-functional8 | discard | 4.67 | 1.0000 | 1.0000 | 0.8438 |
| synthetic-functional8 | likelihood_top5000 | 5.50 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_top5000 | 6.45 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_geneaware | 5.31 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | discard | 8.18 | 0.9688 | 0.9583 | 0.8542 |
| synthetic-functional8x6 | likelihood_top5000 | 17.18 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_top5000 | 26.08 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_geneaware | 14.90 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-difficult5 | discard | 3.77 | 0.9750 | 0.9750 | 0.8750 |
| synthetic-difficult5 | likelihood_top5000 | 7.24 | 0.9750 | 0.9750 | 0.9250 |
| synthetic-difficult5 | enhancedgate_top5000 | 7.30 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_geneaware | 6.15 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | discard | 9.53 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12 | likelihood_top5000 | 13.55 | 0.9250 | 0.9250 | 0.9083 |
| synthetic-difficult5x12 | enhancedgate_top5000 | 26.06 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_geneaware | 19.57 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12-seed5101 | discard | 10.67 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12-seed5101 | likelihood_top5000 | 12.98 | 0.9000 | 0.9000 | 0.8833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_top5000 | 26.92 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_geneaware | 24.58 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5102 | discard | 8.68 | 0.9333 | 0.9250 | 0.8667 |
| synthetic-difficult5x12-seed5102 | likelihood_top5000 | 18.58 | 0.9083 | 0.9000 | 0.8750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_top5000 | 29.07 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_geneaware | 20.37 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5103 | discard | 12.42 | 0.9750 | 0.9750 | 0.9083 |
| synthetic-difficult5x12-seed5103 | likelihood_top5000 | 14.38 | 0.9083 | 0.9083 | 0.8583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_top5000 | 27.91 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_geneaware | 25.43 | 0.9917 | 0.9833 | 0.9417 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| discard | 7 | 8.27 | 0.9693 | 0.9643 | 0.8854 |
| enhancedgate_geneaware | 7 | 16.62 | 0.9866 | 0.9842 | 0.9446 |
| enhancedgate_top5000 | 7 | 21.40 | 0.9866 | 0.9842 | 0.9446 |
| likelihood_top5000 | 7 | 12.77 | 0.9437 | 0.9426 | 0.8946 |

## Enhancedgate Gene-aware Vs Discard

Per-gene 3/5-digit regressions versus discard:

| panel | gene | metric | discard F1 | enhancedgate gene-aware F1 | delta |
|---|---:|---:|---:|---:|---:|
| synthetic-functional8x6 | KIR2DL1 | three_digit | 1.0000 | 0.9167 | -0.0833 |

Per-gene 3/5-digit gains versus discard:

| panel | gene | metric | discard F1 | enhancedgate gene-aware F1 | delta |
|---|---:|---:|---:|---:|---:|
| synthetic-difficult5 | KIR2DS3 | three_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5 | KIR2DS3 | five_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5x12 | KIR2DS3 | five_digit | 0.9167 | 0.9583 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | three_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | five_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9167 | 0.9583 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | three_digit | 0.7500 | 0.8750 | 0.1250 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | five_digit | 0.7500 | 0.8750 | 0.1250 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-functional8x6 | KIR2DS1 | three_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-functional8x6 | KIR2DS1 | five_digit | 0.7500 | 1.0000 | 0.2500 |

Remaining enhancedgate gene-aware 3/5-digit errors:

| panel | gene | metric | F1 |
|---|---:|---:|---:|
| synthetic-difficult5x12 | KIR2DL5A | five_digit | 0.9744 |
| synthetic-difficult5x12 | KIR2DL5A | three_digit | 0.9744 |
| synthetic-difficult5x12 | KIR2DL5B | five_digit | 0.8889 |
| synthetic-difficult5x12 | KIR2DL5B | three_digit | 0.8889 |
| synthetic-difficult5x12 | KIR2DS3 | five_digit | 0.9583 |
| synthetic-difficult5x12 | KIR2DS3 | three_digit | 0.9583 |
| synthetic-difficult5x12 | KIR2DS5 | five_digit | 0.9583 |
| synthetic-difficult5x12 | KIR2DS5 | three_digit | 0.9583 |
| synthetic-difficult5x12-seed5101 | KIR2DL5A | five_digit | 0.9412 |
| synthetic-difficult5x12-seed5101 | KIR2DL5A | three_digit | 0.9412 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | five_digit | 0.9167 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | three_digit | 0.9167 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9583 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | five_digit | 0.8750 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | three_digit | 0.8750 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | five_digit | 0.9167 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | three_digit | 0.9167 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | five_digit | 0.9583 |
| synthetic-functional8x6 | KIR2DL1 | five_digit | 0.9167 |
| synthetic-functional8x6 | KIR2DL1 | three_digit | 0.9167 |

## Decision

`enhancedgate_geneaware` is the current aggregate synthetic lead, but it is
not a final default. It improves mean 3/5/7-digit F1 and preserves the
gene-aware top-n runtime setting, yet it still has a strict per-gene
`KIR2DL1` 3-digit regression versus discard on `synthetic-functional8x6`.
`likelihood_top5000` alone is not viable because it loses substantial
3/5-digit accuracy on the difficult5x12 seed panels.
The next method work should isolate that KIR2DL1 ambiguity-likelihood
failure while preserving the `KIR2DS3`, `KIR2DS4`, and `KIR2DS5` gains.
