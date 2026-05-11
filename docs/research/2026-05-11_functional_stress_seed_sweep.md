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
| synthetic-functional8 | discard | 4.65 | 1.0000 | 1.0000 | 0.8438 |
| synthetic-functional8 | likelihood_top5000 | 13.81 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_top5000 | 9.25 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_geneaware | 7.47 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1fallback_geneaware | 10.29 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 9.04 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_functionalguard_geneaware | 8.80 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | discard | 13.75 | 0.9688 | 0.9583 | 0.8542 |
| synthetic-functional8x6 | likelihood_top5000 | 17.53 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_top5000 | 25.34 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_geneaware | 21.79 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_kir2dl1fallback_geneaware | 27.75 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 24.96 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_functionalguard_geneaware | 24.94 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-difficult5 | discard | 8.15 | 0.9750 | 0.9750 | 0.8750 |
| synthetic-difficult5 | likelihood_top5000 | 6.06 | 0.9750 | 0.9750 | 0.9250 |
| synthetic-difficult5 | enhancedgate_top5000 | 8.78 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_geneaware | 6.84 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1fallback_geneaware | 10.50 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 8.32 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_functionalguard_geneaware | 9.71 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | discard | 11.56 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12 | likelihood_top5000 | 13.77 | 0.9250 | 0.9250 | 0.9083 |
| synthetic-difficult5x12 | enhancedgate_top5000 | 25.40 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_geneaware | 17.90 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1fallback_geneaware | 17.80 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 25.63 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_functionalguard_geneaware | 29.84 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12-seed5101 | discard | 7.20 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12-seed5101 | likelihood_top5000 | 11.29 | 0.9000 | 0.9000 | 0.8833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_top5000 | 23.48 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_geneaware | 20.37 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1fallback_geneaware | 17.26 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 25.45 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5101 | enhancedgate_functionalguard_geneaware | 31.84 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5102 | discard | 7.38 | 0.9333 | 0.9250 | 0.8667 |
| synthetic-difficult5x12-seed5102 | likelihood_top5000 | 15.58 | 0.9083 | 0.9000 | 0.8750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_top5000 | 27.14 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_geneaware | 19.88 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1fallback_geneaware | 20.83 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 25.61 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_functionalguard_geneaware | 30.26 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5103 | discard | 6.84 | 0.9750 | 0.9750 | 0.9083 |
| synthetic-difficult5x12-seed5103 | likelihood_top5000 | 12.20 | 0.9083 | 0.9083 | 0.8583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_top5000 | 26.68 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_geneaware | 22.84 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1fallback_geneaware | 22.02 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 24.83 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_functionalguard_geneaware | 33.06 | 0.9917 | 0.9917 | 0.9500 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| discard | 7 | 8.50 | 0.9693 | 0.9643 | 0.8854 |
| enhancedgate_functionalguard_geneaware | 7 | 24.06 | 0.9905 | 0.9878 | 0.9485 |
| enhancedgate_geneaware | 7 | 16.73 | 0.9866 | 0.9842 | 0.9446 |
| enhancedgate_kir2dl1_kir2ds5guard_geneaware | 7 | 20.55 | 0.9905 | 0.9866 | 0.9473 |
| enhancedgate_kir2dl1fallback_geneaware | 7 | 18.06 | 0.9881 | 0.9842 | 0.9461 |
| enhancedgate_top5000 | 7 | 20.87 | 0.9866 | 0.9842 | 0.9446 |
| likelihood_top5000 | 7 | 12.89 | 0.9437 | 0.9426 | 0.8946 |

## Functional Guard Gene-aware Vs Discard

No per-gene 3/5-digit regressions versus discard were detected.

Per-gene 3/5-digit gains versus discard:

| panel | gene | metric | discard F1 | candidate F1 | delta |
|---|---:|---:|---:|---:|---:|
| synthetic-difficult5 | KIR2DS3 | three_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5 | KIR2DS3 | five_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5x12 | KIR2DS3 | five_digit | 0.9167 | 0.9583 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | three_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | three_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | five_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9167 | 0.9583 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | three_digit | 0.7500 | 0.8750 | 0.1250 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | five_digit | 0.7500 | 0.8750 | 0.1250 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-functional8x6 | KIR2DS1 | three_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-functional8x6 | KIR2DS1 | five_digit | 0.7500 | 1.0000 | 0.2500 |

Remaining functional-guard 3/5-digit errors:

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
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9583 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | five_digit | 0.8750 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | three_digit | 0.8750 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | five_digit | 0.9167 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | three_digit | 0.9167 |
| synthetic-functional8x6 | KIR2DL1 | five_digit | 0.9167 |

## Decision

`enhancedgate_functionalguard_geneaware` resolves the strict
`synthetic-functional8x6` KIR2DL1 3-digit regression, keeps the
KIR2DS5 promotion guard, adds a narrow KIR2DS3 suballele guard,
and keeps the gene-aware top-n runtime setting.
`likelihood_top5000` alone is not viable because it loses substantial
3/5-digit accuracy on the difficult5x12 seed panels.
KIR2DL1 still has one 5-digit miss matching discard's remaining error;
the next method work should focus on the remaining 5-digit errors in
`KIR2DL5A/B`, `KIR2DS3`, and `KIR2DS5` without reintroducing the
KIR2DL1 3-digit ambiguity-likelihood regression.
