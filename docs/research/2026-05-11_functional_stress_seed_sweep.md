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
| synthetic-functional8 | discard | 3.38 | 1.0000 | 1.0000 | 0.8438 |
| synthetic-functional8 | likelihood_top5000 | 4.95 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_top5000 | 8.82 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_geneaware | 4.84 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1fallback_geneaware | 5.02 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | discard | 7.12 | 0.9688 | 0.9583 | 0.8542 |
| synthetic-functional8x6 | likelihood_top5000 | 16.39 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_top5000 | 16.70 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_geneaware | 16.69 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_kir2dl1fallback_geneaware | 14.31 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-difficult5 | discard | 3.53 | 0.9750 | 0.9750 | 0.8750 |
| synthetic-difficult5 | likelihood_top5000 | 7.86 | 0.9750 | 0.9750 | 0.9250 |
| synthetic-difficult5 | enhancedgate_top5000 | 7.02 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_geneaware | 6.44 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1fallback_geneaware | 6.19 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | discard | 11.40 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12 | likelihood_top5000 | 13.84 | 0.9250 | 0.9250 | 0.9083 |
| synthetic-difficult5x12 | enhancedgate_top5000 | 26.08 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_geneaware | 21.46 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1fallback_geneaware | 18.30 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12-seed5101 | discard | 11.83 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12-seed5101 | likelihood_top5000 | 12.53 | 0.9000 | 0.9000 | 0.8833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_top5000 | 26.03 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_geneaware | 24.33 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1fallback_geneaware | 20.22 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5102 | discard | 11.85 | 0.9333 | 0.9250 | 0.8667 |
| synthetic-difficult5x12-seed5102 | likelihood_top5000 | 12.81 | 0.9083 | 0.9000 | 0.8750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_top5000 | 25.04 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_geneaware | 18.54 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1fallback_geneaware | 21.63 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5103 | discard | 7.72 | 0.9750 | 0.9750 | 0.9083 |
| synthetic-difficult5x12-seed5103 | likelihood_top5000 | 14.93 | 0.9083 | 0.9083 | 0.8583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_top5000 | 24.84 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_geneaware | 18.18 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1fallback_geneaware | 21.49 | 0.9917 | 0.9833 | 0.9417 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| discard | 7 | 8.12 | 0.9693 | 0.9643 | 0.8854 |
| enhancedgate_geneaware | 7 | 15.78 | 0.9866 | 0.9842 | 0.9446 |
| enhancedgate_kir2dl1fallback_geneaware | 7 | 15.31 | 0.9881 | 0.9842 | 0.9461 |
| enhancedgate_top5000 | 7 | 19.22 | 0.9866 | 0.9842 | 0.9446 |
| likelihood_top5000 | 7 | 11.90 | 0.9437 | 0.9426 | 0.8946 |

## KIR2DL1 Fallback Gene-aware Vs Discard

No per-gene 3/5-digit regressions versus discard were detected.

Per-gene 3/5-digit gains versus discard:

| panel | gene | metric | discard F1 | fallback gene-aware F1 | delta |
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

Remaining fallback gene-aware 3/5-digit errors:

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

## Decision

`enhancedgate_kir2dl1fallback_geneaware` resolves the strict
`synthetic-functional8x6` KIR2DL1 3-digit regression while preserving
the enhancedgate KIR2DS3/KIR2DS4/KIR2DS5 gains and the gene-aware
top-n runtime setting.
`likelihood_top5000` alone is not viable because it loses substantial
3/5-digit accuracy on the difficult5x12 seed panels.
KIR2DL1 still has one 5-digit miss matching discard's remaining error;
the next method work should focus on the remaining 5-digit errors in
`KIR2DL5A/B`, `KIR2DS3`, and `KIR2DS5` without reintroducing the
KIR2DL1 3-digit ambiguity-likelihood regression.
