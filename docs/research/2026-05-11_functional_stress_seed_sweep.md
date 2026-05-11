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
| synthetic-functional8 | discard | 4.80 | 1.0000 | 1.0000 | 0.8438 |
| synthetic-functional8 | likelihood_top5000 | 7.61 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_top5000 | 8.75 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_geneaware | 7.30 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1fallback_geneaware | 7.09 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 9.38 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_functionalguard_geneaware | 9.10 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5guard_geneaware | 9.74 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 10.42 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 11.63 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | discard | 12.72 | 0.9688 | 0.9583 | 0.8542 |
| synthetic-functional8x6 | likelihood_top5000 | 16.88 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_top5000 | 21.25 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_geneaware | 16.86 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_kir2dl1fallback_geneaware | 18.69 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 23.62 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_functionalguard_geneaware | 24.12 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5guard_geneaware | 26.20 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 27.09 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 27.44 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-difficult5 | discard | 6.36 | 0.9750 | 0.9750 | 0.8750 |
| synthetic-difficult5 | likelihood_top5000 | 8.27 | 0.9750 | 0.9750 | 0.9250 |
| synthetic-difficult5 | enhancedgate_top5000 | 9.41 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_geneaware | 8.36 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1fallback_geneaware | 7.36 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 9.12 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_functionalguard_geneaware | 11.06 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5guard_geneaware | 11.05 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 10.72 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 9.99 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | discard | 9.12 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12 | likelihood_top5000 | 13.98 | 0.9250 | 0.9250 | 0.9083 |
| synthetic-difficult5x12 | enhancedgate_top5000 | 24.14 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_geneaware | 22.42 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1fallback_geneaware | 25.42 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 33.15 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_functionalguard_geneaware | 42.73 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5guard_geneaware | 34.61 | 0.9833 | 0.9833 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 32.30 | 0.9917 | 0.9917 | 0.9583 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 34.80 | 1.0000 | 1.0000 | 0.9667 |
| synthetic-difficult5x12-seed5101 | discard | 9.65 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12-seed5101 | likelihood_top5000 | 15.33 | 0.9000 | 0.9000 | 0.8833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_top5000 | 26.19 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_geneaware | 21.87 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1fallback_geneaware | 21.49 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 27.08 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5101 | enhancedgate_functionalguard_geneaware | 32.85 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5guard_geneaware | 32.12 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 31.62 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 33.34 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5102 | discard | 8.83 | 0.9333 | 0.9250 | 0.8667 |
| synthetic-difficult5x12-seed5102 | likelihood_top5000 | 14.35 | 0.9083 | 0.9000 | 0.8750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_top5000 | 25.99 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_geneaware | 21.95 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1fallback_geneaware | 22.19 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 28.61 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_functionalguard_geneaware | 33.80 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5guard_geneaware | 34.92 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 35.35 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 38.24 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5103 | discard | 10.47 | 0.9750 | 0.9750 | 0.9083 |
| synthetic-difficult5x12-seed5103 | likelihood_top5000 | 17.36 | 0.9083 | 0.9083 | 0.8583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_top5000 | 30.72 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_geneaware | 22.91 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1fallback_geneaware | 22.10 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 27.71 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_functionalguard_geneaware | 32.84 | 0.9917 | 0.9917 | 0.9500 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5guard_geneaware | 32.53 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 31.92 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 32.31 | 1.0000 | 1.0000 | 0.9583 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| discard | 7 | 8.85 | 0.9693 | 0.9643 | 0.8854 |
| enhancedgate_functionalguard_geneaware | 7 | 26.64 | 0.9905 | 0.9878 | 0.9485 |
| enhancedgate_geneaware | 7 | 17.38 | 0.9866 | 0.9842 | 0.9446 |
| enhancedgate_kir2dl1_kir2ds5guard_geneaware | 7 | 22.67 | 0.9905 | 0.9866 | 0.9473 |
| enhancedgate_kir2dl1fallback_geneaware | 7 | 17.76 | 0.9881 | 0.9842 | 0.9461 |
| enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 7 | 25.63 | 0.9988 | 0.9961 | 0.9568 |
| enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 7 | 26.82 | 1.0000 | 0.9973 | 0.9580 |
| enhancedgate_kir2dl5guard_geneaware | 7 | 25.88 | 0.9940 | 0.9914 | 0.9509 |
| enhancedgate_top5000 | 7 | 20.92 | 0.9866 | 0.9842 | 0.9446 |
| likelihood_top5000 | 7 | 13.40 | 0.9437 | 0.9426 | 0.8946 |

## Current Candidate Vs Discard

No per-gene 3/5-digit regressions versus discard were detected.

Per-gene 3/5-digit gains versus discard:

| panel | gene | metric | discard F1 | candidate F1 | delta |
|---|---:|---:|---:|---:|---:|
| synthetic-difficult5 | KIR2DS3 | three_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5 | KIR2DS3 | five_digit | 0.8750 | 1.0000 | 0.1250 |
| synthetic-difficult5x12 | KIR2DL5A | three_digit | 0.9744 | 1.0000 | 0.0256 |
| synthetic-difficult5x12 | KIR2DL5A | five_digit | 0.9744 | 1.0000 | 0.0256 |
| synthetic-difficult5x12 | KIR2DL5B | three_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12 | KIR2DL5B | five_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS3 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS5 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12 | KIR2DS5 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DL5A | three_digit | 0.9412 | 1.0000 | 0.0588 |
| synthetic-difficult5x12-seed5101 | KIR2DL5A | five_digit | 0.9412 | 1.0000 | 0.0588 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5101 | KIR2DS3 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | three_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5101 | KIR2DS5 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | three_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DL5B | five_digit | 0.8889 | 1.0000 | 0.1111 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9167 | 0.9583 | 0.0417 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | three_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-difficult5x12-seed5102 | KIR2DS5 | five_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | three_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5103 | KIR2DL5A | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS3 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | three_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-difficult5x12-seed5103 | KIR2DS4 | five_digit | 0.9583 | 1.0000 | 0.0417 |
| synthetic-functional8x6 | KIR2DS1 | three_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-functional8x6 | KIR2DS1 | five_digit | 0.7500 | 1.0000 | 0.2500 |

Remaining current-candidate 3/5-digit errors:

| panel | gene | metric | F1 |
|---|---:|---:|---:|
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9583 |
| synthetic-functional8x6 | KIR2DL1 | five_digit | 0.9167 |

## Decision

`enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware`
extends the previous KIR2DL5/KIR2DS5 candidate with a second
rank-wide targeted unsupported-overcall guard for the exact
`KIR2DS3*0020101` selected allele.
It keeps the strict
`synthetic-functional8x6` KIR2DL1 3-digit regression fixed, keeps
the KIR2DS5 promotion and targeted unsupported guards, keeps the
KIR2DL5 unsupported guard, and keeps the gene-aware top-n runtime
setting.
`likelihood_top5000` alone is not viable because it loses substantial
3/5-digit accuracy on the difficult5x12 seed panels.
KIR2DS3 is now 3-digit perfect on this stress sweep, with only the
seed5102 5-digit suballele row left. KIR2DL1 still has one 5-digit
miss matching discard's remaining error, so next method work should
focus on those two suballele-level misses without reintroducing the
KIR2DL1 3-digit ambiguity-likelihood regression.
