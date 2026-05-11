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
| synthetic-functional8 | discard | 3.67 | 1.0000 | 1.0000 | 0.8438 |
| synthetic-functional8 | likelihood_top5000 | 6.42 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_top5000 | 7.47 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_geneaware | 6.66 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1fallback_geneaware | 6.04 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 7.42 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_functionalguard_geneaware | 7.46 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5guard_geneaware | 7.43 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 7.84 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 7.64 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 6.94 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | discard | 7.50 | 0.9688 | 0.9583 | 0.8542 |
| synthetic-functional8x6 | likelihood_top5000 | 14.35 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_top5000 | 18.42 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_geneaware | 15.73 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_kir2dl1fallback_geneaware | 17.14 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 21.10 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_functionalguard_geneaware | 21.81 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5guard_geneaware | 22.33 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 22.79 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 30.71 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 24.39 | 1.0000 | 1.0000 | 0.9271 |
| synthetic-difficult5 | discard | 3.98 | 0.9750 | 0.9750 | 0.8750 |
| synthetic-difficult5 | likelihood_top5000 | 6.01 | 0.9750 | 0.9750 | 0.9250 |
| synthetic-difficult5 | enhancedgate_top5000 | 8.08 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_geneaware | 6.85 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1fallback_geneaware | 6.61 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 8.07 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_functionalguard_geneaware | 9.58 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5guard_geneaware | 9.53 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 9.68 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 10.37 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 10.55 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | discard | 8.76 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12 | likelihood_top5000 | 14.52 | 0.9250 | 0.9250 | 0.9083 |
| synthetic-difficult5x12 | enhancedgate_top5000 | 25.71 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_geneaware | 19.92 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1fallback_geneaware | 19.65 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 25.91 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_functionalguard_geneaware | 30.26 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5guard_geneaware | 30.04 | 0.9833 | 0.9833 | 0.9500 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 32.32 | 0.9917 | 0.9917 | 0.9583 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 33.13 | 1.0000 | 1.0000 | 0.9667 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 40.16 | 1.0000 | 1.0000 | 0.9667 |
| synthetic-difficult5x12-seed5101 | discard | 13.55 | 0.9667 | 0.9583 | 0.9250 |
| synthetic-difficult5x12-seed5101 | likelihood_top5000 | 16.00 | 0.9000 | 0.9000 | 0.8833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_top5000 | 25.80 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_geneaware | 24.49 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1fallback_geneaware | 24.69 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 30.35 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5101 | enhancedgate_functionalguard_geneaware | 36.98 | 0.9917 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5guard_geneaware | 39.54 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 36.16 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 40.70 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 49.33 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5102 | discard | 21.41 | 0.9333 | 0.9250 | 0.8667 |
| synthetic-difficult5x12-seed5102 | likelihood_top5000 | 22.54 | 0.9083 | 0.9000 | 0.8750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_top5000 | 36.78 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_geneaware | 28.69 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1fallback_geneaware | 28.75 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 28.90 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_functionalguard_geneaware | 40.31 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5guard_geneaware | 33.68 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 35.82 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 35.61 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 37.71 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5103 | discard | 15.78 | 0.9750 | 0.9750 | 0.9083 |
| synthetic-difficult5x12-seed5103 | likelihood_top5000 | 15.78 | 0.9083 | 0.9083 | 0.8583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_top5000 | 28.28 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_geneaware | 21.06 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1fallback_geneaware | 20.94 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl1_kir2ds5guard_geneaware | 26.78 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | enhancedgate_functionalguard_geneaware | 37.16 | 0.9917 | 0.9917 | 0.9500 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5guard_geneaware | 41.99 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 37.92 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 40.95 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 35.56 | 1.0000 | 1.0000 | 0.9583 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| discard | 7 | 10.66 | 0.9693 | 0.9643 | 0.8854 |
| enhancedgate_functionalguard_geneaware | 7 | 26.22 | 0.9905 | 0.9878 | 0.9485 |
| enhancedgate_geneaware | 7 | 17.63 | 0.9866 | 0.9842 | 0.9446 |
| enhancedgate_kir2dl1_kir2ds5guard_geneaware | 7 | 21.22 | 0.9905 | 0.9866 | 0.9473 |
| enhancedgate_kir2dl1fallback_geneaware | 7 | 17.69 | 0.9881 | 0.9842 | 0.9461 |
| enhancedgate_kir2dl5_kir2ds5unsupported_geneaware | 7 | 26.08 | 0.9988 | 0.9961 | 0.9568 |
| enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 7 | 28.44 | 1.0000 | 0.9973 | 0.9580 |
| enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 7 | 29.23 | 1.0000 | 0.9988 | 0.9595 |
| enhancedgate_kir2dl5guard_geneaware | 7 | 26.36 | 0.9940 | 0.9914 | 0.9509 |
| enhancedgate_top5000 | 7 | 21.51 | 0.9866 | 0.9842 | 0.9446 |
| likelihood_top5000 | 7 | 13.66 | 0.9437 | 0.9426 | 0.8946 |

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
| synthetic-functional8x6 | KIR2DL1 | five_digit | 0.9167 | 1.0000 | 0.0833 |
| synthetic-functional8x6 | KIR2DS1 | three_digit | 0.7500 | 1.0000 | 0.2500 |
| synthetic-functional8x6 | KIR2DS1 | five_digit | 0.7500 | 1.0000 | 0.2500 |

Remaining current-candidate 3/5-digit errors:

| panel | gene | metric | F1 |
|---|---:|---:|---:|
| synthetic-difficult5x12-seed5102 | KIR2DS3 | five_digit | 0.9583 |

## Decision

`enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware`
extends the KIR2DS3 rank-wide candidate with a discard-evidence
unsupported-overcall guard for the exact `KIR2DL1*00303` selected
allele. The guard preserves the selected 3-digit allele multiset,
so it can correct the 5-digit KIR2DL1 suballele miss without
reintroducing the earlier KIR2DL1 3-digit ambiguity regression.
It also keeps the KIR2DS5 promotion and targeted unsupported guards,
the KIR2DL5 unsupported guard, and the gene-aware top-n runtime setting.
`likelihood_top5000` alone is not viable because it loses substantial
3/5-digit accuracy on the difficult5x12 seed panels.
Current remaining 3/5-digit functional error is limited to the
KIR2DS3 seed5102 5-digit suballele row; next method work should
inspect that unsupported truth-only suballele rather than broadening
the existing guards.
