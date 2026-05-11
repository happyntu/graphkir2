# Functional Stress Seed Sweep

This synthetic-first sweep re-runs only the typing and evaluation stages
on existing generated fixtures. It is intended to expose 3-digit and
5-digit functional typing weaknesses before spending time on real data.

Command:

```bash
python benchmarks/scripts/run_functional_stress_sweep.py
```

Summary TSV: `benchmarks/results/functional-stress-sweep/kir2dl1suballele_summary.tsv`
Per-gene TSV: `benchmarks/results/functional-stress-sweep/kir2dl1suballele_per_gene_summary.tsv`

## Aggregate Results

| panel | method | runtime s | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 7.06 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 6.61 | 1.0000 | 1.0000 | 0.9062 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 19.45 | 1.0000 | 0.9896 | 0.9167 |
| synthetic-functional8x6 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 19.21 | 1.0000 | 1.0000 | 0.9271 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 9.14 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 9.28 | 1.0000 | 1.0000 | 1.0000 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 29.21 | 1.0000 | 1.0000 | 0.9667 |
| synthetic-difficult5x12 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 28.71 | 1.0000 | 1.0000 | 0.9667 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 30.63 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5101 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 29.29 | 1.0000 | 1.0000 | 0.9833 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 30.37 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5102 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 31.99 | 1.0000 | 0.9917 | 0.9750 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 33.41 | 1.0000 | 1.0000 | 0.9583 |
| synthetic-difficult5x12-seed5103 | enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 30.85 | 1.0000 | 1.0000 | 0.9583 |

## Mean By Method

| method | panels | mean runtime s | mean 3-digit F1 | mean 5-digit F1 | mean 7-digit F1 |
|---|---:|---:|---:|---:|---:|
| enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_geneaware | 7 | 22.75 | 1.0000 | 0.9973 | 0.9580 |
| enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware | 7 | 22.28 | 1.0000 | 0.9988 | 0.9595 |

## Current Candidate Vs Discard

No per-gene 3/5-digit regressions versus discard were detected.

Per-gene 3/5-digit gains versus discard:

| panel | gene | metric | discard F1 | candidate F1 | delta |
|---|---:|---:|---:|---:|---:|

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
