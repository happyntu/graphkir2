# Gene-aware top-n seed validation

This validates targeted typing top-n reduction on four difficult5x12 synthetic
seeds. The comparison checks whether lower non-target allele candidate caps keep
the same functional typing quality as the full `top_n=5000` baseline.

Compared variants:

* `full_top5000`: all genes use `top_n=5000`.
* `targeted_base600`: enhancedgate target genes keep `top_n=5000`; other genes use `allele_base_top_n=600`.
* `geneaware_base600_kir2dl1_1000`: same as `targeted_base600`, with `KIR2DL1:1000`.

Command:

```bash
python benchmarks/scripts/profile_targeted_top_n.py \
  --label synthetic-difficult5x12 \
  --label synthetic-difficult5x12-seed5101 \
  --label synthetic-difficult5x12-seed5102 \
  --label synthetic-difficult5x12-seed5103 \
  --results-dir benchmarks/results/targeted-topn-difficult5x12-seed-validation \
  --output-tsv benchmarks/results/targeted-topn-difficult5x12-seed-validation/summary.tsv \
  --output-md docs/research/2026-05-11_geneaware_top_n_seed_validation.md
```

| panel | variant | runtime s | max RSS MB | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|---:|---:|
| synthetic-difficult5x12 | full_top5000 | 26.00 | 606.9 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | targeted_base600 | 21.90 | 524.2 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | geneaware_base600_kir2dl1_1000 | 19.87 | 530.0 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12-seed5101 | full_top5000 | 28.29 | 568.2 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | targeted_base600 | 23.32 | 519.9 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5101 | geneaware_base600_kir2dl1_1000 | 23.97 | 510.3 | 0.9750 | 0.9750 | 0.9667 |
| synthetic-difficult5x12-seed5102 | full_top5000 | 27.21 | 608.8 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | targeted_base600 | 19.29 | 586.6 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5102 | geneaware_base600_kir2dl1_1000 | 22.40 | 592.6 | 0.9750 | 0.9667 | 0.9417 |
| synthetic-difficult5x12-seed5103 | full_top5000 | 27.76 | 593.5 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | targeted_base600 | 19.97 | 587.0 | 0.9917 | 0.9833 | 0.9417 |
| synthetic-difficult5x12-seed5103 | geneaware_base600_kir2dl1_1000 | 23.02 | 587.5 | 0.9917 | 0.9833 | 0.9417 |

## Interpretation

Both reduced top-n variants preserved the aggregate 3-digit, 5-digit, and
7-digit F1 scores across all four difficult5x12 seed panels.

| variant | mean runtime s | speedup vs full | mean max RSS MB | RSS change vs full | per-gene F1 changes vs full |
|---|---:|---:|---:|---:|---:|
| full_top5000 | 27.32 | 1.00x | 594.4 | 0.0% | baseline |
| targeted_base600 | 21.12 | 1.29x | 554.4 | -6.7% | 0 |
| geneaware_base600_kir2dl1_1000 | 22.32 | 1.22x | 555.1 | -6.6% | 0 |

`targeted_base600` is faster on three of these four difficult5x12 seeds, but it
already regressed `functional8x6` in the earlier top-n profile. The safer default
for cross-panel testing remains `geneaware_base600_kir2dl1_1000` because it
preserved the functional8x6 result while still reducing difficult5x12 runtime.

## Per-gene Regression Check

The check compared each bundle's `per_gene.tsv` against `full_top5000` for every
panel, gene, and 3-digit/5-digit/7-digit F1 metric.

| variant | panels | genes per panel | metric comparisons | regressions | any F1 change |
|---|---:|---:|---:|---:|---:|
| targeted_base600 | 4 | 6 | 72 | 0 | 0 |
| geneaware_base600_kir2dl1_1000 | 4 | 6 | 72 | 0 | 0 |

No per-gene regressions were detected. The difficult5x12 seed validation therefore
supports using reduced non-target allele caps for synthetic runtime iteration,
while keeping the gene-aware override as the safer real-data-facing setting.

Generated artifact paths are listed in the TSV summary.
