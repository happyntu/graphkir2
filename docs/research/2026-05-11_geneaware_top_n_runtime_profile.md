# Targeted top-n runtime profile

This profiles `rerun_typing_private_support.py` on existing synthetic fixtures.
`full_top5000` keeps every gene at `top_n=5000`; targeted variants keep
the enhancedgate target genes at `top_n=5000` and use lower `base_top_n`
values for non-target genes; gene-aware variants override individual genes.

| panel | variant | runtime s | max RSS MB | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|---:|---:|
| synthetic-functional8x6 | full_top5000 | 25.17 | 1145.7 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | targeted_base600 | 17.92 | 1154.8 | 0.9792 | 0.9792 | 0.8854 |
| synthetic-functional8x6 | geneaware_base600_kir2dl1_1000 | 14.75 | 1146.5 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-difficult5x12 | full_top5000 | 25.32 | 606.4 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | targeted_base600 | 21.84 | 523.9 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | geneaware_base600_kir2dl1_1000 | 18.96 | 531.8 | 0.9750 | 0.9750 | 0.9500 |

## Interpretation

* `synthetic-functional8x6`: fastest non-regressing targeted variant was `geneaware_base600_kir2dl1_1000` at 14.75s vs 25.17s (1.71x speedup).
* `synthetic-functional8x6`: `targeted_base600` changed aggregate F1; inspect its bundle before using it as a default.
* `synthetic-difficult5x12`: fastest non-regressing targeted variant was `geneaware_base600_kir2dl1_1000` at 18.96s vs 25.32s (1.34x speedup).

Generated artifact paths are listed in the TSV summary.
