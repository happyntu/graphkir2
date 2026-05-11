# Targeted top-n runtime profile

Goal: quantify whether `allele_base_top_n` can reduce enhancedgate typing cost
without changing functional accuracy.

Setup:

* `top_n = 5000` for enhancedgate target genes
* non-target genes use `base_top_n = 0, 600, 1000, 2000, 3000`
* existing synthetic fixtures, no remapping or data regeneration
* runtime/RSS measured around `rerun_typing_private_support.py`

## Results

| panel | base_top_n | runtime s | max RSS MB | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|---:|---:|
| synthetic-difficult5x12 | 0 | 24.90 | 606.3 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | 600 | 21.37 | 531.2 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | 1000 | 19.47 | 531.8 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | 2000 | 24.26 | 541.5 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-difficult5x12 | 3000 | 24.22 | 606.0 | 0.9750 | 0.9750 | 0.9500 |
| synthetic-functional8x6 | 0 | 20.34 | 1167.5 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | 600 | 13.68 | 1146.9 | 0.9792 | 0.9792 | 0.8854 |
| synthetic-functional8x6 | 1000 | 17.45 | 1155.7 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | 2000 | 14.24 | 1176.0 | 0.9896 | 0.9896 | 0.9062 |
| synthetic-functional8x6 | 3000 | 18.74 | 1161.4 | 0.9896 | 0.9896 | 0.9062 |

## Findings

* `base_top_n = 600` is too aggressive for the functional panel: the aggregate
  `synthetic-functional8x6` score drops from `0.9896 / 0.9896 / 0.9062` to
  `0.9792 / 0.9792 / 0.8854`.
* The `base_top_n = 600` functional regression is localized to `KIR2DL1`: 3/5
  digit F1 drops from `0.9167` to `0.8333`, and 7-digit F1 drops from `0.8333`
  to `0.6667`.
* `base_top_n = 1000` is the smallest tested value that recovered the full
  `top_n=5000` aggregate F1 on both profiled panels.
* Runtime is noisy at this scale, but `base_top_n = 1000` reduced runtime and
  peak RSS on `synthetic-difficult5x12`, and preserved functional panel
  accuracy.

## Decision

Use `allele_base_top_n = 1000` as the full-gene smoke / real-data sanity OOM
guard. Keep `allele_base_top_n = 0` for reduced synthetic lead configs so the
current synthetic benchmark numbers remain comparable to previous runs.

Generated TSV artifacts:

* `benchmarks/results/targeted-topn-difficult5x12-sweep/summary.tsv`
* `benchmarks/results/targeted-topn-functional8x6-sweep/summary.tsv`
