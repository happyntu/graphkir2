# 2026-05-11 Synthetic Method Comparison

## Scope

This report compares the current synthetic multi-map candidates across the
`synthetic-functional8`, `synthetic-functional8x6`, and `synthetic-difficult5`
panels. It is generated from existing benchmark bundles and can be rebuilt with:

```bash
python benchmarks/scripts/summarize_synthetic_methods.py
```

## Method Comparison

| panel | method | 3-digit F1 | 5-digit F1 | 7-digit F1 |
|---|---:|---:|---:|---:|
| synthetic-functional8 | discard | 1 | 1 | 0.84375 |
| synthetic-functional8 | best-only+exon2 | 1 | 1 | 0.90625 |
| synthetic-functional8 | margin+exon2 | 0.96875 | 0.96875 | 0.90625 |
| synthetic-functional8 | likelihood+exon2 | 1 | 1 | 0.90625 |
| synthetic-functional8 | likelihood+exon2+minfrac0.7+top5000 | 1 | 1 | 0.90625 |
| synthetic-functional8 | functional-target+kir2ds3-private | 1 | 1 | 0.90625 |
| synthetic-functional8x6 | discard | 0.96875 | 0.95833 | 0.85417 |
| synthetic-functional8x6 | best-only+exon2 | 0.95833 | 0.95833 | 0.86458 |
| synthetic-functional8x6 | margin+exon2 | 0.9375 | 0.9375 | 0.83333 |
| synthetic-functional8x6 | likelihood+exon2 | 0.96875 | 0.96875 | 0.875 |
| synthetic-functional8x6 | likelihood+exon2+minfrac0.7+top5000 | 0.98958 | 0.98958 | 0.90625 |
| synthetic-functional8x6 | functional-target+kir2ds3-private | 0.98958 | 0.98958 | 0.90625 |
| synthetic-difficult5 | discard | 0.86111 | 0.86111 | 0.77778 |
| synthetic-difficult5 | best-only+exon2 | 0.86111 | 0.86111 | 0.80556 |
| synthetic-difficult5 | margin+exon2 | 0.86111 | 0.86111 | 0.80556 |
| synthetic-difficult5 | likelihood+exon2 | 0.86111 | 0.86111 | 0.80556 |
| synthetic-difficult5 | likelihood+exon2+minfrac0.7+top5000 | 0.975 | 0.975 | 0.925 |
| synthetic-difficult5 | functional-target+kir2ds3-private | 1 | 1 | 1 |

## Candidate Wrong-Call Summary

The wrong-call table uses the current candidate baseline:
`likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 + top_n=5000`.

Full comparison TSV: `benchmarks/results/synthetic-method-comparison.tsv`
Full wrong-call TSV: `benchmarks/results/synthetic-likelihood-wrong-calls.tsv`

| panel | gene | resolution | wrong-call rows |
|---|---:|---:|---:|
| synthetic-difficult5 | KIR2DS3 | 3 | 1 |
| synthetic-difficult5 | KIR2DS3 | 5 | 1 |
| synthetic-difficult5 | KIR2DS3 | 7 | 1 |
| synthetic-difficult5 | KIR2DS4 | 7 | 2 |
| synthetic-functional8 | KIR2DL2 | 7 | 1 |
| synthetic-functional8 | KIR2DL3 | 7 | 1 |
| synthetic-functional8 | KIR2DS2 | 7 | 1 |
| synthetic-functional8x6 | KIR2DL1 | 3 | 1 |
| synthetic-functional8x6 | KIR2DL1 | 5 | 1 |
| synthetic-functional8x6 | KIR2DL1 | 7 | 1 |
| synthetic-functional8x6 | KIR2DL2 | 7 | 2 |
| synthetic-functional8x6 | KIR2DL3 | 7 | 1 |
| synthetic-functional8x6 | KIR2DS1 | 7 | 1 |
| synthetic-functional8x6 | KIR2DS2 | 7 | 1 |
| synthetic-functional8x6 | KIR3DS1 | 7 | 1 |

## Current Decision

`likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 + top_n=5000` is the current synthetic candidate baseline.
It is the first tested multi-map method that improves `5-digit` and `7-digit`
on `synthetic-functional8x6` while preserving the discard baseline `3-digit`.
Raising `top_n` from `600` to `5000` improves `synthetic-functional8x6`
without regressing `synthetic-functional8` or `synthetic-difficult5`.
On `synthetic-difficult5`, the same candidate now also normalizes
`KIR2DL5A/B` CN hints onto the merged `KIR2DL5*BACKBONE` graph key,
which removes the previous missing-call failure mode.

A separate functional-target candidate adds targeted `KIR2DS3` private-support
reranking plus directional `KIR2DS3/KIR2DS5` cross-gene neutralization.
It preserves `synthetic-functional8` and `synthetic-functional8x6` at
the same `3/5-digit` scores as the top5000 baseline and raises
`synthetic-difficult5` to `3/5/7-digit = 1.0` when combined with
the `KIR2DS4` same-5-digit highest-suffix tie-break. Treat it as the
current functional target.

Remaining failure focus:

* balanced baseline: `KIR2DL1` on `synthetic-functional8x6`
* balanced baseline: `KIR2DS3` functional typing on `synthetic-difficult5`
* functional target: no remaining synthetic-difficult5 error in this panel
