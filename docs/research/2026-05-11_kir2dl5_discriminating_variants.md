# KIR2DL5 Discriminating Variant Audit

This report expands the remaining KIR2DL5 failures to genotype-level
truth-only and candidate-only variants using the current likelihood
multi-map weights.

Command:

```bash
python benchmarks/scripts/inspect_kir2dl5_discriminating_variants.py
```

Full TSV: `benchmarks/results/functional-stress-sweep/kir2dl5_discriminating_variants.tsv`

## Sample Summary

| panel | sample | root cause | truth-only variants | candidate-only variants | truth-only exon | candidate-only exon | truth net | candidate net | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|

## Direction Counts

| support direction | variants |
|---|---:|

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|

## Interpretation

* No KIR2DL5 discriminating-variant rows are present because the current candidate has no remaining KIR2DL5A/B functional errors in this stress sweep.
* The unsupported-overcall guard should remain constrained to KIR2DL5 and to cases with multiple unsupported selected-only variants.
