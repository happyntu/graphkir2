# KIR2DL1 Remaining Failure Audit

This report audits KIR2DL1 functional-stress rows for the current
lead. When KIR2DL1 rows remain, it compares likelihood evidence against
discard evidence because the KIR2DL1 rescue path is discard-evidence
aware and must preserve the 3-digit functional class.

Command:

```bash
python benchmarks/scripts/inspect_kir2dl1_remaining_failures.py
```

Sample TSV: `benchmarks/results/functional-stress-sweep/kir2dl1_failure_audit.tsv`
Variant TSV: `benchmarks/results/functional-stress-sweep/kir2dl1_discriminating_variants.tsv`

## Summary

| root cause | samples |
|---|---:|

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | likelihood truth rank | likelihood candidate rank | discard truth rank | discard candidate rank | likelihood truth support | likelihood candidate support | discard truth support | discard candidate support |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|

## Variant Evidence By Mode

| evidence mode | truth-only variants | candidate-only variants | truth-only exon | candidate-only exon | truth-only unsupported | candidate-only unsupported | truth net | candidate net | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|

## Direction Counts

| evidence mode / support direction | variants |
|---|---:|

## Top Discriminating Variants

| evidence mode | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|

## Interpretation

* No KIR2DL1 3-digit or 5-digit remaining errors are present for the current candidate.
* Keep the KIR2DL1 functional fallback and exact 3-digit suballele guard narrow unless broader stress panels reveal a candidate regression.
