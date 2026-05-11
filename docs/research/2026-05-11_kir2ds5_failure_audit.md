# KIR2DS5 Remaining Failure Audit

This report audits the remaining KIR2DS5 functional errors for the
current functional-stress lead and expands each miss into likelihood
rank/support diagnostics plus truth-only vs candidate-only variant evidence.

Command:

```bash
python benchmarks/scripts/inspect_kir2ds5_remaining_failures.py
```

Sample TSV: `benchmarks/results/functional-stress-sweep/kir2ds5_failure_audit.tsv`
Variant TSV: `benchmarks/results/functional-stress-sweep/kir2ds5_discriminating_variants.tsv`

## Summary

| root cause | samples |
|---|---:|

| variant signal | samples |
|---|---:|

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | truth rank | candidate rank | truth gap | candidate gap | truth support | candidate support | truth-only variants | candidate-only variants | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|

## Direction Counts

| support direction | variants |
|---|---:|

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|

## Interpretation

* No KIR2DS5 3-digit or 5-digit remaining errors are present for the current candidate.
* Keep the existing KIR2DS5 promotion guard and targeted `KIR2DS5*027/KIR2DS5*010` unsupported-overcall guard unchanged unless broader stress panels reveal a candidate regression.
