# KIR2DL5 Remaining Failure Audit

This report audits the remaining KIR2DL5A/B functional errors for the
current functional-stress lead at the merged KIR2DL5 typing-backbone scope.

Command:

```bash
python benchmarks/scripts/inspect_kir2dl5_remaining_failures.py
```

Full TSV: `benchmarks/results/functional-stress-sweep/kir2dl5_failure_audit.tsv`

## Summary

| root cause | samples |
|---|---:|

## Detail

| panel | sample | root cause | truth | candidate | truth CN | candidate CN | normalized CN | likelihood selected | truth rank | truth gap | truth support | candidate support | discard selected | discard truth rank |
|---|---:|---|---|---|---|---|---:|---|---:|---:|---:|---:|---|---:|

## Interpretation

* No KIR2DL5A/B 3-digit or 5-digit remaining errors are present for the current candidate.
* The KIR2DL5-only unsupported-overcall guard resolved the previous A/B placement and KIR2DL5A*001-vs-*012 overcall cases without leaving a KIR2DL5 functional miss in this stress sweep.
* Keep the guard narrow: it should continue to require multiple unsupported selected-only variants before replacing a likelihood-selected KIR2DL5 call.
