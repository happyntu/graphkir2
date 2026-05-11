# KIR2DS3 Remaining Failure Audit

This report audits the remaining KIR2DS3 functional errors for the
current functional-stress lead and expands each miss into likelihood
rank/support diagnostics plus truth-only vs candidate-only variant evidence.

Command:

```bash
python benchmarks/scripts/inspect_kir2ds3_remaining_failures.py
```

Sample TSV: `benchmarks/results/functional-stress-sweep/kir2ds3_failure_audit.tsv`
Variant TSV: `benchmarks/results/functional-stress-sweep/kir2ds3_discriminating_variants.tsv`

## Summary

| root cause | samples |
|---|---:|
| all_methods_disagree_or_shifted | 1 |

| variant signal | samples |
|---|---:|
| variant_signal_favors_candidate | 1 |

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | truth rank | candidate rank | truth gap | candidate gap | truth support | candidate support | truth-only variants | candidate-only variants | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| synthetic-difficult5x12-seed5102 | 00 | all_methods_disagree_or_shifted | KIR2DS3*0010306_KIR2DS3*0011302 | KIR2DS3*0010301_KIR2DS3*0010306 | KIR2DS3*00108_KIR2DS3*0011302 | KIR2DS3*0010306_KIR2DS3*0020101 | 575 | 426 | 160.310 | 154.603 | -13.750 | 6.900 | 1 | 0 | -37.500 | 0.000 | variant_signal_favors_candidate |

## Direction Counts

| support direction | variants |
|---|---:|
| supports_candidate | 1 |

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| synthetic-difficult5x12-seed5102 | 00 | truth_only | hv1432 | 14454 | yes | unsupported | supports_candidate | 0.000 | 37.500 | -37.500 | 0.000 | KIR2DS3*0011302 |  |

## Interpretation

* Variant evidence favors truth in 0 sample(s), candidate in 1 sample(s), and is mixed in 0 sample(s).
* Candidate-regression KIR2DS3 samples: 0.
* Unsupported candidate-only `KIR2DS3*002` variants appear in 0 row(s).
* Unsupported truth-only variants appear in 1 row(s), so the residual set is not uniformly truth-rescuable.
* Treat `shared_with_discard` rows separately from all-method shifted rows; shared rows need new evidence, not rollback of current guards.
* Any next KIR2DS3 method should be narrow and preserve the current KIR2DS5/KIR2DL5 gains.
* Do not add a KIR2DS3 guard from aggregate F1 alone; the sample-level variant signal is not uniformly truth-favoring.
