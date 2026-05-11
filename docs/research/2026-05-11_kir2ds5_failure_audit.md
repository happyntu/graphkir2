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
| shared_with_discard | 2 |
| unresolved_likelihood_pattern | 2 |

| variant signal | samples |
|---|---:|
| variant_signal_favors_truth | 4 |

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | truth rank | candidate rank | truth gap | candidate gap | truth support | candidate support | truth-only variants | candidate-only variants | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| synthetic-difficult5x12 | 05 | shared_with_discard | KIR2DS5*0020503_KIR2DS5*024 | KIR2DS5*0020503_KIR2DS5*0270102 | KIR2DS5*0020503_KIR2DS5*0270102 | KIR2DS5*0020503_KIR2DS5*0270102 | 2 | 1 | 4.544 | 0.000 | 13.250 | 11.750 | 1 | 1 | 19.000 | 1.000 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5102 | 05 | unresolved_likelihood_pattern | KIR2DS5*0020129_KIR2DS5*0210101 | KIR2DS5*0210102_KIR2DS5*0270102 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*0210102_KIR2DS5*0270102 | 11 | 1 | 9.128 | 0.000 | 3.150 | -23.100 | 1 | 3 | 72.500 | 1.000 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5102 | 07 | shared_with_discard | KIR2DS5*0020116_KIR2DS5*024 | KIR2DS5*0020116_KIR2DS5*0270102 | KIR2DS5*0020116_KIR2DS5*0270102 | KIR2DS5*0020116_KIR2DS5*0270102 | 2 | 1 | 5.393 | 0.000 | 16.400 | 8.000 | 1 | 1 | 51.000 | 0.852 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5102 | 08 | unresolved_likelihood_pattern | KIR2DS5*0020113_KIR2DS5*024 | KIR2DS5*0020113_KIR2DS5*0270102 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*0020113_KIR2DS5*0270102 | 19 | 1 | 24.826 | 0.000 | -3.350 | -4.850 | 1 | 1 | 29.000 | 1.000 | variant_signal_favors_truth |

## Direction Counts

| support direction | variants |
|---|---:|
| mixed | 1 |
| supports_candidate | 3 |
| supports_truth | 6 |

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| synthetic-difficult5x12-seed5102 | 07 | candidate_only | hv1797 | 4309 | yes | unsupported | supports_truth | 0.000 | 51.000 | -51.000 | 0.000 |  | KIR2DS5*0270102 |
| synthetic-difficult5x12-seed5102 | 08 | candidate_only | hv1797 | 4309 | yes | unsupported | supports_truth | 0.000 | 45.000 | -45.000 | 0.000 |  | KIR2DS5*0270102 |
| synthetic-difficult5x12-seed5102 | 05 | candidate_only | hv1797 | 4309 | yes | unsupported | supports_truth | 0.000 | 42.000 | -42.000 | 0.000 |  | KIR2DS5*0270102 |
| synthetic-difficult5x12 | 05 | candidate_only | hv1797 | 4309 | yes | unsupported | supports_truth | 0.000 | 37.000 | -37.000 | 0.000 |  | KIR2DS5*0270102 |
| synthetic-difficult5x12-seed5102 | 05 | candidate_only | hv1750 | 1392 | no | unsupported | supports_truth | 0.000 | 24.000 | -24.000 | 0.000 |  | KIR2DS5*0210102 |
| synthetic-difficult5x12-seed5102 | 05 | candidate_only | hv1946 | 15054 | no | unsupported | supports_truth | 0.000 | 19.500 | -19.500 | 0.000 |  | KIR2DS5*0210102 |
| synthetic-difficult5x12 | 05 | truth_only | hv1925 | 13962 | yes | negative_lean | supports_candidate | 13.500 | 31.500 | -18.000 | 1.000 | KIR2DS5*024 |  |
| synthetic-difficult5x12-seed5102 | 08 | truth_only | hv1925 | 13962 | yes | negative_lean | supports_candidate | 9.000 | 25.000 | -16.000 | 1.000 | KIR2DS5*024 |  |
| synthetic-difficult5x12-seed5102 | 05 | truth_only | hv1943 | 14731 | no | negative_lean | supports_candidate | 11.000 | 24.000 | -13.000 | 1.000 | KIR2DS5*0020129 |  |
| synthetic-difficult5x12-seed5102 | 07 | truth_only | hv1925 | 13962 | yes | supported | mixed | 13.500 | 13.500 | 0.000 | 0.852 | KIR2DS5*024 |  |

## Interpretation

* Variant evidence favors truth in 4 sample(s), candidate in 0 sample(s), and is mixed in 0 sample(s).
* Candidate-only unsupported variants appear in 6 row(s); 4 row(s) are carried by the selected `KIR2DS5*027` allele.
* `shared_with_discard` means the current candidate is not introducing the miss; fixing those cases requires new evidence, not just undoing the functional guard.
* `unresolved_likelihood_pattern` means likelihood itself prefers the current candidate over the truth while discard is also wrong, so broad KIR2DS5 rollback would be unsafe unless variant evidence is strongly truth-favoring.
* A narrow follow-up guard may be worth testing, but only against these sample-level variant margins and with KIR2DS3/KIR2DS4 gain preservation checked.
