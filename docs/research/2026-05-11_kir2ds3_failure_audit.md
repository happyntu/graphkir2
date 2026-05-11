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
| shared_with_discard | 1 |

| variant signal | samples |
|---|---:|
| variant_signal_favors_candidate | 1 |
| variant_signal_favors_truth | 1 |

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | truth rank | candidate rank | truth gap | candidate gap | truth support | candidate support | truth-only variants | candidate-only variants | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| synthetic-difficult5x12 | 03 | shared_with_discard | KIR2DS3*011_KIR2DS3*016 | KIR2DS3*0020101_KIR2DS3*011 | KIR2DS3*0020101_KIR2DS3*011 | KIR2DS3*0020101_KIR2DS3*011 | 692 | 1 | 383.156 | 0.000 | 4.400 | -126.950 | 1 | 38 | 964.000 | 0.672 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5102 | 00 | all_methods_disagree_or_shifted | KIR2DS3*0010306_KIR2DS3*0011302 | KIR2DS3*0010301_KIR2DS3*0010306 | KIR2DS3*00108_KIR2DS3*0011302 | KIR2DS3*0010306_KIR2DS3*0020101 | 575 | 426 | 160.310 | 154.603 | -13.750 | 6.900 | 1 | 0 | -37.500 | 0.000 | variant_signal_favors_candidate |

## Direction Counts

| support direction | variants |
|---|---:|
| no_observed_signal | 1 |
| supports_candidate | 3 |
| supports_truth | 36 |

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| synthetic-difficult5x12 | 03 | candidate_only | hv1409 | 12515 | no | unsupported | supports_truth | 0.000 | 64.000 | -64.000 | 0.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1398 | 11648 | no | unsupported | supports_truth | 0.000 | 60.500 | -60.500 | 0.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1411 | 12631 | no | negative_lean | supports_truth | 11.000 | 71.000 | -60.000 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1395 | 11395 | no | unsupported | supports_truth | 0.000 | 59.500 | -59.500 | 0.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1442 | 14731 | no | negative_lean | supports_truth | 10.500 | 66.000 | -55.500 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1394 | 11381 | no | unsupported | supports_truth | 0.000 | 54.000 | -54.000 | 0.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1392 | 11316 | no | unsupported | supports_truth | 0.000 | 49.000 | -49.000 | 0.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1396 | 11476 | no | negative_lean | supports_truth | 12.500 | 51.000 | -38.500 | 0.440 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12-seed5102 | 00 | truth_only | hv1432 | 14454 | yes | unsupported | supports_candidate | 0.000 | 37.500 | -37.500 | 0.000 | KIR2DS3*0011302 |  |
| synthetic-difficult5x12 | 03 | candidate_only | hv1397 | 11621 | no | negative_lean | supports_truth | 16.000 | 49.000 | -33.000 | 0.438 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1412 | 12753 | no | negative_lean | supports_truth | 21.000 | 52.000 | -31.000 | 0.571 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1408 | 12472 | no | negative_lean | supports_truth | 23.000 | 50.000 | -27.000 | 0.391 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1399 | 11714 | no | negative_lean | supports_truth | 19.000 | 44.000 | -25.000 | 0.526 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1415 | 13170 | no | negative_lean | supports_truth | 21.000 | 46.000 | -25.000 | 0.952 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1418 | 13448 | no | negative_lean | supports_truth | 16.500 | 41.500 | -25.000 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1404 | 11883 | no | negative_lean | supports_truth | 17.500 | 42.000 | -24.500 | 0.600 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1413 | 13054 | no | negative_lean | supports_truth | 17.500 | 41.000 | -23.500 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1400 | 11727 | no | negative_lean | supports_truth | 21.000 | 44.000 | -23.000 | 0.571 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1410 | 12622 | no | negative_lean | supports_truth | 30.500 | 53.000 | -22.500 | 0.377 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1416 | 13333 | no | negative_lean | supports_truth | 13.500 | 36.000 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1424 | 13862 | no | negative_lean | supports_truth | 18.500 | 41.000 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1426 | 13912 | no | negative_lean | supports_truth | 16.000 | 38.500 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1407 | 12419 | no | negative_lean | supports_truth | 20.500 | 42.000 | -21.500 | 0.463 |  | KIR2DS3*0020101 |
| synthetic-difficult5x12 | 03 | candidate_only | hv1403 | 11826 | no | negative_lean | supports_truth | 18.000 | 38.500 | -20.500 | 0.778 |  | KIR2DS3*0020101 |

## Interpretation

* Variant evidence favors truth in 1 sample(s), candidate in 1 sample(s), and is mixed in 0 sample(s).
* Candidate-regression KIR2DS3 samples: 0.
* Unsupported candidate-only `KIR2DS3*002` variants appear in 5 row(s).
* Unsupported truth-only variants appear in 1 row(s), so the residual set is not uniformly truth-rescuable.
* Treat `shared_with_discard` rows separately from all-method shifted rows; shared rows need new evidence, not rollback of current guards.
* Any next KIR2DS3 method should be narrow and preserve the current KIR2DS5/KIR2DL5 gains.
* Do not add a KIR2DS3 guard from aggregate F1 alone; the sample-level variant signal is not uniformly truth-favoring.
