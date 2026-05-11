# KIR2DL1 Remaining Failure Audit

This report audits the remaining KIR2DL1 5-digit suballele miss for
the current functional-stress lead. It compares likelihood evidence
against discard evidence because the current candidate is produced by
the KIR2DL1 functional fallback and matches discard at 5-digit scope.

Command:

```bash
python benchmarks/scripts/inspect_kir2dl1_remaining_failures.py
```

Sample TSV: `benchmarks/results/functional-stress-sweep/kir2dl1_failure_audit.tsv`
Variant TSV: `benchmarks/results/functional-stress-sweep/kir2dl1_discriminating_variants.tsv`

## Summary

| root cause | samples |
|---|---:|
| shared_with_discard | 1 |

## Sample Detail

| panel | sample | root cause | truth | candidate | discard | likelihood selected | likelihood truth rank | likelihood candidate rank | discard truth rank | discard candidate rank | likelihood truth support | likelihood candidate support | discard truth support | discard candidate support |
|---|---:|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| synthetic-functional8x6 | 05 | shared_with_discard | KIR2DL1*0030229_KIR2DL1*0030242 | KIR2DL1*0030242_KIR2DL1*00303 | KIR2DL1*0030242_KIR2DL1*00303 | KIR2DL1*0030205_KIR2DL1*0040110 | not in top-n | not in top-n | 277 | 144 | -43.550 | -49.033 | -8.200 | -19.400 |

## Variant Evidence By Mode

| evidence mode | truth-only variants | candidate-only variants | truth-only exon | candidate-only exon | truth-only unsupported | candidate-only unsupported | truth net | candidate net | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| likelihood | 1 | 9 | 0 | 1 | 1 | 8 | -54.500 | -393.500 | 339.000 | 0.818 | variant_signal_favors_truth |
| discard | 1 | 9 | 0 | 1 | 0 | 9 | -17.000 | -303.000 | 286.000 | 0.000 | variant_signal_favors_truth |

## Direction Counts

| evidence mode / support direction | variants |
|---|---:|
| discard:supports_candidate | 1 |
| discard:supports_truth | 9 |
| likelihood:supports_candidate | 1 |
| likelihood:supports_truth | 9 |

## Top Discriminating Variants

| evidence mode | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| likelihood | candidate_only | hv174 | 7433 | no | unsupported | supports_truth | 0.000 | 59.500 | -59.500 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv193 | 8118 | no | unsupported | supports_truth | 0.000 | 58.500 | -58.500 | 0.000 |  | KIR2DL1*00303 |
| likelihood | truth_only | hv61 | 3384 | no | unsupported | supports_candidate | 0.000 | 54.500 | -54.500 | 0.000 | KIR2DL1*0030229 |  |
| likelihood | candidate_only | hv130 | 6026 | no | unsupported | supports_truth | 0.000 | 52.000 | -52.000 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv257 | 10917 | no | unsupported | supports_truth | 0.000 | 50.000 | -50.000 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv355 | 13456 | yes | unsupported | supports_truth | 0.000 | 48.500 | -48.500 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv193 | 8118 | no | unsupported | supports_truth | 0.000 | 46.000 | -46.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv56 | 3195 | no | unsupported | supports_truth | 0.000 | 45.000 | -45.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv174 | 7433 | no | unsupported | supports_truth | 0.000 | 45.000 | -45.000 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv239 | 10124 | no | unsupported | supports_truth | 0.000 | 44.000 | -44.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv257 | 10917 | no | unsupported | supports_truth | 0.000 | 36.000 | -36.000 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv344 | 12952 | no | unsupported | supports_truth | 0.000 | 35.500 | -35.500 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv130 | 6026 | no | unsupported | supports_truth | 0.000 | 34.000 | -34.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv239 | 10124 | no | unsupported | supports_truth | 0.000 | 33.000 | -33.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv355 | 13456 | yes | unsupported | supports_truth | 0.000 | 31.000 | -31.000 | 0.000 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv56 | 3195 | no | negative_lean | supports_truth | 22.000 | 45.500 | -23.500 | 0.818 |  | KIR2DL1*00303 |
| likelihood | candidate_only | hv91 | 4774 | no | unsupported | supports_truth | 0.000 | 22.000 | -22.000 | 0.000 |  | KIR2DL1*00303 |
| discard | candidate_only | hv91 | 4774 | no | unsupported | supports_truth | 0.000 | 22.000 | -22.000 | 0.000 |  | KIR2DL1*00303 |
| discard | truth_only | hv61 | 3384 | no | negative_lean | supports_candidate | 17.000 | 34.000 | -17.000 | 0.000 | KIR2DL1*0030229 |  |
| discard | candidate_only | hv344 | 12952 | no | unsupported | supports_truth | 0.000 | 11.000 | -11.000 | 0.000 |  | KIR2DL1*00303 |

## Interpretation

* The remaining KIR2DL1 row is shared with discard, so it is not introduced by the current graphkir2 guards.
* Because the error is 5-digit-only and KIR2DL1 3-digit is already fixed, any next method must preserve the `KIR2DL1*003` functional class.
* Treat a KIR2DL1 suballele rescue as unsafe unless truth-only 5-digit evidence is supported in both likelihood and discard evidence modes, or unless a tie-break can be shown not to alter 3-digit behavior.
