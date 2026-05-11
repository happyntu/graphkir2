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
| ab_assignment_mismatch_on_merged_backbone | 1 |
| shared_allele_substitution | 2 |

## Detail

| panel | sample | root cause | truth | candidate | truth CN | candidate CN | normalized CN | likelihood selected | truth rank | truth gap | truth support | candidate support | discard selected | discard truth rank |
|---|---:|---|---|---|---|---|---:|---|---:|---:|---:|---:|---|---:|
| synthetic-difficult5x12 | 09 | ab_assignment_mismatch_on_merged_backbone | KIR2DL5A*0010901_KIR2DL5A*028 | KIR2DL5A*028_KIR2DL5B*01301 | KIR2DL5A=2;KIR2DL5B=0;merged=2 | KIR2DL5A=1;KIR2DL5B=1;merged=2 | 2 | KIR2DL5A*028_KIR2DL5B*01301 | 3 | 7.499 | 5.500 | -41.400 | KIR2DL5A*028_KIR2DL5B*01303 | 3 |
| synthetic-difficult5x12-seed5101 | 11 | shared_allele_substitution | KIR2DL5A*00107_KIR2DL5A*00107 | KIR2DL5A*00107_KIR2DL5A*01201 | KIR2DL5A=2;KIR2DL5B=0;merged=2 | KIR2DL5A=2;KIR2DL5B=0;merged=2 | 2 | KIR2DL5A*00107_KIR2DL5A*01201 | 15 | 20.997 | 71.200 | 35.700 | KIR2DL5A*00107_KIR2DL5A*01201 | 15 |
| synthetic-difficult5x12-seed5103 | 04 | shared_allele_substitution | KIR2DL5A*0010101_KIR2DL5A*0230101 | KIR2DL5A*01201_KIR2DL5A*0230101 | KIR2DL5A=2;KIR2DL5B=0;merged=2 | KIR2DL5A=2;KIR2DL5B=0;merged=2 | 2 | KIR2DL5A*01201_KIR2DL5A*0230101 | 11 | 14.998 | -0.400 | 19.500 | KIR2DL5A*01201_KIR2DL5A*0230101 | 11 |

## Interpretation

* The `synthetic-difficult5x12` sample 09 failure is not merged KIR2DL5 CN loss: the truth has KIR2DL5A=2/KIR2DL5B=0 and normalized merged CN=2, while the candidate places one copy on KIR2DL5B. This is an A/B placement problem on the merged KIR2DL5 backbone.
* The seed5101 sample 11 and seed5103 sample 04 failures keep the KIR2DL5A copy count but substitute KIR2DL5A*001 with KIR2DL5A*012. Because the candidate agrees with discard and likelihood here, this is shared allele discrimination, not an enhancedgate regression.
* All KIR2DL5 truths are still inside the current top-n search (`rank 3`, `15`, and `11`), so the immediate blocker is likelihood/reranking, not candidate truncation.
* Private-support evidence is mixed: it favors truth for sample 09 and seed5101 sample 11, but favors the current candidate for seed5103 sample 04. A broad KIR2DL5 private-support rescue would therefore need a strict condition rather than a global switch.
* Next method work should inspect KIR2DL5A*001-vs-*012 and KIR2DL5A-vs-KIR2DL5B discriminating variants before changing KIR2DS3/KIR2DS5 guards; the current KIR2DS3/KIR2DS4/KIR2DS5 gains should remain isolated.
