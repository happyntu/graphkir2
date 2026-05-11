# Remaining Functional Error Triage

This report expands the functional stress sweep into sample-level
3-digit and 5-digit wrong calls for the current candidate.

Command:

```bash
python benchmarks/scripts/inspect_remaining_functional_errors.py
```

Candidate method: `enhancedgate_functionalguard_geneaware`
Full TSV: `benchmarks/results/functional-stress-sweep/remaining_functional_errors.tsv`

## Summary

| gene | resolution | cause hint | wrong-call rows |
|---|---:|---:|---:|
| KIR2DL1 | 5 | shared_with_discard | 1 |
| KIR2DL5A | 3 | copy_number_count_mismatch | 1 |
| KIR2DL5A | 3 | shared_with_discard | 2 |
| KIR2DL5A | 5 | copy_number_count_mismatch | 1 |
| KIR2DL5A | 5 | shared_with_discard | 2 |
| KIR2DL5B | 3 | copy_number_count_mismatch | 1 |
| KIR2DL5B | 5 | copy_number_count_mismatch | 1 |
| KIR2DS3 | 3 | shared_with_discard | 1 |
| KIR2DS3 | 5 | all_methods_disagree_or_shifted | 1 |
| KIR2DS3 | 5 | shared_with_discard | 1 |
| KIR2DS5 | 3 | shared_with_discard | 2 |
| KIR2DS5 | 3 | unresolved_likelihood_pattern | 2 |
| KIR2DS5 | 5 | shared_with_discard | 2 |
| KIR2DS5 | 5 | unresolved_likelihood_pattern | 2 |

## Detail

| panel | sample | gene | res | cause hint | truth | candidate | missing | extra | discard | likelihood | enhancedgate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| synthetic-difficult5x12 | 03 | KIR2DS3 | 3 | shared_with_discard | KIR2DS3*011_KIR2DS3*016 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*016 | KIR2DS3*002 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 |
| synthetic-difficult5x12 | 03 | KIR2DS3 | 5 | shared_with_discard | KIR2DS3*011_KIR2DS3*016 | KIR2DS3*00201_KIR2DS3*011 | KIR2DS3*016 | KIR2DS3*00201 | KIR2DS3*00201_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 |
| synthetic-difficult5x12 | 05 | KIR2DS5 | 3 | shared_with_discard | KIR2DS5*002_KIR2DS5*024 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*024 | KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 |
| synthetic-difficult5x12 | 05 | KIR2DS5 | 5 | shared_with_discard | KIR2DS5*00205_KIR2DS5*024 | KIR2DS5*00205_KIR2DS5*02701 | KIR2DS5*024 | KIR2DS5*02701 | KIR2DS5*00205_KIR2DS5*02701 | KIR2DS5*00205_KIR2DS5*02701 | KIR2DS5*00205_KIR2DS5*02701 |
| synthetic-difficult5x12 | 09 | KIR2DL5A | 3 | copy_number_count_mismatch | KIR2DL5A*001_KIR2DL5A*028 | KIR2DL5A*028 | KIR2DL5A*001 |  | KIR2DL5A*028 | KIR2DL5A*028 | KIR2DL5A*028 |
| synthetic-difficult5x12 | 09 | KIR2DL5A | 5 | copy_number_count_mismatch | KIR2DL5A*00109_KIR2DL5A*028 | KIR2DL5A*028 | KIR2DL5A*00109 |  | KIR2DL5A*028 | KIR2DL5A*028 | KIR2DL5A*028 |
| synthetic-difficult5x12 | 09 | KIR2DL5B | 3 | copy_number_count_mismatch |  | KIR2DL5B*013 |  | KIR2DL5B*013 | KIR2DL5B*013 | KIR2DL5B*013 | KIR2DL5B*013 |
| synthetic-difficult5x12 | 09 | KIR2DL5B | 5 | copy_number_count_mismatch |  | KIR2DL5B*01301 |  | KIR2DL5B*01301 | KIR2DL5B*01303 | KIR2DL5B*01301 | KIR2DL5B*01301 |
| synthetic-difficult5x12-seed5101 | 11 | KIR2DL5A | 3 | shared_with_discard | KIR2DL5A*001_KIR2DL5A*001 | KIR2DL5A*001_KIR2DL5A*012 | KIR2DL5A*001 | KIR2DL5A*012 | KIR2DL5A*001_KIR2DL5A*012 | KIR2DL5A*001_KIR2DL5A*012 | KIR2DL5A*001_KIR2DL5A*012 |
| synthetic-difficult5x12-seed5101 | 11 | KIR2DL5A | 5 | shared_with_discard | KIR2DL5A*00107_KIR2DL5A*00107 | KIR2DL5A*00107_KIR2DL5A*01201 | KIR2DL5A*00107 | KIR2DL5A*01201 | KIR2DL5A*00107_KIR2DL5A*01201 | KIR2DL5A*00107_KIR2DL5A*01201 | KIR2DL5A*00107_KIR2DL5A*01201 |
| synthetic-difficult5x12-seed5102 | 00 | KIR2DS3 | 5 | all_methods_disagree_or_shifted | KIR2DS3*00103_KIR2DS3*00113 | KIR2DS3*00103_KIR2DS3*00103 | KIR2DS3*00113 | KIR2DS3*00103 | KIR2DS3*00108_KIR2DS3*00113 | KIR2DS3*00103_KIR2DS3*00201 | KIR2DS3*00103_KIR2DS3*00103 |
| synthetic-difficult5x12-seed5102 | 05 | KIR2DS5 | 3 | unresolved_likelihood_pattern | KIR2DS5*002_KIR2DS5*021 | KIR2DS5*021_KIR2DS5*027 | KIR2DS5*002 | KIR2DS5*027 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*021_KIR2DS5*027 | KIR2DS5*021_KIR2DS5*027 |
| synthetic-difficult5x12-seed5102 | 05 | KIR2DS5 | 5 | unresolved_likelihood_pattern | KIR2DS5*00201_KIR2DS5*02101 | KIR2DS5*02101_KIR2DS5*02701 | KIR2DS5*00201 | KIR2DS5*02701 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*02101_KIR2DS5*02701 | KIR2DS5*02101_KIR2DS5*02701 |
| synthetic-difficult5x12-seed5102 | 07 | KIR2DS5 | 3 | shared_with_discard | KIR2DS5*002_KIR2DS5*024 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*024 | KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 |
| synthetic-difficult5x12-seed5102 | 07 | KIR2DS5 | 5 | shared_with_discard | KIR2DS5*00201_KIR2DS5*024 | KIR2DS5*00201_KIR2DS5*02701 | KIR2DS5*024 | KIR2DS5*02701 | KIR2DS5*00201_KIR2DS5*02701 | KIR2DS5*00201_KIR2DS5*02701 | KIR2DS5*00201_KIR2DS5*02701 |
| synthetic-difficult5x12-seed5102 | 08 | KIR2DS5 | 3 | unresolved_likelihood_pattern | KIR2DS5*002_KIR2DS5*024 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*024 | KIR2DS5*027 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*002_KIR2DS5*027 | KIR2DS5*002_KIR2DS5*027 |
| synthetic-difficult5x12-seed5102 | 08 | KIR2DS5 | 5 | unresolved_likelihood_pattern | KIR2DS5*00201_KIR2DS5*024 | KIR2DS5*00201_KIR2DS5*02701 | KIR2DS5*024 | KIR2DS5*02701 | KIR2DS5*034_KIR2DS5*034 | KIR2DS5*00201_KIR2DS5*02701 | KIR2DS5*00201_KIR2DS5*02701 |
| synthetic-difficult5x12-seed5103 | 04 | KIR2DL5A | 3 | shared_with_discard | KIR2DL5A*001_KIR2DL5A*023 | KIR2DL5A*012_KIR2DL5A*023 | KIR2DL5A*001 | KIR2DL5A*012 | KIR2DL5A*012_KIR2DL5A*023 | KIR2DL5A*012_KIR2DL5A*023 | KIR2DL5A*012_KIR2DL5A*023 |
| synthetic-difficult5x12-seed5103 | 04 | KIR2DL5A | 5 | shared_with_discard | KIR2DL5A*00101_KIR2DL5A*02301 | KIR2DL5A*01201_KIR2DL5A*02301 | KIR2DL5A*00101 | KIR2DL5A*01201 | KIR2DL5A*01201_KIR2DL5A*02301 | KIR2DL5A*01201_KIR2DL5A*02301 | KIR2DL5A*01201_KIR2DL5A*02301 |
| synthetic-functional8x6 | 05 | KIR2DL1 | 5 | shared_with_discard | KIR2DL1*00302_KIR2DL1*00302 | KIR2DL1*00302_KIR2DL1*00303 | KIR2DL1*00302 | KIR2DL1*00303 | KIR2DL1*00302_KIR2DL1*00303 | KIR2DL1*00201_KIR2DL1*00302 | KIR2DL1*00302_KIR2DL1*00401 |

## Interpretation

* `shared_with_discard` means the current candidate did not introduce the error; fixing it probably requires new gene-specific evidence handling rather than undoing enhancedgate.
* `candidate_regression` is the highest-priority blocker because discard is correct and the candidate is wrong. Current candidate regressions: none.
* If candidate regressions are `none`, next work should target shared or unresolved failure patterns rather than adding broader regression guards.
* Current remaining errors are concentrated in KIR2DL5A/B, KIR2DS3, KIR2DS5, and one KIR2DL1 5-digit case.
* The KIR2DL1 3-digit regression is absent here; only a 5-digit suballele miss remains, matching discard behavior.

## Recommended Next Method Work

* If KIR2DS5 candidate regressions remain, tighten the KIR2DS5*027 promotion guard rather than broadening KIR2DS3/KIR2DS5 neutralization.
* If KIR2DS3 candidate regressions remain, inspect suballele-specific private support before changing the broader KIR2DS3 rescue gate.
* Keep KIR2DL5A/B separate from the KIR2DS3/KIR2DS5 gate work because its errors include copy-number or A/B placement mismatches shared by all current methods.
