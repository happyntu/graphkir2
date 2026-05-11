# Remaining Functional Error Triage

This report expands the functional stress sweep into sample-level
3-digit and 5-digit wrong calls for the current candidate.

Command:

```bash
python benchmarks/scripts/inspect_remaining_functional_errors.py
```

Candidate method: `enhancedgate_kir2dl5_kir2ds5unsupported_kir2ds3rankwide_kir2dl1suballele_geneaware`
Full TSV: `benchmarks/results/functional-stress-sweep/kir2dl1suballele_remaining_functional_errors.tsv`

## Summary

| gene | resolution | cause hint | wrong-call rows |
|---|---:|---:|---:|
| KIR2DS3 | 5 | all_methods_disagree_or_shifted | 1 |

## Detail

| panel | sample | gene | res | cause hint | truth | candidate | missing | extra | discard | likelihood | enhancedgate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| synthetic-difficult5x12-seed5102 | 00 | KIR2DS3 | 5 | all_methods_disagree_or_shifted | KIR2DS3*00103_KIR2DS3*00113 | KIR2DS3*00103_KIR2DS3*00103 | KIR2DS3*00113 | KIR2DS3*00103 |  |  |  |

## Interpretation

* `shared_with_discard` means the current candidate did not introduce the error; fixing it probably requires new gene-specific evidence handling rather than undoing enhancedgate.
* `candidate_regression` is the highest-priority blocker because discard is correct and the candidate is wrong. Current candidate regressions: none.
* If candidate regressions are `none`, next work should target shared or unresolved failure patterns rather than adding broader regression guards.
* Current remaining functional errors are in: KIR2DS3.
* KIR2DL1 is absent from the current 3/5-digit remaining rows; keep the suballele guard narrow unless broader panels reveal a regression.

## Recommended Next Method Work

* KIR2DS5 has no remaining 3/5-digit functional rows for the current candidate; keep the targeted KIR2DS5 guard narrow unless broader panels reveal a regression.
* Inspect the remaining KIR2DS3 rows at suballele/private-support level before changing the broader KIR2DS3 rescue gate.
* Do not broaden the KIR2DL1 suballele guard unless new failures preserve the selected 3-digit genotype multiset.
* Keep any future KIR2DL5A/B work separate from the KIR2DS3/KIR2DS5 gate work.
