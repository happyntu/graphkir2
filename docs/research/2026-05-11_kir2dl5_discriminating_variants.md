# KIR2DL5 Discriminating Variant Audit

This report expands the remaining KIR2DL5 failures to genotype-level
truth-only and candidate-only variants using the current likelihood
multi-map weights.

Command:

```bash
python benchmarks/scripts/inspect_kir2dl5_discriminating_variants.py
```

Full TSV: `benchmarks/results/functional-stress-sweep/kir2dl5_discriminating_variants.tsv`

## Sample Summary

| panel | sample | root cause | truth-only variants | candidate-only variants | truth-only exon | candidate-only exon | truth net | candidate net | truth-candidate margin | ambiguous positive ratio | signal |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| synthetic-difficult5x12 | 09 | ab_assignment_mismatch_on_merged_backbone | 0 | 17 | 0 | 9 | 0.000 | -297.000 | 297.000 | 0.000 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5101 | 11 | shared_allele_substitution | 0 | 10 | 0 | 8 | 0.000 | -49.000 | 49.000 | 0.000 | variant_signal_favors_truth |
| synthetic-difficult5x12-seed5103 | 04 | shared_allele_substitution | 0 | 10 | 0 | 8 | 0.000 | -43.000 | 43.000 | 1.000 | variant_signal_favors_truth |

## Direction Counts

| support direction | variants |
|---|---:|
| no_observed_signal | 24 |
| supports_truth | 13 |

## Top Discriminating Variants

| panel | sample | side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | candidate carriers |
|---|---:|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| synthetic-difficult5x12-seed5101 | 11 | candidate_only | hv889 | 1479 | yes | unsupported | supports_truth | 0.000 | 41.000 | -41.000 | 0.000 |  | KIR2DL5A*01201 |
| synthetic-difficult5x12 | 09 | candidate_only | hv816 | 443 | no | unsupported | supports_truth | 0.000 | 41.000 | -41.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv818 | 463 | no | unsupported | supports_truth | 0.000 | 40.000 | -40.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv817 | 450 | no | unsupported | supports_truth | 0.000 | 39.000 | -39.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv823 | 562 | yes | unsupported | supports_truth | 0.000 | 37.000 | -37.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv917 | 3458 | yes | unsupported | supports_truth | 0.000 | 36.000 | -36.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv813 | 393 | no | unsupported | supports_truth | 0.000 | 35.000 | -35.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv820 | 524 | no | unsupported | supports_truth | 0.000 | 35.000 | -35.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12-seed5103 | 04 | candidate_only | hv889 | 1479 | yes | unsupported | supports_truth | 0.000 | 31.000 | -31.000 | 0.000 |  | KIR2DL5A*01201 |
| synthetic-difficult5x12 | 09 | candidate_only | hv807 | 332 | no | unsupported | supports_truth | 0.000 | 18.000 | -18.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv988 | 9315 | yes | unsupported | supports_truth | 0.000 | 16.000 | -16.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12-seed5103 | 04 | candidate_only | hv988 | 9315 | yes | negative_lean | supports_truth | 1.000 | 13.000 | -12.000 | 1.000 |  | KIR2DL5A*01201 |
| synthetic-difficult5x12-seed5101 | 11 | candidate_only | hv988 | 9315 | yes | unsupported | supports_truth | 0.000 | 8.000 | -8.000 | 0.000 |  | KIR2DL5A*01201 |
| synthetic-difficult5x12 | 09 | candidate_only | hv824 | 581 | yes | unobserved | no_observed_signal | 0.000 | 0.000 | 0.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv890 | 1489 | yes | unobserved | no_observed_signal | 0.000 | 0.000 | 0.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv909 | 2536 | yes | unobserved | no_observed_signal | 0.000 | 0.000 | 0.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv927 | 3705 | yes | unobserved | no_observed_signal | 0.000 | 0.000 | 0.000 | 0.000 |  | KIR2DL5B*01301 |
| synthetic-difficult5x12 | 09 | candidate_only | hv946 | 5928 | yes | unobserved | no_observed_signal | 0.000 | 0.000 | 0.000 | 0.000 |  | KIR2DL5B*01301 |

## Interpretation

* KIR2DL5 failures are not caused by missing truth-only variants: in these three cases the truth genotype has no variants absent from the candidate genotype, while the wrong candidate introduces extra candidate-only variants.
* The observed candidate-only variants are mostly unsupported or explicitly contradicted by negative reads. The positive truth-candidate margins mean this discriminating-variant evidence favors the truth genotype despite the current likelihood rank-1 candidate.
* Seed5103 has ambiguous positive signal on a candidate-only exon variant, but the stronger candidate-only negative evidence still favors truth. A guard should therefore penalize unsupported candidate-only variants rather than simply reward any candidate private positive signal.
* Use this report to design a narrow KIR2DL5A/B placement or KIR2DL5A*001-vs-*012 overcall guard before changing any KIR2DS3/KIR2DS5 logic.
