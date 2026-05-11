# KIR2DS3 Rank-Failure Audit

This report focuses on the `synthetic-difficult5x12` sample `03` miss where
`KIR2DS3*011_KIR2DS3*016` is the truth but current likelihood-based typing
selects `KIR2DS3*011_KIR2DS3*0020101`.

Command:

```bash
python benchmarks/scripts/inspect_kir2ds3_rank_failure.py
```

Rank TSV: `benchmarks/results/functional-stress-sweep/kir2ds3_rank_failure.tsv`
Variant TSV: `benchmarks/results/functional-stress-sweep/kir2ds3_rank_failure_variants.tsv`

## Sample Context

| field | value |
|---|---|
| config | benchmarks/configs/synthetic-difficult5x12-conditional-kir2ds3-enhancedgate.json |
| truth | KIR2DS3*011_KIR2DS3*016 |
| current | KIR2DS3*0020101_KIR2DS3*011 |
| discard | KIR2DS3*0020101_KIR2DS3*011 |
| likelihood | KIR2DS3*0020101_KIR2DS3*011 |
| copy_number | 2 |
| top_n | 5000 |
| likelihood_reads | 2494 |
| likelihood_informative_reads | 1883 |
| variants | 1536 |
| exon_variants | 37 |
| positive_variant_ids | 42 |
| negative_variant_ids | 119 |

## Key Rank Rows

| rank | kind | genotype | 3-digit | 5-digit | value | gap | private support | truth-only variants | candidate-only variants | candidate-only unsupported | truth-candidate margin | signal |
|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | top;current;discard;likelihood | KIR2DS3*0020101_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -126.950 | 1 | 38 | 5 | 964.000 | variant_signal_favors_truth |
| 2 | top | KIR2DS3*0020103_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -140.650 | 1 | 38 | 6 | 1001.000 | variant_signal_favors_truth |
| 3 | top | KIR2DS3*0020105_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -154.350 | 1 | 39 | 7 | 1038.000 | variant_signal_favors_truth |
| 4 | top | KIR2DS3*011_KIR2DS3*024 | KIR2DS3*011_KIR2DS3*024 | KIR2DS3*011_KIR2DS3*024 | -384.430 | 0.000 | -154.450 | 1 | 39 | 7 | 1039.000 | variant_signal_favors_truth |
| 5 | top | KIR2DS3*0030101N_KIR2DS3*011 | KIR2DS3*003_KIR2DS3*011 | KIR2DS3*00301_KIR2DS3*011 | -384.430 | 0.000 | -154.850 | 1 | 39 | 7 | 1043.000 | variant_signal_favors_truth |
| 6 | top | KIR2DS3*0020107_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -154.250 | 1 | 39 | 7 | 1037.000 | variant_signal_favors_truth |
| 7 | top | KIR2DS3*0020109_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -167.450 | 1 | 40 | 8 | 1069.000 | variant_signal_favors_truth |
| 8 | top | KIR2DS3*0020104_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -168.550 | 1 | 40 | 8 | 1080.000 | variant_signal_favors_truth |
| 9 | top | KIR2DS3*0020102_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -168.650 | 1 | 41 | 8 | 1081.000 | variant_signal_favors_truth |
| 10 | top | KIR2DS3*0020106_KIR2DS3*011 | KIR2DS3*002_KIR2DS3*011 | KIR2DS3*00201_KIR2DS3*011 | -384.430 | 0.000 | -170.150 | 1 | 40 | 8 | 1096.000 | variant_signal_favors_truth |
| 692 | truth | KIR2DS3*011_KIR2DS3*016 | KIR2DS3*011_KIR2DS3*016 | KIR2DS3*011_KIR2DS3*016 | -767.586 | 383.156 | 4.400 | 0 | 0 | 0 | 0.000 | variant_signal_mixed |

## Selected-vs-Truth Variant Summary

| truth-only variants | selected-only variants | truth-only exon | selected-only exon | truth-only unsupported | selected-only unsupported | truth net | selected net | truth-selected margin | ambiguous positive ratio | signal |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 38 | 1 | 1 | 0 | 5 | -3.000 | -967.000 | 964.000 | 0.672 | variant_signal_favors_truth |

## Top Selected-vs-Truth Variants

| side | variant | pos | exon | class | direction | positive | negative | net | ambiguous positive ratio | truth carriers | selected carriers |
|---|---|---:|---|---|---|---:|---:|---:|---:|---|---|
| candidate_only | hv1409 | 12515 | no | unsupported | supports_truth | 0.000 | 64.000 | -64.000 | 0.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1398 | 11648 | no | unsupported | supports_truth | 0.000 | 60.500 | -60.500 | 0.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1411 | 12631 | no | negative_lean | supports_truth | 11.000 | 71.000 | -60.000 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1395 | 11395 | no | unsupported | supports_truth | 0.000 | 59.500 | -59.500 | 0.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1442 | 14731 | no | negative_lean | supports_truth | 10.500 | 66.000 | -55.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1394 | 11381 | no | unsupported | supports_truth | 0.000 | 54.000 | -54.000 | 0.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1392 | 11316 | no | unsupported | supports_truth | 0.000 | 49.000 | -49.000 | 0.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1396 | 11476 | no | negative_lean | supports_truth | 12.500 | 51.000 | -38.500 | 0.440 |  | KIR2DS3*0020101 |
| candidate_only | hv1397 | 11621 | no | negative_lean | supports_truth | 16.000 | 49.000 | -33.000 | 0.438 |  | KIR2DS3*0020101 |
| candidate_only | hv1412 | 12753 | no | negative_lean | supports_truth | 21.000 | 52.000 | -31.000 | 0.571 |  | KIR2DS3*0020101 |
| candidate_only | hv1408 | 12472 | no | negative_lean | supports_truth | 23.000 | 50.000 | -27.000 | 0.391 |  | KIR2DS3*0020101 |
| candidate_only | hv1399 | 11714 | no | negative_lean | supports_truth | 19.000 | 44.000 | -25.000 | 0.526 |  | KIR2DS3*0020101 |
| candidate_only | hv1415 | 13170 | no | negative_lean | supports_truth | 21.000 | 46.000 | -25.000 | 0.952 |  | KIR2DS3*0020101 |
| candidate_only | hv1418 | 13448 | no | negative_lean | supports_truth | 16.500 | 41.500 | -25.000 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1404 | 11883 | no | negative_lean | supports_truth | 17.500 | 42.000 | -24.500 | 0.600 |  | KIR2DS3*0020101 |
| candidate_only | hv1413 | 13054 | no | negative_lean | supports_truth | 17.500 | 41.000 | -23.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1400 | 11727 | no | negative_lean | supports_truth | 21.000 | 44.000 | -23.000 | 0.571 |  | KIR2DS3*0020101 |
| candidate_only | hv1410 | 12622 | no | negative_lean | supports_truth | 30.500 | 53.000 | -22.500 | 0.377 |  | KIR2DS3*0020101 |
| candidate_only | hv1416 | 13333 | no | negative_lean | supports_truth | 13.500 | 36.000 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1424 | 13862 | no | negative_lean | supports_truth | 18.500 | 41.000 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1426 | 13912 | no | negative_lean | supports_truth | 16.000 | 38.500 | -22.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1407 | 12419 | no | negative_lean | supports_truth | 20.500 | 42.000 | -21.500 | 0.463 |  | KIR2DS3*0020101 |
| candidate_only | hv1403 | 11826 | no | negative_lean | supports_truth | 18.000 | 38.500 | -20.500 | 0.778 |  | KIR2DS3*0020101 |
| candidate_only | hv1401 | 11760 | no | negative_lean | supports_truth | 19.000 | 39.000 | -20.000 | 0.737 |  | KIR2DS3*0020101 |
| candidate_only | hv1425 | 13870 | no | negative_lean | supports_truth | 19.000 | 38.000 | -19.000 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1428 | 13960 | no | negative_lean | supports_truth | 16.000 | 33.000 | -17.000 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1402 | 11803 | no | negative_lean | supports_truth | 18.500 | 35.000 | -16.500 | 0.838 |  | KIR2DS3*0020101 |
| candidate_only | hv1405 | 12062 | no | negative_lean | supports_truth | 21.500 | 37.000 | -15.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1422 | 13796 | yes | negative_lean | supports_truth | 18.500 | 33.000 | -14.500 | 1.000 |  | KIR2DS3*0020101 |
| candidate_only | hv1427 | 13954 | no | negative_lean | supports_truth | 17.500 | 32.000 | -14.500 | 1.000 |  | KIR2DS3*0020101 |

## Interpretation

* Truth rank is 692 with likelihood gap 383.156 and private support 4.400.
* Current selected rank is 1 with private support -126.950.
* The selected-vs-truth variant margin is 964.000; positive values favor truth at variant-evidence level.
* The miss is therefore not a near-window selection issue; `KIR2DS3*016` is too far below rank 1 for the existing top-window gates to rescue.
* Any next experiment should test a very narrow `KIR2DS3*0020101` overcall penalty or pre-likelihood candidate pruning, then verify it does not regress the seed5102 sample where variant evidence favors the non-truth candidate.
