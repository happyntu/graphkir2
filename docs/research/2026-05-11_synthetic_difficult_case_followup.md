# 2026-05-11 Synthetic Difficult-Case Follow-up

## Goal

Check whether the current candidate baseline

* `multi_map_policy = best-only`
* `allele_exon_weight = 2.0`

still helps outside the original `synthetic-functional8` benchmark.

## Datasets

### 1. `synthetic-functional8x6`

Gene panel:

* `KIR2DL1`
* `KIR2DS1`
* `KIR2DL2`
* `KIR2DL3`
* `KIR2DS2`
* `KIR3DL1`
* `KIR3DS1`
* `KIR3DL2`

Setup:

* `6` samples
* `20x` ART depth
* `msa_type = ab`
* exon-only alleles disabled

### 2. `synthetic-difficult5`

Gene panel:

* `KIR2DL4`
* `KIR2DL5`
* `KIR2DS3`
* `KIR2DS4`
* `KIR2DS5`

Setup:

* `4` samples
* `20x` ART depth
* `msa_type = ab`
* exon-only alleles disabled

## Results

### `synthetic-functional8x6`

Discard baseline:

* `3-digit = 0.96875`
* `5-digit = 0.95833`
* `7-digit = 0.85417`

Candidate baseline (`best-only + exon_weight=2.0`):

* `3-digit = 0.95833`
* `5-digit = 0.95833`
* `7-digit = 0.86458`

Interpretation:

* `7-digit` improved slightly
* `5-digit` stayed flat
* `3-digit` regressed slightly

Per-gene shift:

* `KIR2DS1` improved strongly at `3/5/7-digit`
* `KIR2DL1` regressed strongly at `3/5/7-digit`

This means the current candidate baseline redistributes ambiguity rather than
cleanly removing it.

### `synthetic-difficult5`

Discard baseline:

* `3-digit = 0.86111`
* `5-digit = 0.86111`
* `7-digit = 0.77778`

Candidate baseline (`best-only + exon_weight=2.0`):

* `3-digit = 0.86111`
* `5-digit = 0.86111`
* `7-digit = 0.80556`

Interpretation:

* `3-digit` and `5-digit` did not improve
* `7-digit` improved modestly
* `KIR2DL5A/B` still fail completely in this synthetic setup

This panel is useful because it shows the current candidate baseline does not
solve all difficult loci, especially `KIR2DL5`.

## Decision

Do **not** promote `best-only + exon_weight=2.0` to the default global baseline
yet.

Current status:

* it is clearly better than `keep-all`
* it can improve `7-digit`
* it can rescue specific `3/5-digit` failures such as `KIR2DL3` on the original
  `synthetic-functional8`
* but it is not uniformly better on expanded difficult-case panels

## Next Recommended Work

1. Add gene-aware scoring rather than one global exon weight.
2. Focus first on the confusion pairs:
   * `KIR2DL1` vs `KIR2DS1`
   * `KIR2DL2` / `KIR2DL3` / `KIR2DS2`
   * `KIR2DL5A/B` vs nearby short-read-confusable loci
3. Keep `best-only` as the preferred multi-map policy scaffold for now.

## Gene-aware follow-up

Implemented a first gene-aware weighting surface with
`GENE:WEIGHT,GENE:WEIGHT` overrides.

Tested on `synthetic-functional8x6` with:

* `KIR2DS1:2.0,KIR2DL3:2.0`
* `KIR2DS1:1.5,KIR2DL3:2.0`
* `KIR2DS1:1.25,KIR2DL3:2.0`

Observed outcome:

* all three settings produced the same cohort summary as the earlier
  `best-only + exon_weight=2.0` candidate
* `KIR2DS1` stayed improved
* `KIR2DL1` stayed degraded

Interpretation:

* the current ambiguity is not fixed by a simple per-gene exon-weight override
* the next likely step is not more scalar weighting, but a different scoring rule
  for confusing gene pairs or a better multi-map evidence model

## Margin-aware multi-map follow-up

Implemented a `margin` multi-map mode with this rule:

* use `best-only` when the top alignment is clearly ahead
* if `top1 - top2` is within a score gap, keep only the near-tie alignments
  and normalize their weights

Default tested settings:

* `margin_gap = 5.0`
* `margin_scale = 2.0`
* plus a small sweep on `synthetic-functional8x6`:
  * `gap=3.0, scale=2.0`
  * `gap=5.0, scale=1.0`
  * `gap=8.0, scale=2.0`

Observed outcome:

### `synthetic-functional8x6`

All tested margin settings produced the same cohort summary:

* `3-digit = 0.9375`
* `5-digit = 0.9375`
* `7-digit = 0.83333`

Interpretation:

* margin-aware sharing is worse than `best-only + exon_weight=2.0`
* it is also worse than the discard baseline on this panel
* the top1-top2 margin heuristic does not resolve the `KIR2DL1/2DS1` and
  `KIR2DL2/2DL3/2DS2` ambiguity in a useful way

### `synthetic-functional8`

`margin + exon_weight=2.0` dropped back below the earlier
`best-only + exon_weight=2.0` result:

* `3-digit = 0.96875`
* `5-digit = 0.96875`
* `7-digit = 0.90625`

Interpretation:

* margin mode loses the `3/5-digit = 1.0` improvement that
  `best-only + exon_weight=2.0` achieved on the original smaller panel
* this makes it closer to the earlier weighted-like behavior than to the
  current best candidate

### `synthetic-difficult5`

`margin + exon_weight=2.0` also matched the earlier
`best-only + exon_weight=2.0` result:

* `3-digit = 0.86111`
* `5-digit = 0.86111`
* `7-digit = 0.80556`

Interpretation:

* margin mode does not rescue `KIR2DL5A/B`
* it does not provide a meaningful new direction on this panel either

## Updated Decision

Do **not** promote the current `margin` heuristic.

Current status:

* it is reproducible and now available as a benchmark surface
* but it does not improve the target `3-digit/5-digit` functional typing
* and on the harder `synthetic-functional8x6` panel it regresses all three
  allele-resolution metrics

The next likely step should move away from preselecting or pre-sharing mappings
at the alignment stage and instead handle ambiguity deeper inside the typing
likelihood or candidate-selection logic.

## Ambiguity-likelihood follow-up

Implemented a `likelihood` multi-map mode that defers alignment uncertainty into
the allele likelihood:

* keep all multi-mapped alignments with score-derived alignment priors
* collapse alignments with the same physical read name inside the typing model
* score one physical read with a mixture likelihood:
  `P(read | allele) = sum(P(alignment) * P(read evidence | allele, alignment))`
* any non-target alignment mass is treated as neutral for the current gene

This differs from the earlier `weighted` mode. `weighted` uses a weighted log
product across alignments, which effectively asks one allele to explain multiple
candidate alignments. `likelihood` instead treats the alignment as a latent
choice and marginalizes it.

Default tested settings:

* `multi_map_policy = likelihood`
* `allele_exon_weight = 2.0`
* `allele_ambiguity_neutral_prob = 0.999`

### `synthetic-functional8`

`likelihood + exon_weight=2.0`:

* `3-digit = 1.0`
* `5-digit = 1.0`
* `7-digit = 0.90625`

Interpretation:

* matches the best result previously seen with `best-only + exon_weight=2.0`
* avoids the `3/5-digit` regression seen with `weighted` and `margin`

### `synthetic-functional8x6`

`likelihood + exon_weight=2.0`:

* `3-digit = 0.96875`
* `5-digit = 0.96875`
* `7-digit = 0.875`

Comparison:

* discard baseline: `0.96875 / 0.95833 / 0.85417`
* `best-only + exon_weight=2.0`: `0.95833 / 0.95833 / 0.86458`
* `margin + exon_weight=2.0`: `0.9375 / 0.9375 / 0.83333`

Interpretation:

* this is the first tested multi-map enhancement that improves `5-digit` and
  `7-digit` on `synthetic-functional8x6` while preserving the discard baseline
  `3-digit`
* `KIR2DL1` remains the main weak gene, but `KIR2DL3`, `KIR2DS1`, and `KIR2DS2`
  are all clean at `3/5-digit`

Neutral-probability sweep on this panel:

* `0.9`
* `0.99`
* `0.999`
* `0.9999`

All tested settings produced the same cohort-level summary, so the observed
gain appears to come from the mixture-likelihood structure rather than a fragile
neutral-probability constant.

### `synthetic-difficult5`

`likelihood + exon_weight=2.0`:

* `3-digit = 0.86111`
* `5-digit = 0.86111`
* `7-digit = 0.80556`

Interpretation:

* matches `best-only + exon_weight=2.0`
* still does not rescue `KIR2DL5A/B`

## Updated Candidate

Promote `likelihood + exon_weight=2.0` to the current synthetic candidate
baseline for multi-map handling.

It is not a complete solution, but it is better aligned with the functional
typing objective than `weighted`, `margin`, or confusion-pair reranking:

* keeps `synthetic-functional8` at `3/5-digit = 1.0`
* improves `synthetic-functional8x6` over the previous candidate
* does not harm `synthetic-difficult5`

Next work should focus on `KIR2DL1` and `KIR2DL5A/B`, likely by improving
candidate selection or adding gene-specific structural/functional evidence
rather than adding another alignment-stage heuristic.

## Candidate fraction reranking follow-up

The `KIR2DL1` wrong calls in `synthetic-functional8x6` showed a consistent
pattern:

* wrong candidates often had an imbalanced allele fraction near `0.69 / 0.31`
* correct functional candidates existed deeper in the candidate list with more
  balanced fractions

Implemented a candidate selection parameter:

* `select_min_fraction_ratio`
* default legacy-compatible value: `0.5`
* tested candidate value: `0.7`

For CN=2, this changes the selected allele minimum fraction from `0.25` to
`0.35`.

### Results

`likelihood + exon_weight=2.0 + select_min_fraction_ratio=0.7`:

* `synthetic-functional8`: `1.0 / 1.0 / 0.90625`
* `synthetic-functional8x6`: `0.97917 / 0.97917 / 0.88542`
* `synthetic-difficult5`: `0.86111 / 0.86111 / 0.80556`

Interpretation:

* improves `synthetic-functional8x6` over `likelihood + exon_weight=2.0`
* does not regress the original `synthetic-functional8` panel
* does not fix or worsen `KIR2DL5A/B`

Updated current synthetic candidate baseline:

* `multi_map_policy = likelihood`
* `allele_exon_weight = 2.0`
* `allele_select_min_fraction_ratio = 0.7`

## KIR2DL5 backbone normalization follow-up

Inspection showed that the remaining `KIR2DL5A/B` failures were not caused by
the ambiguity likelihood surface. The synthetic CN truth used
`KIR2DL5A*BACKBONE` / `KIR2DL5B*BACKBONE`, while the `ab` HISAT graph stores the
merged locus as `KIR2DL5*BACKBONE`. Typing therefore looked up A/B-specific
backbone buckets that had no reads or variants.

Implemented conservative CN-key normalization in `TypingWithPosNegAllele`:

* keep exact CN keys when the corresponding backbone exists
* only map `KIR2DL5A*BACKBONE` and `KIR2DL5B*BACKBONE` to
  `KIR2DL5*BACKBONE` when the exact key is absent and the merged key is present
* sum A/B CN values before allele typing on the merged KIR2DL5 graph

Re-run candidate:

* `multi_map_policy = likelihood`
* `allele_exon_weight = 2.0`
* `allele_select_min_fraction_ratio = 0.7`

Updated `synthetic-difficult5` result:

* `3-digit = 0.975`
* `5-digit = 0.975`
* `7-digit = 0.925`
* `KIR2DL5A = 1.0 / 1.0 / 1.0`
* `KIR2DL5B = 1.0 / 1.0 / 1.0`

Remaining difficult5 failures are now:

* `KIR2DS3` at `3/5/7-digit`
* `KIR2DS4` at `7-digit`

## KIR2DS3 / KIR2DS5 cross-gene ambiguity follow-up

The remaining `KIR2DS3` functional error in `synthetic-difficult5` is sample
`02`:

* truth: `KIR2DS3*0011201 + KIR2DS3*0010301`
* prediction: `KIR2DS3*0011201 + KIR2DS3*0020101`

Diagnostic result:

* `KIR2DS3*0010301` is present in the graph candidate set
* the false-positive `KIR2DS3*0020101` variants are mostly supported by
  synthetic reads whose read names start with `KIR2DS5*0020105` or
  `KIR2DS5*0020109`
* these reads enter `KIR2DS3` through cross-gene ambiguous alignments with
  target weight around `0.5`

Tested but not promoted:

* `filter_minor=True`: rescues neither the target case nor global performance;
  it regresses multiple genes
* individual-support threshold sweeps: not enough to fix `KIR2DS3`
* ambiguity target-weight powers: no cohort-level improvement on difficult5
* neutralizing all cross-gene ambiguous reads: does not fix `KIR2DS3` and
  regresses `7-digit`

Useful finding:

* increasing typing `top_n` from `600` to `5000` improves
  `synthetic-functional8x6` from `0.97917 / 0.97917 / 0.88542` to
  `0.98958 / 0.98958 / 0.90625`
* `synthetic-functional8` remains `1.0 / 1.0 / 0.90625`
* `synthetic-difficult5` remains `0.975 / 0.975 / 0.925`

Updated current synthetic candidate baseline:

* `multi_map_policy = likelihood`
* `allele_exon_weight = 2.0`
* `allele_select_min_fraction_ratio = 0.7`
* `top_n = 5000`

Next method work should target cross-gene ambiguous read assignment between
`KIR2DS3` and `KIR2DS5`, not generic candidate filtering.

## Dosage-aware reranking follow-up

Tested a minimal dosage-aware reranking idea for CN=2 candidates:

* score each allele pair as a mixture:
  `P(read | pair) = f * P(read | allele1) + (1 - f) * P(read | allele2)`
* grid-search `f` to estimate allele dosage/fraction
* add optional balance prior to penalize extreme allele fractions
* combine cross-gene neutralization with individual-support filtering

Result on the difficult `KIR2DS3` sample:

* plain mixture reranking penalizes the false `KIR2DS3*0020101` candidate, but
  ranks other near-neighbor alleles such as `KIR2DS3*0220103`,
  `KIR2DS3*00109`, and `KIR2DS3*00108` above the truth
* strong dosage-balance priors move the truth closer, but still do not make it
  the top functional call
* cross-gene neutralization plus individual-support filtering does not fix the
  `KIR2DS3` functional error and regresses `7-digit` on the synthetic panels

Cohort-level ablation outcome:

* cross-gene neutralization alone:
  * `synthetic-functional8`: `1.0 / 1.0 / 0.84375`
  * `synthetic-functional8x6`: `0.97917 / 0.97917 / 0.86458`
  * `synthetic-difficult5`: `0.975 / 0.975 / 0.85`
* cross-gene neutralization plus individual-support filters did not improve
  `synthetic-difficult5` functional F1 beyond the current `0.975 / 0.975`

Decision:

Do not promote these reranking/filtering heuristics. The remaining failure is a
joint assignment problem: `KIR2DS5` reads that ambiguously align to `KIR2DS3`
must be assigned across genes before per-gene allele selection, rather than
counted independently inside both genes.

## Targeted private-support follow-up

Tested a narrower version of the joint-assignment idea:

* keep the current `likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 +
  top_n=5000` candidate surface
* neutralize only the explicit cross-gene ambiguity group `KIR2DS3/KIR2DS5`
* rerank only `KIR2DS3` with private-variant support
* use `private_support_lambda = 10`
* use `private_support_window = 50`

This differs from the earlier broad neutralization tests. It does not remove all
cross-gene ambiguous reads globally, and it does not apply private-support
reranking to every gene.

### Results

`synthetic-functional8`:

* top5000 baseline: `1.0 / 1.0 / 0.90625`
* targeted private-support: `1.0 / 1.0 / 0.90625`

`synthetic-functional8x6`:

* top5000 baseline: `0.98958 / 0.98958 / 0.90625`
* targeted private-support: `0.98958 / 0.98958 / 0.90625`

`synthetic-difficult5`:

* top5000 baseline: `0.975 / 0.975 / 0.925`
* targeted private-support, broad neutralization: `1.0 / 1.0 / 0.9`
* targeted private-support, directional neutralization: `1.0 / 1.0 / 0.95`

Per-gene `synthetic-difficult5` result for the targeted candidate:

* `KIR2DL4 = 1.0 / 1.0 / 1.0`
* `KIR2DL5A = 1.0 / 1.0 / 1.0`
* `KIR2DL5B = 1.0 / 1.0 / 1.0`
* `KIR2DS3 = 1.0 / 1.0 / 1.0`
* `KIR2DS4 = 1.0 / 1.0 / 0.75`
* `KIR2DS5 = 1.0 / 1.0 / 1.0`

The difficult sample `02` now calls the expected
`KIR2DS3*0010301 + KIR2DS3*0011201`. The remaining tradeoff is that one
`KIR2DS4` 7-digit call remains unresolved.

Directional-neutralization finding:

* broad neutralization rescued `KIR2DS3` but also removed useful `KIR2DS5`
  ambiguous evidence, lowering `KIR2DS5` 7-digit
* directional neutralization zeros only the target `KIR2DS3` side of the
  `KIR2DS3/KIR2DS5` ambiguous evidence
* `KIR2DS5` keeps its own evidence and returns to `7-digit = 1.0`

## Current interpretation

Keep two candidate labels:

* balanced synthetic baseline:
  `likelihood + exon_weight=2.0 + min_fraction_ratio=0.7 + top_n=5000`
* functional-target candidate:
  balanced baseline plus targeted `KIR2DS3` private-support reranking and
  directional `KIR2DS3/KIR2DS5` ambiguity neutralization

For the stated `graphkir2` objective, the functional-target candidate is the
more relevant lead because it is the first tested method to reach
`synthetic-difficult5` `3/5-digit = 1.0` without regressing the other synthetic
functional panels. It also improves `synthetic-difficult5` `7-digit` above the
balanced baseline after directional neutralization.

Implementation follow-up:

* the private-support and cross-gene neutralization helpers now live under
  `src/graphkir2/typing/private_support.py`
* the relevant tuning knobs are carried by `TypingConfig`, the CLI, benchmark
  presets, and typing plans
* `benchmarks/configs/synthetic-difficult5-functional-target-kir2ds3-private.json`
  captures the current difficult5 functional-target settings
* running `rerun_typing_private_support.py` with that preset and no private
  support flags reproduces the previous hand-specified
  `private_support_kir2ds3_group_neutral.allele.tsv` output
