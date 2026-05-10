# GraphKir2 Functional Typing Improvement Plan

## Objective

The working method target for `graphkir2` is:

* outperform `Geny` on `3-digit` and `5-digit` functional typing

This target is stricter than preserving the current Graph-KIR strengths in copy number
and `7-digit`. Those remain important guardrails, but they are not the only success
criteria for the refactor.

## Why this needs a separate plan

The legacy Graph-KIR story is strongest in:

* copy number
* `7-digit` full-resolution typing
* runtime

That is not enough for the new method goal. If the refactor keeps optimizing only for
those axes, it may still fail the practical question:

* can short-read functional typing beat `Geny`?

## Primary hypotheses

### 1. Multi-mapped read discard hurts functional typing too early

Current legacy behavior discards multi-mapped reads aggressively. That may be acceptable
for reducing false support in some loci, but it likely removes informative evidence
needed to separate functional groups at `3-digit` and `5-digit`.

Expected intervention:

* compare `discard` vs `best-only` vs `weighted`

Expected benefit:

* better read evidence retention for functionally distinct but highly homologous alleles

### 2. Current typing likelihood is probably too coarse for functional typing

The current likelihood surface appears more aligned with global variant matching than
with explicitly winning the functional-typing task.

Expected intervention:

* separate coding-region and non-coding-region evidence weights
* increase influence of functional-region evidence when evaluating `3-digit/5-digit`
* keep `7-digit` as a separate report target, not the only optimization objective

Expected benefit:

* better discrimination on coding-relevant allele groups

### 3. Confidence-aware downgrade may improve practical `5-digit` accuracy

If a call is weak at full resolution, forcing a `7-digit` answer can still harm the
functional-typing result through unstable upstream ranking.

Expected intervention:

* explicit downgrade path:
  * `7-digit`
  * `5-digit`
  * warning / unresolved

Expected benefit:

* better effective `5-digit` performance and more trustworthy output behavior

### 4. Structural and novel second pass should be isolated from functional scoring

Structural and novel cases matter, but they should not dominate the first iteration of
functional-typing improvement. The initial objective is to improve routine short-read
functional typing without regressing on hard cases.

Expected intervention:

* first improve core functional typing on ordinary cases
* then add structural/novel refinement as a second pass

## Execution order

1. Add benchmark outputs for:
   * `3_digit_f1`
   * `5_digit_f1`
   * `7_digit_f1`
   * `copy_number_f1`
2. Implement mapping-stage multi-map policy execution:
   * `discard`
   * `best-only`
   * `weighted`
3. Add functional-region-aware typing scoring in `graphkir2`
4. Add confidence-aware downgrade behavior
5. Re-run difficult-case and HPRC sanity comparisons

## Decision rules

Prefer a change if it satisfies all of the following:

* improves `3-digit` or `5-digit`
* does not materially break CN
* does not create a severe runtime regression

Strongly prefer a change if it:

* improves both `3-digit` and `5-digit`
* keeps `7-digit` neutral or better
* reduces per-gene inconsistency on known difficult loci

## Benchmark outputs needed next

The benchmark harness should soon emit:

* per-gene `3-digit` F1
* per-gene `5-digit` F1
* overall `3-digit` / `5-digit` comparison against `Geny`
* multi-map-policy-tagged result summaries
* notes describing whether a change helps ordinary cases, difficult cases, or both
