# GraphKir2 Synthetic-First Benchmark Strategy

## Position

`graphkir2` should tune performance and method behavior on synthetic data first, then use
real data for sanity and generalization checks.

This is the right order for the current stage because:

* the refactor is still changing stage boundaries
* runtime and artifact sizes are not yet stable
* `3-digit/5-digit` functional typing is the primary optimization target
* repeated full real-data runs are too expensive for broad parameter search

## Why synthetic data first

Synthetic cohorts give three things that real-data-first tuning does not:

1. cheap iteration on runtime and memory behavior
2. direct truth for `3-digit`, `5-digit`, and `7-digit` comparison
3. controlled difficult cases for multi-map and structural ablation

This means synthetic data should be used for:

* profiling stage-level runtime
* testing multi-map policy changes
* testing likelihood-weighting changes
* testing confidence downgrade logic

Real data should be used later for:

* confirming gains are not synthetic-only artifacts
* checking structural and novel behavior on realistic samples
* validating practical call stability

## Recommended execution ladder

### Level 0: planning smoke

Purpose:

* verify manifests
* verify artifact naming
* verify compare bundle output

Current preset family:

* `smoke_examples`

### Level 1: synthetic runtime tuning

Purpose:

* find the expensive stages first
* reduce execution time before broader algorithm search

Questions:

* is mapping still dominating runtime?
* does multi-map handling increase artifact size too much?
* does `best-only` capture most of the benefit of `weighted` at much lower cost?

### Level 2: synthetic functional-typing tuning

Purpose:

* optimize `3-digit/5-digit` functional typing against truth

Questions:

* which policy improves `3-digit` most?
* which policy improves `5-digit` most?
* does exon-weighted or CDS-aware scoring help without damaging CN?

### Level 3: synthetic difficult-case ablation

Purpose:

* isolate multi-map, structural, and novel failure modes

Questions:

* are gains limited to ordinary cases?
* do difficult loci regress?
* which interventions help routine typing vs hard cases?

### Level 4: real-data sanity

Purpose:

* confirm that synthetic gains transfer to HPRC-like data

Questions:

* does the synthetic winner still help `3-digit/5-digit` on real samples?
* is runtime still acceptable?
* do warning rates or unresolved calls become impractical?

## Decision rule

Do not promote a method based on real data alone if synthetic profiling still shows:

* large stage-level inefficiencies
* unstable runtime scaling
* unresolved `3-digit/5-digit` regressions

Conversely, do not promote a synthetic winner to default status until:

* it has been checked on real data
* CN remains acceptable
* runtime remains practical

## Immediate next benchmark deliverables

1. add runtime fields to the result bundle at stage granularity
2. add synthetic functional-typing presets with explicit truth files
3. add synthetic difficult-case presets for multi-map and structural ablation
4. run real-data sanity only for the best synthetic candidates
