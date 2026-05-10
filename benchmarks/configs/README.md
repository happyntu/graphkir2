# Benchmark Configs

This directory stores benchmark presets for `graphkir2`.

Current state:

* presets are JSON because the current compare runner is planning-only
* each preset maps to the `BenchmarkPreset` schema used by `src/graphkir2/benchmark/`
* labels should stay stable because result artifact names depend on them
* several non-smoke preset families currently reuse `examples/cohort.csv` as a runnable
  placeholder until dedicated Graph-KIR sample manifests are prepared

Suggested preset families:

* smoke
* simulated runtime tuning
* simulated CN
* simulated typing
* simulated difficult cases
* real-data sanity
* multi-map rescue ablation
* structural / novel ablation
* thread scaling
* index reuse

Recommended usage order:

1. smoke
2. simulated runtime tuning
3. simulated typing / CN tuning
4. simulated difficult cases
5. real-data sanity

Current note:

* the checked-in `simulated_*` presets are still scaffold presets until dedicated
  generated synthetic manifests are wired into the repo
* `benchmarks/scripts/prepare_synthetic_fixture.py` can generate a small runnable
  synthetic manifest, truth TSV, and matching preset JSON
* use `--gene-whitelist` to create a smaller functional-typing panel when the goal is
  early `3-digit/5-digit` tuning rather than full-locus stress testing
