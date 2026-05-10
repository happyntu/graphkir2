# Repository Layout

Graph-KIR now follows a clearer top-level split inspired by a benchmark-oriented project layout:

* `graphkir/`: main Graph-KIR package
* `kir/`: wrappers for Graph-KIR and third-party KIR pipelines
* `src/graphkir2/`: next-generation implementation under refactor
* `examples/`: small runnable example inputs and expected outputs
* `data/`
  * `cohorts/`: cohort/sample manifests
  * `reference/`: static reference-side helper inputs
  * `groundtruth/`: evaluation truth tables
* `research/`: paper and benchmark scripts
* `benchmarks/`: compare presets, results, and benchmark-oriented scripts
* `docs/`: architecture, roadmap, developer notes, ADRs, research notes, and manuscript material

Design intent:

* Keep installable code separate from paper workflow code.
* Keep static datasets out of the repository root.
* Make it obvious where to add future benchmark assets without polluting `research/`.
* Keep accepted design decisions separate from exploratory benchmark notes.
