# Repository Layout

Graph-KIR now follows a clearer top-level split inspired by a benchmark-oriented project layout:

* `graphkir/`: main Graph-KIR package
* `kir/`: wrappers for Graph-KIR and third-party KIR pipelines
* `examples/`: small runnable example inputs and expected outputs
* `data/`
  * `cohorts/`: cohort/sample manifests
  * `reference/`: static reference-side helper inputs
  * `groundtruth/`: evaluation truth tables
* `research/`: paper and benchmark scripts
* `docs/`: documentation and manuscript material

Design intent:

* Keep installable code separate from paper workflow code.
* Keep static datasets out of the repository root.
* Make it obvious where to add future benchmark assets without polluting `research/`.
