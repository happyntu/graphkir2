# Benchmarks

This directory is reserved for old-vs-new comparison work.

Suggested structure:

* `benchmarks/configs/` — benchmark manifests and run settings
* `benchmarks/results/` — committed summary tables and plots
* `benchmarks/scripts/` — orchestration helpers, if they grow beyond `research/`

Benchmark goals:

* compare legacy `graphkir` against refactored `graphkir2`
* track runtime, memory, copy-number accuracy, and allele-typing accuracy
* protect against regressions during refactor work

Current scaffold:

* `benchmarks/configs/examples_skip_extraction.json` — sample preset for the bundled `examples/cohort.csv`
* `benchmarks/scripts/prepare_synthetic_fixture.py` — generates a small synthetic cohort, truth TSV, manifest CSV, and preset JSON
* `benchmarks/scripts/run_compare.py` — prints, executes, or evaluates an old-vs-new benchmark artifact from a preset
* `benchmarks/scripts/rerun_typing.py` — re-runs only the legacy typing stage from an existing preset for multi-map ablation work
* `benchmarks/results/` — target folder for canonical compare artifacts

Artifact schema:

* `metadata` — schema version, label, preset path, manifest path, timestamp, git commit, status
* `comparison` — `graphkir2` planning summary plus the equivalent legacy `graphkir` command
* `metrics` — runtime and accuracy fields collected from execution and evaluation
  * `three_digit_f1`
  * `five_digit_f1`
  * `seven_digit_f1`
  * `copy_number_f1`
  * `allele_f1`
* `collector` — notes about collection mode and pending metrics

Examples:

```bash
python benchmarks/scripts/run_compare.py --config benchmarks/configs/examples_skip_extraction.json
python benchmarks/scripts/run_compare.py --config benchmarks/configs/examples_skip_extraction.json --json
python benchmarks/scripts/run_compare.py --config benchmarks/configs/examples_skip_extraction.json --output benchmarks/results/examples-skip-extraction.compare.json
python benchmarks/scripts/run_compare.py --config benchmarks/configs/synthetic-functional4.json --execute-legacy --evaluate --output-dir benchmarks/results/synthetic-functional4-warm-bundle
python benchmarks/scripts/prepare_synthetic_fixture.py --label synthetic-smoke3 --samples 2 --depth 1
python benchmarks/scripts/run_compare.py --config benchmarks/configs/smoke_examples.json --output-dir benchmarks/results/smoke-examples
python benchmarks/scripts/run_compare.py --config benchmarks/configs/hprc_real_sanity.json --evaluate --output-dir benchmarks/results/hprc-real-sanity
python benchmarks/scripts/rerun_typing.py --config benchmarks/configs/synthetic-functional8.json --multi-map-mode keep-all --output-tsv benchmarks/results/synthetic-functional8/keepmulti.allele.tsv
python benchmarks/scripts/rerun_typing.py --config benchmarks/configs/synthetic-functional8.json --multi-map-mode best-only --output-tsv benchmarks/results/synthetic-functional8/bestonly.allele.tsv
python benchmarks/scripts/rerun_typing.py --config benchmarks/configs/synthetic-functional8.json --multi-map-mode weighted --output-tsv benchmarks/results/synthetic-functional8/weighted.allele.tsv
python benchmarks/scripts/rerun_typing.py --config benchmarks/configs/synthetic-functional8.json --multi-map-mode best-only --exon-weight 2.0 --output-tsv benchmarks/results/synthetic-functional8/bestonly_exon2.allele.tsv
graphkir2 --input-csv benchmarks/generated/synthetic-functional8/manifest.csv --index-folder benchmarks/generated/synthetic-functional8/index_ab --msa-type ab --msa-no-exon-only-allele --multi-map-policy best-only --allele-exon-weight 2.0 --output-cohort-name benchmarks/results/synthetic-functional8/cohort --print-compare-summary
graphkir2 --input-csv examples/cohort.csv --step-skip-extraction --output-cohort-name example_data/cohort --print-compare-summary
```
