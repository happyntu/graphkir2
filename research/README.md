## Research Layout

This directory contains paper-specific and benchmark-specific scripts.

Boundaries:

* `research/*.py`: orchestration, evaluation, and exploratory scripts
* `research/*.template`: SLURM templates used by the paper workflow
* `data/`: static inputs consumed by research scripts, including cohort manifests,
  reference-side helper files, and HPRC ground truth

This separation keeps the installable packages (`graphkir/`, `kir/`) independent from
paper artifacts while preserving the existing research workflow.
