# ADRs

This directory stores accepted architectural decisions for the `graphkir2` refactor.

Use an ADR when:

* a design choice materially changes how `graphkir2` is structured
* a benchmark interpretation rule should become a durable project policy
* a method decision needs to remain stable across refactor iterations

Planned first ADRs:

* `ADR-001`: keep legacy `graphkir` as the benchmark baseline during refactor
* `ADR-002`: separate benchmark planning artifacts from executed benchmark results
* `ADR-003`: prioritize copy number and `7-digit` performance over a single aggregate score

Do not use ADRs for exploratory notes or one-off benchmark observations. Put those in
`docs/research/`.
