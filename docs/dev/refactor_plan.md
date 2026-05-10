# Refactor Plan

## Goal

Refactor Graph-KIR without losing the current implementation as a correctness and
performance baseline.

## Repository strategy

* Keep legacy implementation in:
  * `graphkir/`
  * `kir/`
* Build the new implementation in:
  * `src/graphkir2/`
* Compare both in:
  * `benchmarks/`

## Why this split

* avoids import and CLI conflicts during the transition
* preserves a stable baseline for performance comparisons
* makes it possible to incrementally reimplement subsystems

## Recommended migration order

1. define shared domain models and configuration objects in `src/graphkir2/`
2. extract stage boundaries for mapping, CN, and typing
3. add focused unit tests for new code
4. benchmark new stages against legacy outputs
5. only replace the legacy CLI after parity is demonstrated
