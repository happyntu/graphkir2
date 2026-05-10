"""Benchmark helpers for graphkir2."""

from .execute import LegacyExecutionResult, execute_legacy_command
from .eval import AlleleEvaluationSummary, evaluate_allele_calls
from .legacy import LegacyCommandPlan, build_legacy_graphkir_command
from .preset import BenchmarkPreset, load_benchmark_preset
from .results import (
    BenchmarkResultArtifact,
    apply_allele_evaluation,
    apply_execution_metrics,
    build_result_artifact,
    write_result_bundle,
)
from .runner import BenchmarkLabel, BenchmarkSummary, ComparisonSummary

__all__ = [
    "BenchmarkLabel",
    "BenchmarkSummary",
    "ComparisonSummary",
    "LegacyCommandPlan",
    "LegacyExecutionResult",
    "BenchmarkPreset",
    "BenchmarkResultArtifact",
    "AlleleEvaluationSummary",
    "apply_allele_evaluation",
    "apply_execution_metrics",
    "build_legacy_graphkir_command",
    "build_result_artifact",
    "execute_legacy_command",
    "evaluate_allele_calls",
    "load_benchmark_preset",
    "write_result_bundle",
]
