"""Structured benchmark result artifacts."""

from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .execute import LegacyExecutionResult
from .eval import AlleleEvaluationSummary
from .runner import ComparisonSummary


@dataclass(frozen=True)
class BenchmarkArtifactMetadata:
    """Metadata stored alongside a benchmark comparison artifact."""

    schema_version: str
    benchmark_label: str
    preset_path: str
    manifest_path: str
    generated_at_utc: str
    git_commit: str
    status: str


@dataclass(frozen=True)
class BenchmarkMetricPlaceholders:
    """Metric slots reserved for future executed benchmark runs."""

    runtime_seconds: float | None = None
    max_rss_mb: float | None = None
    three_digit_f1: float | None = None
    five_digit_f1: float | None = None
    seven_digit_f1: float | None = None
    copy_number_f1: float | None = None
    allele_f1: float | None = None


@dataclass(frozen=True)
class BenchmarkCollectorNotes:
    """Notes describing what is and is not collected in the current scaffold."""

    mode: str
    pending_metrics: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class BenchmarkResultArtifact:
    """Canonical JSON artifact for old-vs-new benchmark comparisons."""

    metadata: BenchmarkArtifactMetadata
    comparison: ComparisonSummary
    metrics: BenchmarkMetricPlaceholders
    collector: BenchmarkCollectorNotes

    def to_json(self) -> str:
        """Serialize the result artifact as formatted JSON."""
        return json.dumps(asdict(self), indent=2)

    def write_json(self, path: str) -> None:
        """Write the result artifact to disk."""
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.to_json() + "\n", encoding="utf-8")


def _pending_metric_names(metrics: BenchmarkMetricPlaceholders) -> list[str]:
    """List metrics that are still unset in the artifact."""
    pending: list[str] = []
    if metrics.runtime_seconds is None:
        pending.append("runtime_seconds")
    if metrics.max_rss_mb is None:
        pending.append("max_rss_mb")
    if metrics.three_digit_f1 is None:
        pending.append("three_digit_f1")
    if metrics.five_digit_f1 is None:
        pending.append("five_digit_f1")
    if metrics.seven_digit_f1 is None:
        pending.append("seven_digit_f1")
    if metrics.copy_number_f1 is None:
        pending.append("copy_number_f1")
    if metrics.allele_f1 is None:
        pending.append("allele_f1")
    return pending


def _metric_text(value: float | None) -> str:
    """Render an optional numeric metric for TSV output."""
    if value is None:
        return ""
    return str(value)


def build_summary_tsv(artifact: BenchmarkResultArtifact) -> str:
    """Build a one-row benchmark summary TSV."""
    header = [
        "benchmark_label",
        "status",
        "manifest_path",
        "multi_map_policy",
        "typing_strategy",
        "legacy_exact_match",
        "runtime_seconds",
        "max_rss_mb",
        "three_digit_f1",
        "five_digit_f1",
        "seven_digit_f1",
        "copy_number_f1",
        "allele_f1",
    ]
    row = [
        artifact.metadata.benchmark_label,
        artifact.metadata.status,
        artifact.metadata.manifest_path,
        artifact.comparison.benchmark.mapping.multi_map_policy,
        artifact.comparison.benchmark.typing.strategy,
        str(artifact.comparison.legacy_command.exact_config_match),
        _metric_text(artifact.metrics.runtime_seconds),
        _metric_text(artifact.metrics.max_rss_mb),
        _metric_text(artifact.metrics.three_digit_f1),
        _metric_text(artifact.metrics.five_digit_f1),
        _metric_text(artifact.metrics.seven_digit_f1),
        _metric_text(artifact.metrics.copy_number_f1),
        _metric_text(artifact.metrics.allele_f1),
    ]
    return "\t".join(header) + "\n" + "\t".join(row) + "\n"


def build_per_gene_tsv(
    evaluation: AlleleEvaluationSummary | None = None,
) -> str:
    """Build a per-gene TSV, optionally populated from allele evaluation."""
    header = [
        "gene",
        "truth_total_3d",
        "pred_total_3d",
        "tp_3d",
        "three_digit_f1",
        "truth_total_5d",
        "pred_total_5d",
        "tp_5d",
        "five_digit_f1",
        "truth_total_7d",
        "pred_total_7d",
        "tp_7d",
        "seven_digit_f1",
    ]
    lines = ["\t".join(header)]
    if evaluation is None:
        return "\n".join(lines) + "\n"
    for gene in evaluation.per_gene:
        lines.append(
            "\t".join(
                [
                    gene.gene,
                    str(gene.three_digit.truth_total),
                    str(gene.three_digit.pred_total),
                    str(gene.three_digit.true_positive),
                    str(gene.three_digit.f1),
                    str(gene.five_digit.truth_total),
                    str(gene.five_digit.pred_total),
                    str(gene.five_digit.true_positive),
                    str(gene.five_digit.f1),
                    str(gene.seven_digit.truth_total),
                    str(gene.seven_digit.pred_total),
                    str(gene.seven_digit.true_positive),
                    str(gene.seven_digit.f1),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def build_failure_modes_tsv(artifact: BenchmarkResultArtifact) -> str:
    """Build a failure-mode TSV skeleton with expected categories."""
    header = ["benchmark_label", "category", "metric", "value", "notes"]
    rows = [
        [
            artifact.metadata.benchmark_label,
            "multi_map",
            "policy",
            artifact.comparison.benchmark.mapping.multi_map_policy,
            "planning placeholder",
        ],
        [
            artifact.metadata.benchmark_label,
            "structural",
            "warning_rate",
            "",
            "planning placeholder",
        ],
        [
            artifact.metadata.benchmark_label,
            "novel",
            "novel_like_rate",
            "",
            "planning placeholder",
        ],
    ]
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    return "\n".join(lines) + "\n"


def build_notes_markdown(artifact: BenchmarkResultArtifact) -> str:
    """Build a human-readable notes file for the benchmark result bundle."""
    lines = [
        f"# {artifact.metadata.benchmark_label}",
        "",
        "## Status",
        "",
        f"* status: `{artifact.metadata.status}`",
        f"* mode: `{artifact.collector.mode}`",
        f"* manifest: `{artifact.metadata.manifest_path}`",
        f"* preset: `{artifact.metadata.preset_path}`",
        f"* multi-map policy: `{artifact.comparison.benchmark.mapping.multi_map_policy}`",
        f"* typing strategy: `{artifact.comparison.benchmark.typing.strategy}`",
        f"* legacy exact match: `{artifact.comparison.legacy_command.exact_config_match}`",
        f"* runtime seconds: `{_metric_text(artifact.metrics.runtime_seconds)}`",
        f"* max RSS MB: `{_metric_text(artifact.metrics.max_rss_mb)}`",
        "",
        "## Pending Metrics",
        "",
    ]
    lines.extend(f"* `{metric}`" for metric in artifact.collector.pending_metrics)
    if artifact.comparison.legacy_command.notes:
        lines.extend(
            [
                "",
                "## Legacy Notes",
                "",
            ]
        )
        lines.extend(f"* {note}" for note in artifact.comparison.legacy_command.notes)
    return "\n".join(lines) + "\n"


def write_result_bundle(
    output_dir: str,
    artifact: BenchmarkResultArtifact,
    evaluation: AlleleEvaluationSummary | None = None,
    execution: LegacyExecutionResult | None = None,
) -> None:
    """Write the canonical result bundle for a benchmark label."""
    result_dir = Path(output_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    artifact.write_json(str(result_dir / "compare.json"))
    (result_dir / "summary.tsv").write_text(build_summary_tsv(artifact), encoding="utf-8")
    (result_dir / "per_gene.tsv").write_text(
        build_per_gene_tsv(evaluation),
        encoding="utf-8",
    )
    (result_dir / "failure_modes.tsv").write_text(
        build_failure_modes_tsv(artifact),
        encoding="utf-8",
    )
    (result_dir / "notes.md").write_text(build_notes_markdown(artifact), encoding="utf-8")
    if execution is not None:
        (result_dir / "legacy.stdout.log").write_text(execution.stdout, encoding="utf-8")
        (result_dir / "legacy.stderr.log").write_text(execution.stderr, encoding="utf-8")


def apply_allele_evaluation(
    artifact: BenchmarkResultArtifact,
    evaluation: AlleleEvaluationSummary,
) -> BenchmarkResultArtifact:
    """Return a copy of the artifact with allele metrics filled in."""
    metrics = BenchmarkMetricPlaceholders(
        runtime_seconds=artifact.metrics.runtime_seconds,
        max_rss_mb=artifact.metrics.max_rss_mb,
        three_digit_f1=evaluation.three_digit.f1,
        five_digit_f1=evaluation.five_digit.f1,
        seven_digit_f1=evaluation.seven_digit.f1,
        copy_number_f1=artifact.metrics.copy_number_f1,
        allele_f1=evaluation.seven_digit.f1,
    )
    return BenchmarkResultArtifact(
        metadata=BenchmarkArtifactMetadata(
            schema_version=artifact.metadata.schema_version,
            benchmark_label=artifact.metadata.benchmark_label,
            preset_path=artifact.metadata.preset_path,
            manifest_path=artifact.metadata.manifest_path,
            generated_at_utc=artifact.metadata.generated_at_utc,
            git_commit=artifact.metadata.git_commit,
            status="evaluated",
        ),
        comparison=artifact.comparison,
        metrics=metrics,
        collector=BenchmarkCollectorNotes(
            mode=(
                "legacy-executed+allele-evaluated"
                if artifact.metrics.runtime_seconds is not None
                else "allele-evaluated"
            ),
            pending_metrics=_pending_metric_names(metrics),
        ),
    )


def apply_execution_metrics(
    artifact: BenchmarkResultArtifact,
    execution: LegacyExecutionResult,
) -> BenchmarkResultArtifact:
    """Return a copy of the artifact with runtime and memory metrics filled in."""
    metrics = BenchmarkMetricPlaceholders(
        runtime_seconds=execution.runtime_seconds,
        max_rss_mb=execution.max_rss_mb,
        three_digit_f1=artifact.metrics.three_digit_f1,
        five_digit_f1=artifact.metrics.five_digit_f1,
        seven_digit_f1=artifact.metrics.seven_digit_f1,
        copy_number_f1=artifact.metrics.copy_number_f1,
        allele_f1=artifact.metrics.allele_f1,
    )
    return BenchmarkResultArtifact(
        metadata=BenchmarkArtifactMetadata(
            schema_version=artifact.metadata.schema_version,
            benchmark_label=artifact.metadata.benchmark_label,
            preset_path=artifact.metadata.preset_path,
            manifest_path=artifact.metadata.manifest_path,
            generated_at_utc=artifact.metadata.generated_at_utc,
            git_commit=artifact.metadata.git_commit,
            status="executed",
        ),
        comparison=artifact.comparison,
        metrics=metrics,
        collector=BenchmarkCollectorNotes(
            mode="legacy-executed",
            pending_metrics=_pending_metric_names(metrics),
        ),
    )


def detect_git_commit(repo_root: Path) -> str:
    """Resolve the current git commit when available."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def build_result_artifact(
    comparison: ComparisonSummary,
    preset_path: str,
    repo_root: Path,
) -> BenchmarkResultArtifact:
    """Build a planning-only benchmark artifact from the comparison summary."""
    metadata = BenchmarkArtifactMetadata(
        schema_version="1",
        benchmark_label=comparison.benchmark.label.name,
        preset_path=preset_path,
        manifest_path=comparison.benchmark.manifest_path,
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        git_commit=detect_git_commit(repo_root),
        status="planned",
    )
    collector = BenchmarkCollectorNotes(
        mode="planning-only",
        pending_metrics=_pending_metric_names(BenchmarkMetricPlaceholders()),
    )
    return BenchmarkResultArtifact(
        metadata=metadata,
        comparison=comparison,
        metrics=BenchmarkMetricPlaceholders(),
        collector=collector,
    )
