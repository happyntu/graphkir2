"""Benchmark evaluation helpers for graphkir2."""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from graphkir.utils import getGeneName, limitAlleleField


Resolution = int
AlleleCounter = Counter[str]
GeneAlleleCounters = dict[str, AlleleCounter]


@dataclass(frozen=True)
class ResolutionMetrics:
    """Precision/recall/F1 counts for one resolution."""

    truth_total: int
    pred_total: int
    true_positive: int
    precision: float
    recall: float
    f1: float


@dataclass(frozen=True)
class GeneMetrics:
    """Per-gene multi-resolution metrics."""

    gene: str
    three_digit: ResolutionMetrics
    five_digit: ResolutionMetrics
    seven_digit: ResolutionMetrics


@dataclass(frozen=True)
class AlleleEvaluationSummary:
    """Cohort-level allele evaluation summary."""

    truth_path: str
    pred_path: str
    truth_samples: int
    pred_samples: int
    matched_samples: int
    missing_truth_samples: tuple[str, ...]
    missing_pred_samples: tuple[str, ...]
    three_digit: ResolutionMetrics
    five_digit: ResolutionMetrics
    seven_digit: ResolutionMetrics
    per_gene: tuple[GeneMetrics, ...]


def _calc_metrics(truth_total: int, pred_total: int, true_positive: int) -> ResolutionMetrics:
    """Build precision/recall/F1 metrics from raw counts."""
    precision = true_positive / pred_total if pred_total else 0.0
    recall = true_positive / truth_total if truth_total else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return ResolutionMetrics(
        truth_total=truth_total,
        pred_total=pred_total,
        true_positive=true_positive,
        precision=precision,
        recall=recall,
        f1=f1,
    )


def _normalize_allele(raw: str) -> str:
    """Normalize an allele token for benchmark comparison."""
    token = raw.strip()
    if not token or "*" not in token:
        return ""
    token = re.sub(r"[\$\#\+\=\-]+$", "", token)
    gene, field = token.split("*", 1)
    field = field.split("e", 1)[0]
    field = re.sub(r"[^0-9A-Za-z]+$", "", field)
    if not gene or not field:
        return ""
    return f"{gene}*{field}"


def _read_allele_rows(path: str) -> dict[str, list[str]]:
    """Read a merged allele-style TSV into sample -> allele list."""
    rows: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Allele file is missing a header: {path}")
        for row in reader:
            sample_id = (
                row.get("id")
                or row.get("sample_id")
                or row.get("name")
                or ""
            ).strip()
            if not sample_id:
                continue
            if "name" in row and not row.get("id") and not row.get("sample_id"):
                raw_name = Path(sample_id).name
                match = re.search(r"\.(\d{2})\.", raw_name)
                sample_id = match.group(1) if match else raw_name
            allele_field = row.get("alleles", "").strip()
            if not allele_field:
                rows[sample_id] = []
                continue
            alleles = [_normalize_allele(token) for token in allele_field.split("_")]
            rows[sample_id] = [allele for allele in alleles if allele]
    return rows


def _build_gene_counters(
    sample_alleles: dict[str, list[str]],
    resolution: Resolution,
) -> dict[str, GeneAlleleCounters]:
    """Convert sample allele strings into per-gene resolution counters."""
    cohort: dict[str, GeneAlleleCounters] = {}
    for sample_id, alleles in sample_alleles.items():
        gene_counters: dict[str, AlleleCounter] = defaultdict(Counter)
        for allele in alleles:
            gene = getGeneName(allele)
            gene_counters[gene][limitAlleleField(allele, resolution)] += 1
        cohort[sample_id] = dict(gene_counters)
    return cohort


def _evaluate_resolution(
    truth: dict[str, GeneAlleleCounters],
    pred: dict[str, GeneAlleleCounters],
    gene_filter: str | None = None,
) -> ResolutionMetrics:
    """Evaluate one resolution over the cohort, optionally for one gene only."""
    truth_total = 0
    pred_total = 0
    true_positive = 0
    sample_ids = set(truth) | set(pred)
    for sample_id in sample_ids:
        truth_genes = truth.get(sample_id, {})
        pred_genes = pred.get(sample_id, {})
        genes = set(truth_genes) | set(pred_genes)
        if gene_filter is not None:
            genes = {gene for gene in genes if gene == gene_filter}
        for gene in genes:
            truth_counter = truth_genes.get(gene, Counter())
            pred_counter = pred_genes.get(gene, Counter())
            truth_total += sum(truth_counter.values())
            pred_total += sum(pred_counter.values())
            true_positive += sum((truth_counter & pred_counter).values())
    return _calc_metrics(truth_total, pred_total, true_positive)


def evaluate_allele_calls(truth_path: str, pred_path: str) -> AlleleEvaluationSummary:
    """Evaluate merged allele predictions at 3/5/7-digit resolution."""
    truth_rows = _read_allele_rows(truth_path)
    pred_rows = _read_allele_rows(pred_path)

    truth3 = _build_gene_counters(truth_rows, 3)
    pred3 = _build_gene_counters(pred_rows, 3)
    truth5 = _build_gene_counters(truth_rows, 5)
    pred5 = _build_gene_counters(pred_rows, 5)
    truth7 = _build_gene_counters(truth_rows, 7)
    pred7 = _build_gene_counters(pred_rows, 7)

    all_genes = sorted(
        {
            gene
            for cohort in (truth7, pred7)
            for sample in cohort.values()
            for gene in sample
        }
    )
    per_gene = tuple(
        GeneMetrics(
            gene=gene,
            three_digit=_evaluate_resolution(truth3, pred3, gene_filter=gene),
            five_digit=_evaluate_resolution(truth5, pred5, gene_filter=gene),
            seven_digit=_evaluate_resolution(truth7, pred7, gene_filter=gene),
        )
        for gene in all_genes
    )
    truth_ids = set(truth_rows)
    pred_ids = set(pred_rows)
    return AlleleEvaluationSummary(
        truth_path=truth_path,
        pred_path=pred_path,
        truth_samples=len(truth_rows),
        pred_samples=len(pred_rows),
        matched_samples=len(truth_ids & pred_ids),
        missing_truth_samples=tuple(sorted(pred_ids - truth_ids)),
        missing_pred_samples=tuple(sorted(truth_ids - pred_ids)),
        three_digit=_evaluate_resolution(truth3, pred3),
        five_digit=_evaluate_resolution(truth5, pred5),
        seven_digit=_evaluate_resolution(truth7, pred7),
        per_gene=per_gene,
    )
