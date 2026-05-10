"""Private-variant support helpers for functional typing reranking."""

from __future__ import annotations

from collections import defaultdict
from typing import Any


VariantSupport = defaultdict[str, float]
GeneGroups = tuple[frozenset[str], ...]


def read_pair_name(read: object) -> str:
    """Extract a physical read-pair name from a PairRead-like object."""
    if read.l_sam:  # type: ignore[attr-defined]
        return read.l_sam.split("\t", 1)[0]  # type: ignore[attr-defined]
    if read.r_sam:  # type: ignore[attr-defined]
        return read.r_sam.split("\t", 1)[0]  # type: ignore[attr-defined]
    return ""


def pure_gene(backbone: str) -> str:
    """Extract the pure gene from a backbone or allele name."""
    return backbone.split("*", 1)[0]


def limit_allele_resolution(allele: str, resolution: int) -> str:
    """Limit an allele name to a fixed numeric field length."""
    gene, field = allele.split("*", 1)
    return f"{gene}*{field[:resolution]}"


def parse_gene_groups(spec: str) -> GeneGroups:
    """Parse comma-separated slash-delimited gene groups."""
    groups: list[frozenset[str]] = []
    for field in spec.split(","):
        group = frozenset(gene.strip() for gene in field.split("/") if gene.strip())
        if len(group) >= 2:
            groups.append(group)
    return tuple(groups)


def parse_gene_set(spec: str) -> frozenset[str]:
    """Parse a comma-separated pure-gene set."""
    return frozenset(gene.strip() for gene in spec.split(",") if gene.strip())


def neutralize_cross_gene_reads(
    reads: list[Any],
    gene_groups: GeneGroups = (),
    target_genes: frozenset[str] = frozenset(),
) -> None:
    """Remove allele-specific target evidence from cross-gene ambiguous reads."""
    reads_by_name: dict[str, list[Any]] = defaultdict(list)
    for read in reads:
        reads_by_name[read_pair_name(read)].append(read)
    for read_group in reads_by_name.values():
        genes = {pure_gene(read.backbone) for read in read_group}
        if len(genes) <= 1:
            continue
        if gene_groups and not any(len(genes & group) >= 2 for group in gene_groups):
            continue
        for read in read_group:
            if target_genes and pure_gene(read.backbone) not in target_genes:
                continue
            read.weight = 0.0
            read.ambiguous_weight = 1.0


def collect_variant_support(reads: list[Any]) -> tuple[VariantSupport, VariantSupport]:
    """Collect weighted positive and negative support per variant ID."""
    positive: VariantSupport = defaultdict(float)
    negative: VariantSupport = defaultdict(float)
    for read in reads:
        weight = float(read.weight)
        for variant_id in read.lpv + read.rpv:
            positive[variant_id] += weight
        for variant_id in read.lnv + read.rnv:
            negative[variant_id] += weight
    return positive, negative


def private_support_score(
    allele_names: list[str],
    allele_variants: dict[str, set[str]],
    positive: VariantSupport,
    negative: VariantSupport,
) -> float:
    """Score whether each allele has supported private variants."""
    score = 0.0
    if len(set(allele_names)) < len(allele_names):
        score -= 50.0
    for allele in allele_names:
        other_variants = set().union(
            *(allele_variants[other] for other in allele_names if other != allele)
        )
        private_variants = allele_variants[allele] - other_variants
        supported = sum(
            1
            for variant_id in private_variants
            if positive[variant_id] >= 2.0 and positive[variant_id] >= negative[variant_id]
        )
        unsupported = sum(
            1
            for variant_id in private_variants
            if negative[variant_id] >= 5.0 and positive[variant_id] == 0.0
        )
        net_support = sum(
            positive[variant_id] - negative[variant_id]
            for variant_id in private_variants
        )
        score += 5.0 * supported - 10.0 * unsupported + 0.1 * net_support
    return score


def select_with_private_support(
    result: object,
    allele_variants: dict[str, set[str]],
    positive: VariantSupport,
    negative: VariantSupport,
    support_lambda: float,
    support_window: float,
    min_fraction_ratio: float,
) -> list[str]:
    """Select a candidate using likelihood plus private-variant support."""
    if not support_lambda or result.isFail():  # type: ignore[attr-defined]
        return result.selectBest(min_fraction_ratio=min_fraction_ratio)  # type: ignore[attr-defined]

    max_value = float(result.value[0])  # type: ignore[attr-defined]
    best_score = float("-inf")
    best_index = 0
    for index, value in enumerate(result.value):  # type: ignore[attr-defined]
        if max_value - float(value) > support_window:
            continue
        allele_names = result.allele_name[index]  # type: ignore[attr-defined]
        support = private_support_score(
            allele_names,
            allele_variants,
            positive,
            negative,
        )
        score = float(value) + support_lambda * support
        if score > best_score:
            best_score = score
            best_index = index
    return result.allele_name[best_index]  # type: ignore[attr-defined]


def select_with_highest_suffix_tie_break(
    result: object,
    selected: list[str],
    min_fraction_ratio: float,
    resolution: int = 5,
    value_epsilon: float = 1e-9,
) -> list[str]:
    """Prefer the highest full allele suffix only for exact same-resolution ties."""
    if result.isFail():  # type: ignore[attr-defined]
        return selected

    expect_prob = 1 / int(result.n)  # type: ignore[attr-defined]
    selected_key = sorted(limit_allele_resolution(allele, resolution) for allele in selected)
    selected_index = 0
    for index, allele_names in enumerate(result.allele_name):  # type: ignore[attr-defined]
        if allele_names == selected:
            selected_index = index
            break
    selected_value = float(result.value[selected_index])  # type: ignore[attr-defined]
    candidates: list[list[str]] = []
    for index, allele_names in enumerate(result.allele_name):  # type: ignore[attr-defined]
        if abs(float(result.value[index]) - selected_value) > value_epsilon:  # type: ignore[attr-defined]
            continue
        if not all(
            fraction >= expect_prob * min_fraction_ratio
            for fraction in result.fraction[index]  # type: ignore[attr-defined]
        ):
            continue
        candidate_key = sorted(
            limit_allele_resolution(allele, resolution)
            for allele in allele_names
        )
        if candidate_key != selected_key:
            continue
        candidates.append(allele_names)
    if not candidates:
        return selected
    return max(candidates, key=lambda alleles: sorted(alleles))
