"""Private-variant support helpers for functional typing reranking."""

from __future__ import annotations

from collections import defaultdict
from typing import Any


VariantSupport = defaultdict[str, float]
GeneGroups = tuple[frozenset[str], ...]
NameSet = frozenset[str]


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


def parse_name_set(spec: str) -> NameSet:
    """Parse a comma-separated set of allele or gene names."""
    return frozenset(name.strip() for name in spec.split(",") if name.strip())


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


def cross_gene_read_names(reads: list[Any], gene_groups: GeneGroups = ()) -> set[str]:
    """Return physical read names that span configured cross-gene groups."""
    reads_by_name: dict[str, list[Any]] = defaultdict(list)
    for read in reads:
        reads_by_name[read_pair_name(read)].append(read)

    cross_names: set[str] = set()
    for name, read_group in reads_by_name.items():
        genes = {pure_gene(read.backbone) for read in read_group}
        if len(genes) <= 1:
            continue
        if gene_groups and not any(len(genes & group) >= 2 for group in gene_groups):
            continue
        cross_names.add(name)
    return cross_names


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


def collect_private_variants(
    allele_names: list[str],
    allele_variants: dict[str, set[str]],
) -> set[str]:
    """Collect variants private to each selected allele within a genotype."""
    private_variants: set[str] = set()
    for allele in allele_names:
        other_variants = set().union(
            *(allele_variants[other] for other in allele_names if other != allele)
        )
        private_variants.update(allele_variants[allele] - other_variants)
    return private_variants


def private_positive_cross_gene_ratio(
    reads: list[Any],
    allele_names: list[str],
    allele_variants: dict[str, set[str]],
    gene_groups: GeneGroups = (),
) -> float:
    """Estimate how much selected-private positive support is cross-gene ambiguous."""
    private_variants = collect_private_variants(allele_names, allele_variants)
    if not private_variants:
        return 0.0

    cross_names = cross_gene_read_names(reads, gene_groups)
    cross_support = 0.0
    total_support = 0.0
    for read in reads:
        weight = float(read.weight)
        if weight <= 0.0:
            continue
        read_support = sum(
            weight
            for variant_id in read.lpv + read.rpv
            if variant_id in private_variants
        )
        total_support += read_support
        if read_pair_name(read) in cross_names:
            cross_support += read_support
    if total_support <= 0.0:
        return 0.0
    return cross_support / total_support


def selected_has_name_prefix(allele_names: list[str], prefixes: NameSet) -> bool:
    """Check whether a selected allele starts with one of the configured prefixes."""
    return any(
        allele.startswith(prefix)
        for allele in allele_names
        for prefix in prefixes
    )


def should_apply_conditional_private_support(
    selected: list[str],
    reads: list[Any],
    allele_variants: dict[str, set[str]],
    condition_alleles: NameSet,
    gene_groups: GeneGroups,
    min_cross_gene_ratio: float,
) -> bool:
    """Decide whether private-support rescue is justified for a selected call."""
    if condition_alleles and not selected_has_name_prefix(selected, condition_alleles):
        return False
    if min_cross_gene_ratio <= 0.0:
        return True
    ratio = private_positive_cross_gene_ratio(
        reads,
        selected,
        allele_variants,
        gene_groups,
    )
    return ratio >= min_cross_gene_ratio


def should_use_discard_fallback(
    selected: list[str],
    base_selected: list[str],
    residual_alleles: NameSet,
    introduced_alleles: NameSet,
    cross_gene_ratio: float,
    introduced_max_cross_gene_ratio: float,
    private_support: float,
    max_private_support: float,
    residual_min_cross_gene_ratio: float = 0.0,
) -> bool:
    """Decide whether a rescued call should fall back to discard-style evidence."""
    if private_support > max_private_support:
        return False
    if (
        residual_alleles
        and selected_has_name_prefix(selected, residual_alleles)
        and (
            cross_gene_ratio >= residual_min_cross_gene_ratio
            or selected_has_name_prefix(base_selected, introduced_alleles)
        )
    ):
        return True
    if not introduced_alleles:
        return False
    introduced = any(
        allele.startswith(prefix)
        and not any(base_allele.startswith(prefix) for base_allele in base_selected)
        for allele in selected
        for prefix in introduced_alleles
    )
    if not introduced:
        return False
    return (
        introduced_max_cross_gene_ratio <= 0.0
        or cross_gene_ratio < introduced_max_cross_gene_ratio
    )


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
