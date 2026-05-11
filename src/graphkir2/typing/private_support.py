"""Private-variant support helpers for functional typing reranking."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any


VariantSupport = defaultdict[str, float]
GeneGroups = tuple[frozenset[str], ...]
NameSet = frozenset[str]
GeneTopN = dict[str, int]


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


def parse_gene_top_n_spec(spec: str) -> GeneTopN:
    """Parse comma-separated `GENE:TOPN` overrides."""
    gene_top_n: GeneTopN = {}
    for field in spec.split(","):
        item = field.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid gene top-n field: {item}")
        gene, value = (part.strip() for part in item.split(":", 1))
        if not gene:
            raise ValueError(f"Invalid empty gene in top-n field: {item}")
        top_n = int(value)
        if top_n <= 0:
            raise ValueError(f"Gene top-n must be positive: {item}")
        gene_top_n[gene] = top_n
    return gene_top_n


def choose_targeted_top_n(
    gene_name: str,
    top_n: int,
    base_top_n: int | None,
    target_genes: NameSet,
    gene_base_top_ns: GeneTopN | None = None,
) -> int:
    """Use high top-n only for target genes when a lower base top-n is supplied."""
    if gene_name in target_genes:
        return top_n
    if gene_base_top_ns and gene_name in gene_base_top_ns:
        return min(top_n, gene_base_top_ns[gene_name])
    if base_top_n is None:
        return top_n
    return base_top_n


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


def _functional_resolution_key(
    allele_names: list[str],
    resolution: int,
) -> tuple[str, ...]:
    """Build a sorted allele key at one functional resolution."""
    if not allele_names or any("*" not in allele for allele in allele_names):
        return ()
    return tuple(sorted(limit_allele_resolution(allele, resolution) for allele in allele_names))


def should_use_functional_discard_fallback(
    selected: list[str],
    discard_selected: list[str],
    selected_private_support: float,
    discard_private_support: float,
    resolution: int = 3,
    min_score_delta: float = 0.0,
) -> bool:
    """Use discard evidence when likelihood changes a weak functional call."""
    selected_key = _functional_resolution_key(selected, resolution)
    discard_key = _functional_resolution_key(discard_selected, resolution)
    if not selected_key or not discard_key or len(selected_key) != len(discard_key):
        return False
    if selected_key == discard_key:
        return False
    return discard_private_support >= selected_private_support + min_score_delta


def _name_prefix_counts(allele_names: list[str], prefixes: NameSet) -> Counter[str]:
    """Count selected alleles by configured allele-name prefix."""
    counts: Counter[str] = Counter()
    for allele in allele_names:
        for prefix in prefixes:
            if allele.startswith(prefix):
                counts[prefix] += 1
    return counts


def apply_functional_promotion_guard(
    selected: list[str],
    discard_selected: list[str],
    promoted_alleles: NameSet,
    protected_alleles: NameSet,
    resolution: int = 3,
) -> list[str]:
    """Replace configured promoted calls with protected discard evidence.

    This guard is intentionally narrower than a full discard fallback. It only
    fires when likelihood introduced a configured promoted allele prefix and
    discard has extra copies of a configured protected prefix. When possible,
    it reuses an already selected protected allele suffix to avoid unnecessary
    7-digit churn while fixing the targeted functional or suballele class.
    """
    if not promoted_alleles or not protected_alleles:
        return selected
    selected_key = _functional_resolution_key(selected, resolution)
    discard_key = _functional_resolution_key(discard_selected, resolution)
    if not selected_key or not discard_key or len(selected_key) != len(discard_key):
        return selected

    selected_promoted = _name_prefix_counts(selected, promoted_alleles)
    discard_promoted = _name_prefix_counts(discard_selected, promoted_alleles)
    selected_protected = _name_prefix_counts(selected, protected_alleles)
    discard_protected = _name_prefix_counts(discard_selected, protected_alleles)

    promoted_shortlist = [
        prefix
        for prefix, count in selected_promoted.items()
        if count > discard_promoted[prefix]
    ]
    protected_shortlist = [
        prefix
        for prefix, count in discard_protected.items()
        if count > selected_protected[prefix]
    ]
    if not promoted_shortlist or not protected_shortlist:
        return selected

    guarded = list(selected)
    promoted_indices: list[int] = []
    remaining_promoted = {
        prefix: selected_promoted[prefix] - discard_promoted[prefix]
        for prefix in promoted_shortlist
    }
    for index, allele in enumerate(selected):
        for prefix in promoted_shortlist:
            if remaining_promoted[prefix] > 0 and allele.startswith(prefix):
                promoted_indices.append(index)
                remaining_promoted[prefix] -= 1
                break

    replacements: list[str] = []
    for prefix in protected_shortlist:
        needed = discard_protected[prefix] - selected_protected[prefix]
        selected_options = [allele for allele in selected if allele.startswith(prefix)]
        discard_options = [allele for allele in discard_selected if allele.startswith(prefix)]
        options = selected_options or discard_options
        if not options:
            continue
        for index in range(needed):
            replacements.append(options[min(index, len(options) - 1)])

    if not promoted_indices or not replacements:
        return selected
    for index, replacement in zip(promoted_indices, replacements):
        guarded[index] = replacement
    return guarded


def _genotype_variant_ids(
    allele_names: list[str],
    allele_variants: dict[str, set[str]],
) -> set[str]:
    """Return the union of variant IDs carried by one allele tuple."""
    variant_ids: set[str] = set()
    for allele in allele_names:
        variant_ids.update(allele_variants.get(allele, set()))
    return variant_ids


def _genotype_variant_ids_for_prefixes(
    allele_names: list[str],
    allele_variants: dict[str, set[str]],
    prefixes: NameSet,
) -> set[str]:
    """Return genotype variant IDs carried by alleles matching prefixes."""
    if not prefixes:
        return _genotype_variant_ids(allele_names, allele_variants)
    variant_ids: set[str] = set()
    for allele in allele_names:
        if any(allele.startswith(prefix) for prefix in prefixes):
            variant_ids.update(allele_variants.get(allele, set()))
    return variant_ids


def _non_prefix_resolution_key(
    allele_names: list[str],
    prefixes: NameSet,
    resolution: int,
) -> tuple[str, ...]:
    """Build a functional key for alleles that are not covered by prefixes."""
    return tuple(
        sorted(
            limit_allele_resolution(allele, resolution)
            for allele in allele_names
            if "*" in allele
            and not any(allele.startswith(prefix) for prefix in prefixes)
        )
    )


def unsupported_candidate_only_evidence(
    candidate: list[str],
    alternative: list[str],
    allele_variants: dict[str, set[str]],
    positive: VariantSupport,
    negative: VariantSupport,
    negative_threshold: float = 5.0,
    max_positive: float = 1.0,
    candidate_allele_prefixes: NameSet = frozenset(),
) -> tuple[int, float]:
    """Measure unsupported variants introduced by `candidate` vs `alternative`.

    The KIR2DL5 failure audit showed wrong calls that add candidate-only
    variants contradicted by negative reads. This helper scores that overcall
    evidence without using sample truth labels.
    """
    candidate_only = _genotype_variant_ids(candidate, allele_variants) - _genotype_variant_ids(
        alternative,
        allele_variants,
    )
    if candidate_allele_prefixes:
        candidate_only &= _genotype_variant_ids_for_prefixes(
            candidate,
            allele_variants,
            candidate_allele_prefixes,
        )
    unsupported = 0
    net_penalty = 0.0
    for variant_id in candidate_only:
        positive_support = positive[variant_id]
        negative_support = negative[variant_id]
        net_penalty += max(0.0, negative_support - positive_support)
        if negative_support >= negative_threshold and positive_support <= max_positive:
            unsupported += 1
    return unsupported, net_penalty


def _passes_fraction_filter(
    result: object,
    index: int,
    min_fraction_ratio: float,
) -> bool:
    """Return whether a result row passes the standard abundance filter."""
    expect_prob = 1 / int(result.n)  # type: ignore[attr-defined]
    return all(
        fraction >= expect_prob * min_fraction_ratio
        for fraction in result.fraction[index]  # type: ignore[attr-defined]
    )


def select_against_unsupported_candidate_only_variants(
    result: object,
    selected: list[str],
    allele_variants: dict[str, set[str]],
    positive: VariantSupport,
    negative: VariantSupport,
    min_fraction_ratio: float,
    max_likelihood_gap: float = 25.0,
    min_unsupported_delta: int = 2,
    min_net_delta: float = 20.0,
    negative_threshold: float = 5.0,
    max_positive: float = 1.0,
    selected_allele_prefixes: NameSet = frozenset(),
    preserve_non_target_resolution: int = 0,
) -> list[str]:
    """Fallback from an overcalled genotype to a nearby less-unsupported one.

    This is intentionally unsupervised: it compares the selected genotype with
    nearby likelihood candidates and only switches when the selected call adds
    substantially more unsupported candidate-only variants than the alternative.
    """
    if result.isFail():  # type: ignore[attr-defined]
        return selected
    if selected_allele_prefixes and not selected_has_name_prefix(
        selected,
        selected_allele_prefixes,
    ):
        return selected
    selected_non_target_counts: Counter[str] = Counter()
    if selected_allele_prefixes and preserve_non_target_resolution > 0:
        selected_non_target_counts = Counter(
            _non_prefix_resolution_key(
                selected,
                selected_allele_prefixes,
                preserve_non_target_resolution,
            )
        )

    selected_key = sorted(selected)
    selected_index = 0
    for index, allele_names in enumerate(result.allele_name):  # type: ignore[attr-defined]
        if sorted(allele_names) == selected_key:
            selected_index = index
            break
    selected_value = float(result.value[selected_index])  # type: ignore[attr-defined]

    best: tuple[float, int, float, list[str]] | None = None
    for index, alternative in enumerate(result.allele_name):  # type: ignore[attr-defined]
        if sorted(alternative) == selected_key:
            continue
        likelihood_gap = selected_value - float(result.value[index])  # type: ignore[attr-defined]
        if likelihood_gap < 0.0 or likelihood_gap > max_likelihood_gap:
            continue
        if not _passes_fraction_filter(result, index, min_fraction_ratio):
            continue
        if selected_non_target_counts:
            alternative_non_target_counts = Counter(
                _non_prefix_resolution_key(
                    alternative,
                    selected_allele_prefixes,
                    preserve_non_target_resolution,
                )
            )
            if any(
                alternative_non_target_counts[key] < count
                for key, count in selected_non_target_counts.items()
            ):
                continue
        selected_unsupported, selected_penalty = unsupported_candidate_only_evidence(
            selected,
            alternative,
            allele_variants,
            positive,
            negative,
            negative_threshold,
            max_positive,
            selected_allele_prefixes,
        )
        alternative_unsupported, alternative_penalty = unsupported_candidate_only_evidence(
            alternative,
            selected,
            allele_variants,
            positive,
            negative,
            negative_threshold,
            max_positive,
            selected_allele_prefixes,
        )
        unsupported_delta = selected_unsupported - alternative_unsupported
        net_delta = selected_penalty - alternative_penalty
        if unsupported_delta < min_unsupported_delta or net_delta < min_net_delta:
            continue
        candidate = (net_delta, unsupported_delta, -likelihood_gap, alternative)
        if best is None or candidate > best:
            best = candidate

    if best is None:
        return selected
    return best[3]


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
