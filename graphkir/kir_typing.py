"""
Typing the sample's alleles using variants file and CN file as input
"""
from typing import Any, Literal
from collections import defaultdict
import json
import re
import math

from .utils import NumpyEncoder, logger
from .msa2hisat import Variant
from .hisat2 import ReadsAndVariantsData, loadReadsAndVariantsData, removeMultipleMapped, PairRead
from .typing_mulit_allele import AlleleTyping, AlleleTypingExonFirst, isHetrozygous
from .typing_em import preprocessHisatReads, hisat2TypingPerGene, printHisatTyping


def groupReads(reads: list[PairRead]) -> dict[str, list[PairRead]]:
    """Group the reads by reference name (i.e. group by gene)"""
    gene_reads = defaultdict(list)
    for read in reads:
        gene_reads[read.backbone].append(read)
    return gene_reads


def groupVariants(variants: list[Variant]) -> dict[str, list[Variant]]:
    """Group the reads by reference name (i.e. group by gene)"""
    gene_variants = defaultdict(list)
    for variant in variants:
        gene_variants[variant.ref].append(variant)
    return gene_variants


def _sam_int_tag(record: str, tag: str) -> int:
    """Extract an integer SAM tag from a record."""
    match = re.search(rf"\t{re.escape(tag)}:i:(-?\d+)", record)
    if match:
        return int(match.group(1))
    return 0


def _sam_mapq(record: str) -> int:
    """Extract MAPQ from a SAM record."""
    fields = record.split("\t")
    if len(fields) > 4:
        try:
            return int(fields[4])
        except ValueError:
            return 0
    return 0


def _pair_read_name(read: PairRead) -> str:
    """Extract the SAM read name for a PairRead."""
    if read.l_sam:
        return read.l_sam.split("\t", 1)[0]
    if read.r_sam:
        return read.r_sam.split("\t", 1)[0]
    return ""


def _pair_read_rank_key(read: PairRead) -> tuple[int, int, int]:
    """Ranking key for selecting the single best alignment of a pair read."""
    align_score = _sam_int_tag(read.l_sam, "AS") + _sam_int_tag(read.r_sam, "AS")
    mapq_score = _sam_mapq(read.l_sam) + _sam_mapq(read.r_sam)
    return (align_score, mapq_score, -read.multiple)


def _pair_read_scalar_score(read: PairRead) -> float:
    """Project the alignment ranking key into a single scalar score."""
    align_score, mapq_score, multiple_penalty = _pair_read_rank_key(read)
    return align_score + mapq_score * 0.1 + multiple_penalty * 0.01


def selectBestMapped(
    reads_data: ReadsAndVariantsData,
    confusion_gene_groups: tuple[frozenset[str], ...] = (),
) -> ReadsAndVariantsData:
    """Keep only one best-scoring alignment per read pair."""
    grouped_reads: dict[str, list[PairRead]] = defaultdict(list)
    variants_map = {str(variant.id): variant for variant in reads_data["variants"]}
    for read in reads_data["reads"]:
        grouped_reads[_pair_read_name(read)].append(read)
    best_reads: list[PairRead] = []
    for group in grouped_reads.values():
        if not group:
            continue
        group_genes = {_backbone_gene(read.backbone) for read in group}
        matched_group = next(
            (
                confusion_group
                for confusion_group in confusion_gene_groups
                if len(group_genes & confusion_group) >= 2
                and group_genes <= confusion_group
            ),
            None,
        )
        if matched_group is not None:
            best_reads.append(
                max(
                    group,
                    key=lambda read: (
                        _read_exon_evidence(read, variants_map),
                        *_pair_read_rank_key(read),
                    ),
                )
            )
        else:
            best_reads.append(max(group, key=_pair_read_rank_key))
    return {
        "variants": reads_data["variants"],
        "reads": best_reads,
    }


def marginWeightMapped(
    reads_data: ReadsAndVariantsData,
    score_gap: float = 5.0,
    score_scale: float = 2.0,
) -> ReadsAndVariantsData:
    """
    Keep the best alignment when the winner is clear, otherwise softly share
    evidence across near-tie alignments.
    """
    grouped_reads: dict[str, list[PairRead]] = defaultdict(list)
    for read in reads_data["reads"]:
        grouped_reads[_pair_read_name(read)].append(read)

    margin_reads: list[PairRead] = []
    for group in grouped_reads.values():
        if not group:
            continue
        if len(group) == 1:
            group[0].weight = 1.0
            margin_reads.append(group[0])
            continue

        backbones = {read.backbone for read in group}
        if len(backbones) == 1:
            best_read = max(group, key=_pair_read_rank_key)
            best_read.weight = 1.0
            margin_reads.append(best_read)
            continue

        scored_group = sorted(
            ((_pair_read_scalar_score(read), read) for read in group),
            key=lambda item: item[0],
            reverse=True,
        )
        best_score = scored_group[0][0]
        second_score = scored_group[1][0]
        if best_score - second_score > score_gap:
            scored_group[0][1].weight = 1.0
            margin_reads.append(scored_group[0][1])
            continue

        near_tie = [
            (score, read)
            for score, read in scored_group
            if best_score - score <= score_gap
        ]
        raw_weights = [
            math.exp((score - best_score) / score_scale) for score, _ in near_tie
        ]
        total = sum(raw_weights)
        if not total:
            uniform = 1.0 / len(near_tie)
            for _, read in near_tie:
                read.weight = uniform
                margin_reads.append(read)
            continue
        for (_, read), weight in zip(near_tie, raw_weights):
            read.weight = weight / total
            margin_reads.append(read)

    return {
        "variants": reads_data["variants"],
        "reads": margin_reads,
    }


def weightMultipleMapped(
    reads_data: ReadsAndVariantsData,
    score_scale: float = 5.0,
) -> ReadsAndVariantsData:
    """Keep all alignments but assign a normalized weight per read pair."""
    grouped_reads: dict[str, list[PairRead]] = defaultdict(list)
    for read in reads_data["reads"]:
        grouped_reads[_pair_read_name(read)].append(read)

    weighted_reads: list[PairRead] = []
    for group in grouped_reads.values():
        if not group:
            continue
        if len(group) == 1:
            group[0].weight = 1.0
            weighted_reads.append(group[0])
            continue
        scores = [
            score_key[0] + score_key[1] * 0.1
            for score_key in map(_pair_read_rank_key, group)
        ]
        max_score = max(scores)
        raw_weights = [math.exp((score - max_score) / score_scale) for score in scores]
        total = sum(raw_weights)
        if not total:
            uniform = 1.0 / len(group)
            for read in group:
                read.weight = uniform
                weighted_reads.append(read)
            continue
        for read, weight in zip(group, raw_weights):
            read.weight = weight / total
            weighted_reads.append(read)
    return {
        "variants": reads_data["variants"],
        "reads": weighted_reads,
    }


def likelihoodAmbiguousMapped(
    reads_data: ReadsAndVariantsData,
    score_scale: float = 5.0,
) -> ReadsAndVariantsData:
    """
    Keep ambiguous alignments and defer their uncertainty to the typing likelihood.

    The per-alignment weight is the prior probability that the physical read
    belongs to that alignment. The typing model later mixes each read group's
    allele probability with neutral probability for the non-target mass.
    """
    grouped_reads: dict[str, list[PairRead]] = defaultdict(list)
    for read in reads_data["reads"]:
        grouped_reads[_pair_read_name(read)].append(read)

    likelihood_reads: list[PairRead] = []
    for group in grouped_reads.values():
        if not group:
            continue
        if len(group) == 1:
            group[0].weight = 1.0
            group[0].ambiguous_weight = 0.0
            likelihood_reads.append(group[0])
            continue

        scores = [_pair_read_scalar_score(read) for read in group]
        max_score = max(scores)
        raw_weights = [math.exp((score - max_score) / score_scale) for score in scores]
        total = sum(raw_weights)
        if not total:
            uniform = 1.0 / len(group)
            for read in group:
                read.weight = uniform
                read.ambiguous_weight = 1.0 - uniform
                likelihood_reads.append(read)
            continue
        for read, weight in zip(group, raw_weights):
            normalized_weight = weight / total
            read.weight = normalized_weight
            read.ambiguous_weight = 1.0 - normalized_weight
            likelihood_reads.append(read)

    return {
        "variants": reads_data["variants"],
        "reads": likelihood_reads,
    }


def parseGeneWeightSpec(spec: str) -> dict[str, float]:
    """Parse a simple `GENE:WEIGHT,...` specification."""
    mapping: dict[str, float] = {}
    if not spec.strip():
        return mapping
    for field in spec.split(","):
        item = field.strip()
        if not item:
            continue
        gene, sep, raw_weight = item.partition(":")
        if not sep:
            raise ValueError(f"Invalid gene weight item: {item}")
        mapping[gene.strip()] = float(raw_weight.strip())
    return mapping


def parseConfusionGeneGroups(spec: str) -> tuple[frozenset[str], ...]:
    """Parse `GENE1/GENE2,GENE3/GENE4/GENE5` into confusion groups."""
    groups: list[frozenset[str]] = []
    if not spec.strip():
        return tuple()
    for group_field in spec.split(","):
        group = frozenset(
            gene.strip() for gene in group_field.split("/") if gene.strip()
        )
        if len(group) >= 2:
            groups.append(group)
    return tuple(groups)


def _backbone_gene(backbone: str) -> str:
    """Convert `KIRX*BACKBONE` into `KIRX`."""
    return backbone.split("*")[0]


def _read_exon_evidence(read: PairRead, variants_map: dict[str, Variant]) -> int:
    """Count exon-tagged variant observations for one aligned pair read."""
    exon_hits = 0
    for variant_id in read.lpv + read.lnv + read.rpv + read.rnv:
        variant = variants_map.get(variant_id)
        if variant is not None and variant.in_exon:
            exon_hits += 1
    return exon_hits


def _typing_backbone_alias(gene: str) -> str:
    """Return the merged graph backbone used for genes represented jointly."""
    if gene.startswith("KIR2DL5A*") or gene.startswith("KIR2DL5B*"):
        return "KIR2DL5*BACKBONE"
    return gene


class Typing:
    """Abstract class for typing allele"""

    def __init__(self) -> None:
        """Read sample's variants or allele abundance result"""
        self._result: dict[str, Any] = {}

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
        """Typing for allele within gene with given CN"""
        raise NotImplementedError

    def typing(self, gene_cn: dict[str, int], min_reads_num: int = 100) -> tuple[list[str], list[str]]:
        """
        Allele Typing for all genes by given CN value per gene (`gene_cn`).

        Returns:
          list of alleles in this sample
        """
        predict_alleles = []
        warning_genes = []
        for gene, cn in gene_cn.items():
            if not cn:
                continue
            # debug
            # if "KIR3DP1" not in gene:
            #     continue
            predict_allele, reads_num = self.typingPerGene(gene, cn)
            predict_alleles.extend(predict_allele)
            if reads_num < min_reads_num:  # 100(reads) * 300(base/read) ~= 5x (min 5k length)
                warning_genes.append(gene)

        return predict_alleles, warning_genes

    def save(self, filename: str) -> None:
        """Save data in file"""
        with open(filename, "w") as f:
            json.dump(self._result, f, cls=NumpyEncoder)

    def getAllPossibleTyping(self) -> list[dict[Any, Any]]:
        """
        Return all possible set of allele typing.
        Call typing() before calling this function
        """
        raise NotImplementedError


class TypingWithPosNegAllele(Typing):
    """Our proposed allele typing method"""

    def __init__(
        self,
        filename_variant_json: str,
        top_n: int = 300,
        multiple: bool = False,
        multi_map_mode: Literal[
            "discard",
            "keep-all",
            "best-only",
            "weighted",
            "margin",
            "likelihood",
        ] = "discard",
        exon_first: bool = False,
        exon_only: bool = False,
        exon_candidate_threshold: float = .9,
        variant_correction: bool = False,
        exon_weight: float = 1.0,
        gene_exon_weights: dict[str, float] | None = None,
        confusion_gene_groups: tuple[frozenset[str], ...] = (),
        margin_score_gap: float = 5.0,
        margin_score_scale: float = 2.0,
        ambiguity_neutral_prob: float = 0.999,
        ambiguity_target_weight_power: float = 1.0,
        select_min_fraction_ratio: float = 0.5,
    ):
        """Read all reads and variants from the json file (From .hisat2.py)"""
        super().__init__()
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        if multiple:
            multi_map_mode = "keep-all"
        if multi_map_mode == "discard":
            reads_data = removeMultipleMapped(reads_data)
        elif multi_map_mode == "best-only":
            reads_data = selectBestMapped(
                reads_data,
                confusion_gene_groups=confusion_gene_groups,
            )
        elif multi_map_mode == "weighted":
            reads_data = weightMultipleMapped(reads_data)
        elif multi_map_mode == "margin":
            reads_data = marginWeightMapped(
                reads_data,
                score_gap=margin_score_gap,
                score_scale=margin_score_scale,
            )
        elif multi_map_mode == "likelihood":
            reads_data = likelihoodAmbiguousMapped(reads_data)
        elif multi_map_mode != "keep-all":
            raise ValueError(f"Unsupported multi_map_mode: {multi_map_mode}")
        self._top_n = top_n
        self._gene_reads = groupReads(reads_data["reads"])
        self._gene_variants = groupVariants(reads_data["variants"])
        self._exon_first = exon_first
        self._exon_only = exon_only
        self._exon_candidate_threshold = exon_candidate_threshold
        self._variant_correction = variant_correction
        self._exon_weight = exon_weight
        self._gene_exon_weights = gene_exon_weights or {}
        self._ambiguity_likelihood = multi_map_mode == "likelihood"
        self._ambiguity_neutral_prob = ambiguity_neutral_prob
        self._ambiguity_target_weight_power = ambiguity_target_weight_power
        self._select_min_fraction_ratio = select_min_fraction_ratio

    def _normalizeGeneCopyNumbers(self, gene_cn: dict[str, int]) -> dict[str, int]:
        """Map CN hints onto the backbones actually present in the typing graph."""
        available_genes = set(self._gene_reads) | set(self._gene_variants)
        normalized_cn: dict[str, int] = defaultdict(int)
        for gene, cn in gene_cn.items():
            if not cn:
                continue
            if gene in available_genes:
                normalized_cn[gene] += cn
                continue

            alias_gene = _typing_backbone_alias(gene)
            if alias_gene in available_genes:
                normalized_cn[alias_gene] += cn
            else:
                normalized_cn[gene] += cn
        return dict(normalized_cn)

    def typing(self, gene_cn: dict[str, int], min_reads_num: int = 100) -> tuple[list[str], list[str]]:
        """Normalize CN keys before running per-backbone allele typing."""
        return super().typing(self._normalizeGeneCopyNumbers(gene_cn), min_reads_num)

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
        """Select reads belonged to the gene and typing it"""
        logger.debug(f"[Allele] {gene=} {cn=}")
        force_homo = False if isHetrozygous(gene) else None
        pure_gene = gene.split("*")[0]
        gene_exon_weight = self._gene_exon_weights.get(pure_gene, self._exon_weight)

        if not self._exon_first and not self._exon_only:
            typ = AlleleTyping(
                self._gene_reads[gene],
                self._gene_variants[gene],
                force_homo=force_homo,
                top_n=self._top_n,
                variant_correction=self._variant_correction,
                exon_weight=gene_exon_weight,
                ambiguity_likelihood=self._ambiguity_likelihood,
                ambiguity_neutral_prob=self._ambiguity_neutral_prob,
                ambiguity_target_weight_power=self._ambiguity_target_weight_power,
            )
        else:
            typ = AlleleTypingExonFirst(
                self._gene_reads[gene],
                self._gene_variants[gene],
                force_homo=force_homo,
                top_n=self._top_n,
                exon_only=self._exon_only,
                candidate_set_threshold=self._exon_candidate_threshold,
                exon_weight=gene_exon_weight,
                ambiguity_likelihood=self._ambiguity_likelihood,
                ambiguity_neutral_prob=self._ambiguity_neutral_prob,
                ambiguity_target_weight_power=self._ambiguity_target_weight_power,
            )
        res = typ.typing(cn)
        self._result[gene] = typ.result
        # return res.selectBest(filter_minor=True)
        alleles = res.selectBest(min_fraction_ratio=self._select_min_fraction_ratio)
        # KIR2DL1*BACKBONE -> KIR2DL1
        alleles = [i if i != "fail" else f"{pure_gene}*" for i in alleles]
        return alleles, typ.getReadsNum()

    def getAllPossibleTyping(self) -> list[dict[Any, Any]]:
        """
        Return all possible set of allele typing.
        Call typing() before calling this function
        """
        possible_list = []
        for gene, result in self._result.items():
            for rank, (value, alleles) in enumerate(result[-1].selectAllPossible(.9)):
                simple_result = {
                    "gene": gene,
                    "rank":  rank,
                    "value": value,
                }
                for i, allele in enumerate(alleles):
                    simple_result[str(i + 1)] = allele
                possible_list.append(simple_result)
        return possible_list


class TypingWithReport(Typing):
    """Typing alleles from report by abundance value"""

    def __init__(self, filename_variant_json: str):
        """Read report files (json)"""
        super().__init__()
        reads_data = loadReadsAndVariantsData(filename_variant_json)
        reads_data = removeMultipleMapped(reads_data)
        self._gene_reads = preprocessHisatReads(reads_data)

    def typingPerGene(self, gene: str, cn: int) -> tuple[list[str], int]:
        """
        Typing for allele within gene with given CN

        Example:
          * case1: a1=0.66 a2=0.33
              * CN=1: a1
              * CN=2: a1 a1
              * CN=3: a1 a1 a2
          * case2: a1=0.90 a2=0.10
              * CN=1: a1
              * CN=2: a1 a1
              * CN=3: a1 a1 a1
        """
        # main typing
        report = hisat2TypingPerGene(self._gene_reads[gene])
        report = sorted(report, key=lambda i: -i.prob)

        # main calling logic wrote here
        est_prob = 1 / cn
        called_alleles = []
        for allele in report:
            pred_count = max(1, round(allele.prob / est_prob))
            for _ in range(min(cn, pred_count)):
                called_alleles.append(allele.allele)
            allele.cn = pred_count

            cn -= pred_count
            if cn <= 0:
                break

        self._result[gene] = report
        return called_alleles, len(self._gene_reads[gene])

    def save(self, filename: str) -> None:
        """save additional report txt"""
        super().save(filename)
        name = filename
        if filename.endswith(".json"):
            name = filename[:-5]
        with open(name + ".txt", "w") as f:
            printHisatTyping(self._result, file=f)


def selectKirTypingModel(
    method: str,
    filename_variant_json: str,
    **kwargs: Any,
) -> Typing:
    """Select and Init typing model"""
    if method == "full":
        return TypingWithPosNegAllele(filename_variant_json, **kwargs)
    if method.startswith("exonfirst"):
        fields = method.split("_")
        threshold = 0.0
        if len(fields) == 2:  # e.g. exonfirst_1.2
            threshold = float(method[len("exonfirst_") :])
        return TypingWithPosNegAllele(
            filename_variant_json,
            exon_first=True,
            exon_candidate_threshold=threshold,
            **kwargs,
        )
    if method == "em":
        return TypingWithReport(filename_variant_json)
    raise NotImplementedError
