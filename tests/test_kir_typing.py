from graphkir.kir_typing import TypingWithPosNegAllele, _typing_backbone_alias


def test_typing_backbone_alias_maps_kir2dl5_subgenes_to_merged_backbone() -> None:
    assert _typing_backbone_alias("KIR2DL5A*BACKBONE") == "KIR2DL5*BACKBONE"
    assert _typing_backbone_alias("KIR2DL5B*BACKBONE") == "KIR2DL5*BACKBONE"
    assert _typing_backbone_alias("KIR2DS3*BACKBONE") == "KIR2DS3*BACKBONE"


def test_normalize_gene_copy_numbers_merges_kir2dl5_when_only_shared_backbone_exists() -> None:
    model = TypingWithPosNegAllele.__new__(TypingWithPosNegAllele)
    model._gene_reads = {"KIR2DL5*BACKBONE": []}
    model._gene_variants = {"KIR2DL5*BACKBONE": []}

    normalized = model._normalizeGeneCopyNumbers(
        {
            "KIR2DL5A*BACKBONE": 2,
            "KIR2DL5B*BACKBONE": 1,
            "KIR2DS3*BACKBONE": 2,
        }
    )

    assert normalized == {
        "KIR2DL5*BACKBONE": 3,
        "KIR2DS3*BACKBONE": 2,
    }


def test_normalize_gene_copy_numbers_keeps_exact_backbone_when_available() -> None:
    model = TypingWithPosNegAllele.__new__(TypingWithPosNegAllele)
    model._gene_reads = {"KIR2DL5A*BACKBONE": [], "KIR2DL5*BACKBONE": []}
    model._gene_variants = {}

    normalized = model._normalizeGeneCopyNumbers(
        {
            "KIR2DL5A*BACKBONE": 2,
        }
    )

    assert normalized == {"KIR2DL5A*BACKBONE": 2}
