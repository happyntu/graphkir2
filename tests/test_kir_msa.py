from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from graphkir import kir_msa


def _record(name: str, sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id=name, description="")


def test_merge_block_to_msa_pads_missing_alleles_before_concat(monkeypatch) -> None:
    monkeypatch.setattr(kir_msa, "kir_block_name", ["5UTR", "exon1", "3UTR"])

    msa = kir_msa.mergeBlockToMsa(
        {
            "5UTR": [
                _record("KIR2DL1*001", "AA"),
                _record("KIR2DS1*001", "CC"),
            ],
            "exon1": [
                _record("KIR2DL1*001", "GG"),
                _record("KIR2DS1*001", "TT"),
            ],
            "3UTR": [
                _record("KIR2DS1*001", "GGG"),
                _record("KIR2DS1*002", "TTT"),
            ],
        }
    )

    assert set(msa.alleles) == {
        "KIR2DL1*001",
        "KIR2DS1*001",
        "KIR2DS1*002",
    }
    assert str(msa.get("KIR2DL1*001")) == "AAGG---"
    assert str(msa.get("KIR2DS1*001")) == "CCTTGGG"
    assert str(msa.get("KIR2DS1*002")) == "----TTT"
