from pathlib import Path
from typing import Any

from graphkir import hisat2


def test_hisat_map_creates_output_parent(monkeypatch, tmp_path: Path) -> None:
    calls: list[tuple[str, list[str] | str]] = []

    def fake_run_tool(name: str, args: list[str], **_: Any) -> None:
        calls.append((name, args))

    def fake_samtobam(name: str, keep: bool = False) -> None:
        calls.append(("samtobam", name))

    monkeypatch.setattr(hisat2, "runTool", fake_run_tool)
    monkeypatch.setattr(hisat2, "samtobam", fake_samtobam)

    output = tmp_path / "nested" / "sample.bam"
    hisat2.hisatMap("index", "reads.1.fq.gz", "reads.2.fq.gz", str(output), threads=3)

    assert output.parent.is_dir()
    assert calls[0][0] == "hisat"
    assert calls[0][1][-2:] == ["-S", str(output.with_suffix(".sam"))]
    assert calls[1] == ("samtobam", str(output.with_suffix("")))
