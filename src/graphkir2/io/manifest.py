"""Sample-manifest loading for graphkir2."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ..model import SampleInput


REQUIRED_COLUMNS = ("name", "r1", "r2")
OPTIONAL_COLUMNS = ("cnfile",)


@dataclass(frozen=True)
class SampleManifest:
    """Collection of parsed sample inputs."""

    path: str
    samples: tuple[SampleInput, ...]

    def describe(self) -> str:
        """Return a compact textual summary of the manifest."""
        lines = [
            f"manifest={self.path}",
            f"samples={len(self.samples)}",
        ]
        lines.extend(
            f"- {sample.sample_id}: {sample.read1} | {sample.read2}"
            for sample in self.samples
        )
        return "\n".join(lines)


def load_sample_manifest(path: str) -> SampleManifest:
    """Load the standard Graph-KIR sample manifest CSV."""
    manifest_path = Path(path)
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Sample manifest is missing a header: {path}")

        columns = set(reader.fieldnames)
        missing = [column for column in REQUIRED_COLUMNS if column not in columns]
        if missing:
            raise ValueError(
                f"Sample manifest {path} is missing required columns: {missing}"
            )

        samples: list[SampleInput] = []
        for row_id, row in enumerate(reader, start=2):
            output_prefix = row["name"].strip()
            read1 = row["r1"].strip()
            read2 = row["r2"].strip()
            if not output_prefix or not read1 or not read2:
                raise ValueError(
                    f"Sample manifest {path} has empty required fields on line {row_id}"
                )
            samples.append(
                SampleInput(
                    sample_id=Path(output_prefix).name,
                    output_prefix=output_prefix,
                    read1=read1,
                    read2=read2,
                    copy_number_hint=row.get("cnfile", "").strip(),
                )
            )
    return SampleManifest(path=str(manifest_path), samples=tuple(samples))
