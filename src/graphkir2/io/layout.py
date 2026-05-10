"""Repository-side layout metadata for graphkir2."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModuleLayout:
    """Named modules expected in the refactored implementation."""

    entries: tuple[str, ...]

    @classmethod
    def default(cls) -> "ModuleLayout":
        """Return the default graphkir2 module layout."""
        return cls(
            entries=(
                "config",
                "io",
                "mapping",
                "cn",
                "typing",
                "benchmark",
                "model",
            )
        )
