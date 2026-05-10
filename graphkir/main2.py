"""Bootstrap entrypoint for the refactored graphkir2 implementation."""

from __future__ import annotations

import sys
from pathlib import Path


def entrypoint() -> None:
    """Add the src directory to sys.path and dispatch to graphkir2."""
    repo_root = Path(__file__).resolve().parent.parent
    src_dir = repo_root / "src"
    sys.path.insert(0, str(src_dir))
    from graphkir2.cli import entrypoint as graphkir2_entrypoint

    graphkir2_entrypoint()
