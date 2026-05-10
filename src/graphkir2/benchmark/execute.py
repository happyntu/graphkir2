"""Executed benchmark collection helpers."""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from .legacy import LegacyCommandPlan


@dataclass(frozen=True)
class LegacyExecutionResult:
    """Captured metrics and logs for one executed legacy benchmark command."""

    command_shell: str
    exit_code: int
    runtime_seconds: float
    max_rss_mb: float | None
    stdout: str
    stderr: str


def _parse_time_metrics(stderr: str) -> tuple[float | None, float | None]:
    """Extract runtime and max RSS from a `/usr/bin/time` stderr payload."""
    runtime_match = re.search(r"__GK_RUNTIME_SECONDS__=([0-9.]+)", stderr)
    rss_match = re.search(r"__GK_MAXRSS_KB__=([0-9]+)", stderr)
    runtime_seconds = float(runtime_match.group(1)) if runtime_match else None
    max_rss_mb = None
    if rss_match:
        max_rss_mb = float(rss_match.group(1)) / 1024.0
    return runtime_seconds, max_rss_mb


def execute_legacy_command(
    plan: LegacyCommandPlan,
    repo_root: Path,
    allow_closest_match: bool = False,
) -> LegacyExecutionResult:
    """Execute the legacy graphkir command and collect runtime / memory metrics."""
    if not plan.exact_config_match and not allow_closest_match:
        raise ValueError(
            "Legacy command is not an exact config match. "
            "Pass allow_closest_match=True to run the closest baseline anyway."
        )

    command = [plan.program, *plan.args]
    command_shell = plan.render_shell()
    time_binary = shutil.which("time")
    if time_binary == "/usr/bin/time":
        timed_command = [
            time_binary,
            "-f",
            "__GK_RUNTIME_SECONDS__=%e\n__GK_MAXRSS_KB__=%M",
            *command,
        ]
        result = subprocess.run(
            timed_command,
            cwd=repo_root,
            text=True,
            capture_output=True,
            check=False,
        )
        runtime_seconds, max_rss_mb = _parse_time_metrics(result.stderr)
        if runtime_seconds is None:
            runtime_seconds = 0.0
        return LegacyExecutionResult(
            command_shell=command_shell,
            exit_code=result.returncode,
            runtime_seconds=runtime_seconds,
            max_rss_mb=max_rss_mb,
            stdout=result.stdout,
            stderr=result.stderr,
        )

    started = time.monotonic()
    result = subprocess.run(
        command,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )
    runtime_seconds = time.monotonic() - started
    return LegacyExecutionResult(
        command_shell=command_shell,
        exit_code=result.returncode,
        runtime_seconds=runtime_seconds,
        max_rss_mb=None,
        stdout=result.stdout,
        stderr=result.stderr,
    )
