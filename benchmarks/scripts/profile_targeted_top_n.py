"""Profile targeted typing top-n on existing synthetic fixtures."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


DEFAULT_LABELS = (
    "synthetic-difficult5x12",
    "synthetic-functional8x6",
)


@dataclass(frozen=True)
class ProfileVariant:
    """One targeted-top-n profiling variant."""

    name: str
    base_top_n: int
    gene_base_top_ns: str = ""


PROFILE_VARIANTS = (
    ProfileVariant("full_top5000", 0),
    ProfileVariant("targeted_base600", 600),
    ProfileVariant("geneaware_base600_kir2dl1_1000", 600, "KIR2DL1:1000"),
)


@dataclass(frozen=True)
class TimedCommandResult:
    """Command runtime and memory measurement."""

    exit_code: int
    runtime_seconds: float
    max_rss_mb: float | None
    stdout: str
    stderr: str


def parse_time_metrics(stderr: str) -> tuple[float | None, float | None]:
    """Parse `/usr/bin/time` metrics from stderr."""
    runtime_match = re.search(r"__GK_RUNTIME_SECONDS__=([0-9.]+)", stderr)
    rss_match = re.search(r"__GK_MAXRSS_KB__=([0-9]+)", stderr)
    runtime_seconds = float(runtime_match.group(1)) if runtime_match else None
    max_rss_mb = float(rss_match.group(1)) / 1024.0 if rss_match else None
    return runtime_seconds, max_rss_mb


def run_timed(command: list[str], cwd: Path) -> TimedCommandResult:
    """Run a command and collect elapsed time plus peak RSS when available."""
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
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
        )
        runtime_seconds, max_rss_mb = parse_time_metrics(result.stderr)
        return TimedCommandResult(
            exit_code=result.returncode,
            runtime_seconds=runtime_seconds if runtime_seconds is not None else 0.0,
            max_rss_mb=max_rss_mb,
            stdout=result.stdout,
            stderr=result.stderr,
        )

    started = time.monotonic()
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    return TimedCommandResult(
        exit_code=result.returncode,
        runtime_seconds=time.monotonic() - started,
        max_rss_mb=None,
        stdout=result.stdout,
        stderr=result.stderr,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the profiling parser."""
    parser = argparse.ArgumentParser(
        description="Profile full vs targeted typing top-n on synthetic panels.",
    )
    parser.add_argument(
        "--label",
        action="append",
        default=[],
        help="Synthetic label to profile. Defaults to a representative pair.",
    )
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument(
        "--base-top-n",
        action="append",
        type=int,
        default=[],
        help=(
            "Base top-n variants to profile. Repeat for multiple values; "
            "0 means full top-n. Defaults to full, base600, and gene-aware base600."
        ),
    )
    parser.add_argument(
        "--gene-base-top-ns",
        action="append",
        default=[],
        help=(
            "Gene-specific overrides paired by position with --base-top-n; "
            "use an empty string for no override."
        ),
    )
    parser.add_argument(
        "--results-dir",
        default="benchmarks/results/targeted-topn-profile",
    )
    parser.add_argument(
        "--output-tsv",
        default="benchmarks/results/targeted-topn-profile/summary.tsv",
    )
    parser.add_argument(
        "--output-md",
        default="docs/research/2026-05-11_targeted_top_n_runtime_profile.md",
    )
    parser.add_argument(
        "--skip-md",
        action="store_true",
        help="Do not write the markdown report.",
    )
    return parser


def build_profile_variants(base_top_ns: list[int]) -> tuple[ProfileVariant, ...]:
    """Build profile variants from requested base top-n values."""
    if not base_top_ns:
        return PROFILE_VARIANTS
    return tuple(
        ProfileVariant(
            "full_top5000" if base_top_n == 0 else f"targeted_base{base_top_n}",
            base_top_n,
        )
        for base_top_n in base_top_ns
    )


def build_profile_variants_with_gene_overrides(
    base_top_ns: list[int],
    gene_base_top_specs: list[str],
) -> tuple[ProfileVariant, ...]:
    """Build profile variants from base top-n values and optional gene overrides."""
    if not base_top_ns:
        return PROFILE_VARIANTS
    if gene_base_top_specs and len(gene_base_top_specs) != len(base_top_ns):
        raise ValueError("--gene-base-top-ns must be provided once per --base-top-n")
    variants: list[ProfileVariant] = []
    for index, base_top_n in enumerate(base_top_ns):
        gene_spec = gene_base_top_specs[index] if gene_base_top_specs else ""
        if base_top_n == 0:
            name = "full_top5000"
        elif gene_spec:
            safe_gene_spec = gene_spec.replace(":", "_").replace(",", "_")
            name = f"geneaware_base{base_top_n}_{safe_gene_spec.lower()}"
        else:
            name = f"targeted_base{base_top_n}"
        variants.append(ProfileVariant(name, base_top_n, gene_spec))
    return tuple(variants)


def config_for_label(label: str) -> Path:
    """Return the committed enhancedgate config for one synthetic label."""
    return Path("benchmarks/configs") / f"{label}-conditional-kir2ds3-enhancedgate.json"


def prediction_path(results_dir: Path, label: str, variant: ProfileVariant) -> Path:
    """Return the prediction path for a profile run."""
    return results_dir / f"{label}.{variant.name}.allele.tsv"


def bundle_path(results_dir: Path, label: str, variant: ProfileVariant) -> Path:
    """Return the evaluation bundle path for a profile run."""
    return results_dir / f"{label}.{variant.name}.bundle"


def run_profile_variant(
    repo_root: Path,
    label: str,
    variant: ProfileVariant,
    top_n: int,
    results_dir: Path,
) -> dict[str, str]:
    """Run one profiling variant and return a summary row."""
    config = config_for_label(label)
    if not config.exists():
        raise FileNotFoundError(f"Missing config for {label}: {config}")

    pred_tsv = prediction_path(results_dir, label, variant)
    bundle_dir = bundle_path(results_dir, label, variant)
    pred_tsv.parent.mkdir(parents=True, exist_ok=True)

    rerun_command = [
        sys.executable,
        "benchmarks/scripts/rerun_typing_private_support.py",
        "--config",
        str(config),
        "--top-n",
        str(top_n),
        "--output-tsv",
        str(pred_tsv),
    ]
    if variant.base_top_n:
        rerun_command.extend(["--base-top-n", str(variant.base_top_n)])
    if variant.gene_base_top_ns:
        rerun_command.extend(["--gene-base-top-ns", variant.gene_base_top_ns])
    timed = run_timed(rerun_command, repo_root)
    if timed.exit_code != 0:
        raise RuntimeError(
            f"Typing profile failed for {label}/{variant.name}\n"
            f"stdout:\n{timed.stdout}\n"
            f"stderr:\n{timed.stderr}"
        )

    compare_command = [
        sys.executable,
        "benchmarks/scripts/run_compare.py",
        "--config",
        str(config),
        "--evaluate",
        "--pred-tsv",
        str(pred_tsv),
        "--output-dir",
        str(bundle_dir),
    ]
    compare = subprocess.run(
        compare_command,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )
    if compare.returncode != 0:
        raise RuntimeError(
            f"Evaluation failed for {label}/{variant.name}\n"
            f"stdout:\n{compare.stdout}\n"
            f"stderr:\n{compare.stderr}"
        )
    artifact = json.loads((bundle_dir / "compare.json").read_text(encoding="utf-8"))
    metrics = artifact["metrics"]
    return {
        "label": label,
        "variant": variant.name,
        "top_n": str(top_n),
        "base_top_n": str(variant.base_top_n),
        "gene_base_top_ns": variant.gene_base_top_ns,
        "runtime_seconds": f"{timed.runtime_seconds:.3f}",
        "max_rss_mb": "" if timed.max_rss_mb is None else f"{timed.max_rss_mb:.1f}",
        "three_digit_f1": f"{float(metrics['three_digit_f1']):.6f}",
        "five_digit_f1": f"{float(metrics['five_digit_f1']):.6f}",
        "seven_digit_f1": f"{float(metrics['seven_digit_f1']):.6f}",
        "prediction_tsv": str(pred_tsv),
        "bundle_dir": str(bundle_dir),
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write profile rows to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "label",
        "variant",
        "top_n",
        "base_top_n",
        "gene_base_top_ns",
        "runtime_seconds",
        "max_rss_mb",
        "three_digit_f1",
        "five_digit_f1",
        "seven_digit_f1",
        "prediction_tsv",
        "bundle_dir",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fmt_float(value: str, digits: int = 3) -> str:
    """Format a string float for markdown."""
    return f"{float(value):.{digits}f}" if value else ""


def write_markdown(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a concise markdown runtime profile report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Targeted top-n runtime profile",
        "",
        "This profiles `rerun_typing_private_support.py` on existing synthetic fixtures.",
        "`full_top5000` keeps every gene at `top_n=5000`; targeted variants keep",
        "the enhancedgate target genes at `top_n=5000` and use lower `base_top_n`",
        "values for non-target genes; gene-aware variants override individual genes.",
        "",
        "| panel | variant | runtime s | max RSS MB | 3-digit F1 | 5-digit F1 | 7-digit F1 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["label"],
                    row["variant"],
                    _fmt_float(row["runtime_seconds"], 2),
                    _fmt_float(row["max_rss_mb"], 1),
                    _fmt_float(row["three_digit_f1"], 4),
                    _fmt_float(row["five_digit_f1"], 4),
                    _fmt_float(row["seven_digit_f1"], 4),
                ]
            )
            + " |"
        )

    by_label: dict[str, dict[str, dict[str, str]]] = {}
    for row in rows:
        by_label.setdefault(row["label"], {})[row["variant"]] = row
    lines.extend(["", "## Interpretation", ""])
    for label, variants in by_label.items():
        full = variants.get("full_top5000")
        if not full:
            continue
        full_runtime = float(full["runtime_seconds"])
        full_accuracy = (
            full["three_digit_f1"],
            full["five_digit_f1"],
            full["seven_digit_f1"],
        )
        non_regressing = [
            row
            for name, row in variants.items()
            if name != "full_top5000"
            and (
                row["three_digit_f1"],
                row["five_digit_f1"],
                row["seven_digit_f1"],
            )
            == full_accuracy
        ]
        if non_regressing:
            fastest = min(non_regressing, key=lambda row: float(row["runtime_seconds"]))
            speedup = full_runtime / float(fastest["runtime_seconds"])
            lines.append(
                f"* `{label}`: fastest non-regressing targeted variant was "
                f"`{fastest['variant']}` at {float(fastest['runtime_seconds']):.2f}s "
                f"vs {full_runtime:.2f}s ({speedup:.2f}x speedup)."
            )
        regressions = [
            row
            for name, row in variants.items()
            if name != "full_top5000"
            and (
                row["three_digit_f1"],
                row["five_digit_f1"],
                row["seven_digit_f1"],
            )
            != full_accuracy
        ]
        for row in regressions:
            lines.append(
                f"* `{label}`: `{row['variant']}` changed aggregate F1; inspect its "
                "bundle before using it as a default."
            )
    lines.extend(
        [
            "",
            "Generated artifact paths are listed in the TSV summary.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run the targeted top-n profile."""
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    labels = tuple(args.label) if args.label else DEFAULT_LABELS
    variants = build_profile_variants_with_gene_overrides(
        args.base_top_n,
        args.gene_base_top_ns,
    )
    results_dir = Path(args.results_dir)
    rows: list[dict[str, str]] = []
    for label in labels:
        for variant in variants:
            row = run_profile_variant(
                repo_root,
                label,
                variant,
                top_n=args.top_n,
                results_dir=results_dir,
            )
            rows.append(row)
            print(
                f"{label}\t{variant.name}\t{row['runtime_seconds']}s\t"
                f"{row['max_rss_mb']} MB\t"
                f"{row['three_digit_f1']}/{row['five_digit_f1']}/{row['seven_digit_f1']}"
            )
    write_tsv(Path(args.output_tsv), rows)
    if not args.skip_md:
        write_markdown(Path(args.output_md), rows)
        print(f"Wrote markdown report to {args.output_md}")
    print(f"Wrote TSV summary to {args.output_tsv}")


if __name__ == "__main__":
    main()
