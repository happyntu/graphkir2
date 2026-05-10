"""Run a preset-backed old-vs-new benchmark planning comparison."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Build the benchmark runner parser."""
    parser = argparse.ArgumentParser(
        description="Run a graphkir benchmark comparison from a JSON preset.",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to a benchmark preset JSON file.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the comparison summary as JSON.",
    )
    parser.add_argument(
        "--output",
        help="Optional path to write a canonical benchmark result JSON artifact.",
    )
    parser.add_argument(
        "--output-dir",
        help="Optional path to write a benchmark result bundle directory.",
    )
    parser.add_argument(
        "--evaluate",
        action="store_true",
        help="Evaluate merged allele output if truth/prediction files are available.",
    )
    parser.add_argument(
        "--execute-legacy",
        action="store_true",
        help="Execute the legacy graphkir command and collect runtime / memory metrics.",
    )
    parser.add_argument(
        "--force-closest-legacy",
        action="store_true",
        help="Allow execution even when the legacy command is only the closest baseline.",
    )
    parser.add_argument(
        "--pred-tsv",
        help="Optional prediction TSV override used for evaluation.",
    )
    parser.add_argument(
        "--multi-map-policy-override",
        choices=["discard", "keep-all", "best-only", "weighted", "margin", "likelihood"],
        help="Override the policy label stored in the result artifact.",
    )
    return parser


def main() -> None:
    """Load the preset and print the old-vs-new comparison summary."""
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root / "src"))

    from graphkir2.benchmark import (
        apply_allele_evaluation,
        apply_execution_metrics,
        build_result_artifact,
        execute_legacy_command,
        evaluate_allele_calls,
        load_benchmark_preset,
        write_result_bundle,
    )
    from graphkir2.core.pipeline import GraphKir2Pipeline

    args = build_parser().parse_args()
    preset = load_benchmark_preset(args.config)
    pipeline = GraphKir2Pipeline(preset.to_run_config())
    comparison = pipeline.build_comparison_summary()
    if args.multi_map_policy_override:
        notes = list(comparison.legacy_command.notes)
        notes.append(
            "result artifact policy label overridden for a re-run typing ablation; "
            "legacy command is no longer an exact execution match"
        )
        comparison = replace(
            comparison,
            benchmark=replace(
                comparison.benchmark,
                mapping=replace(
                    comparison.benchmark.mapping,
                    multi_map_policy=args.multi_map_policy_override,
                ),
            ),
            legacy_command=replace(
                comparison.legacy_command,
                exact_config_match=False,
                notes=tuple(notes),
            ),
        )
    artifact = build_result_artifact(
        comparison=comparison,
        preset_path=args.config,
        repo_root=repo_root,
    )
    execution = None
    if args.execute_legacy:
        execution = execute_legacy_command(
            comparison.legacy_command,
            repo_root=repo_root,
            allow_closest_match=args.force_closest_legacy,
        )
        if execution.exit_code != 0:
            raise SystemExit(
                "Legacy execution failed with exit code "
                f"{execution.exit_code}. See output logs in the result bundle."
            )
        artifact = apply_execution_metrics(artifact, execution)
    allele_evaluation = None
    if args.evaluate:
        truth_path = preset.allele_truth_tsv
        pred_path = args.pred_tsv or comparison.benchmark.typing.merged_allele_tsv
        if not truth_path:
            raise SystemExit("--evaluate requires allele_truth_tsv in the preset")
        if not Path(pred_path).exists():
            raise SystemExit(f"Predicted allele TSV not found: {pred_path}")
        allele_evaluation = evaluate_allele_calls(truth_path, pred_path)
        artifact = apply_allele_evaluation(artifact, allele_evaluation)
    if args.output_dir:
        write_result_bundle(
            args.output_dir,
            artifact,
            evaluation=allele_evaluation,
            execution=execution,
        )
        print(f"Wrote benchmark bundle to {args.output_dir}")
        return
    if args.output:
        artifact.write_json(args.output)
        print(f"Wrote benchmark artifact to {args.output}")
        return
    if args.json:
        print(artifact.to_json())
        return
    print(comparison.describe())


if __name__ == "__main__":
    main()
