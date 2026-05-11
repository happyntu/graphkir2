from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))

from run_functional_stress_sweep import (  # noqa: E402
    DEFAULT_METHODS,
    MethodSpec,
    build_typing_command,
    functional_gains,
    functional_regressions,
    paths_for,
    remaining_functional_errors,
    select_methods,
)


def test_select_methods_preserves_default_order() -> None:
    methods = select_methods(["enhancedgate_geneaware", "discard"])

    assert [method.name for method in methods] == ["discard", "enhancedgate_geneaware"]


def test_paths_for_uses_config_suffix_and_stable_output_names() -> None:
    method = MethodSpec(
        name="enhancedgate_geneaware",
        runner="private_support",
        config_suffix="-conditional-kir2ds3-enhancedgate",
        multi_map_mode="likelihood",
        top_n=5000,
        base_top_n=600,
        gene_base_top_ns="KIR2DL1:1000",
    )

    paths = paths_for("synthetic-x", method, Path("results"))

    assert paths.config == Path(
        "benchmarks/configs/synthetic-x-conditional-kir2ds3-enhancedgate.json"
    )
    assert paths.prediction_tsv == Path("results/synthetic-x.enhancedgate_geneaware.allele.tsv")
    assert paths.bundle_dir == Path("results/synthetic-x.enhancedgate_geneaware.bundle")


def test_build_typing_command_adds_geneaware_overrides() -> None:
    method = [item for item in DEFAULT_METHODS if item.name == "enhancedgate_geneaware"][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert "benchmarks/scripts/rerun_typing_private_support.py" in command
    assert command[-4:] == ["--base-top-n", "600", "--gene-base-top-ns", "KIR2DL1:1000"]


def test_build_typing_command_adds_functional_fallback_overrides() -> None:
    method = [
        item
        for item in DEFAULT_METHODS
        if item.name == "enhancedgate_kir2dl1fallback_geneaware"
    ][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert command[-8:] == [
        "--functional-discard-fallback-genes",
        "KIR2DL1",
        "--functional-discard-fallback-resolution",
        "3",
        "--functional-discard-fallback-max-score",
        "-100.0",
        "--functional-discard-fallback-min-score-delta",
        "20.0",
    ]


def test_build_typing_command_adds_functional_promotion_guard_overrides() -> None:
    method = [
        item
        for item in DEFAULT_METHODS
        if item.name == "enhancedgate_kir2dl1_kir2ds5guard_geneaware"
    ][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert command[-4:] == [
        "--functional-discard-fallback-promoted-alleles",
        "KIR2DS5*027",
        "--functional-discard-fallback-protected-alleles",
        "KIR2DS5*002",
    ]
    assert "--functional-discard-fallback-genes" in command
    assert "KIR2DL1,KIR2DS5" in command


def test_build_typing_command_adds_suballele_guard_overrides() -> None:
    method = [
        item
        for item in DEFAULT_METHODS
        if item.name == "enhancedgate_functionalguard_geneaware"
    ][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert "KIR2DL1,KIR2DS5,KIR2DS3" in command
    assert "KIR2DS5*027,KIR2DS3*00109" in command
    assert "KIR2DS5*002,KIR2DS3*00103" in command


def test_build_typing_command_adds_kir2dl5_overcall_guard_overrides() -> None:
    method = [
        item
        for item in DEFAULT_METHODS
        if item.name == "enhancedgate_kir2dl5guard_geneaware"
    ][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert "--unsupported-overcall-guard-genes" in command
    assert "KIR2DL5" in command
    assert "--unsupported-overcall-guard-min-net-delta" in command
    assert "20.0" in command
    assert "--unsupported-overcall-guard-min-unsupported-delta" in command
    assert "2" in command


def test_build_typing_command_adds_targeted_kir2ds5_overcall_guard_overrides() -> None:
    method = [
        item
        for item in DEFAULT_METHODS
        if item.name == "enhancedgate_kir2dl5_kir2ds5unsupported_geneaware"
    ][0]
    paths = paths_for("synthetic-x", method, Path("results"))

    command = build_typing_command(method, paths)

    assert "--unsupported-overcall-guard-genes" in command
    assert "KIR2DL5" in command
    assert "--targeted-unsupported-overcall-guard-genes" in command
    assert "KIR2DS5" in command
    assert "--targeted-unsupported-overcall-guard-alleles" in command
    assert "KIR2DS5*027,KIR2DS5*010" in command
    assert "--targeted-unsupported-overcall-guard-min-unsupported-delta" in command
    assert "1" in command
    assert "--targeted-unsupported-overcall-guard-preserve-non-target-resolution" in command
    assert "5" in command


def test_functional_regressions_compares_only_three_and_five_digit() -> None:
    rows = [
        {
            "label": "panel",
            "method": "discard",
            "gene": "KIR2DS3",
            "three_digit_f1": "1.000000",
            "five_digit_f1": "0.900000",
            "seven_digit_f1": "0.900000",
        },
        {
            "label": "panel",
            "method": "candidate",
            "gene": "KIR2DS3",
            "three_digit_f1": "0.950000",
            "five_digit_f1": "0.950000",
            "seven_digit_f1": "0.100000",
        },
    ]

    regressions = functional_regressions(rows, "discard", "candidate")

    assert regressions == [
        {
            "label": "panel",
            "gene": "KIR2DS3",
            "metric": "three_digit",
            "baseline_f1": "1.000000",
            "candidate_f1": "0.950000",
            "delta": "-0.050000",
        }
    ]


def test_functional_gains_compares_only_three_and_five_digit() -> None:
    rows = [
        {
            "label": "panel",
            "method": "discard",
            "gene": "KIR2DS3",
            "three_digit_f1": "0.900000",
            "five_digit_f1": "0.900000",
            "seven_digit_f1": "1.000000",
        },
        {
            "label": "panel",
            "method": "candidate",
            "gene": "KIR2DS3",
            "three_digit_f1": "0.950000",
            "five_digit_f1": "0.900000",
            "seven_digit_f1": "0.100000",
        },
    ]

    gains = functional_gains(rows, "discard", "candidate")

    assert gains == [
        {
            "label": "panel",
            "gene": "KIR2DS3",
            "metric": "three_digit",
            "baseline_f1": "0.900000",
            "candidate_f1": "0.950000",
            "delta": "0.050000",
        }
    ]


def test_remaining_functional_errors_ignores_seven_digit() -> None:
    rows = [
        {
            "label": "panel",
            "method": "candidate",
            "gene": "KIR2DL1",
            "three_digit_f1": "1.000000",
            "five_digit_f1": "0.875000",
            "seven_digit_f1": "0.500000",
        }
    ]

    errors = remaining_functional_errors(rows, "candidate")

    assert errors == [
        {
            "label": "panel",
            "gene": "KIR2DL1",
            "metric": "five_digit",
            "f1": "0.875000",
        }
    ]
