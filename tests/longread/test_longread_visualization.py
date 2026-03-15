"""Tests for longread workflow orchestrator and pipeline config not covered by test_longread.py.

Tests LongReadOrchestrator, PipelineStep, PipelineResult, and pipeline configuration.
All tests use real implementations -- NO MOCKING.

Note: Visualization functions are fully tested in test_longread.py (TestVisualization class).
This file tests the workflow orchestration layer that drives those visualizations.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import pytest

from metainformant.longread.workflow.orchestrator import (
    LongReadOrchestrator,
    PipelineResult,
    PipelineStep,
)
from metainformant.longread.workflow.pipelines import (
    get_assembly_pipeline_config,
    get_methylation_pipeline_config,
    get_qc_pipeline_config,
    get_sv_pipeline_config,
    validate_pipeline_config,
)


def _random_dna(length: int, seed: int = 42) -> str:
    """Generate a deterministic random DNA sequence."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def _phred_string(scores: list[int]) -> str:
    """Convert a list of Phred scores to an ASCII quality string (Phred+33)."""
    return "".join(chr(s + 33) for s in scores)


# ---------------------------------------------------------------------------
# PipelineStep dataclass
# ---------------------------------------------------------------------------


class TestPipelineStep:
    """Tests for PipelineStep dataclass."""

    def test_pipeline_step_defaults(self) -> None:
        """PipelineStep has correct defaults."""
        step = PipelineStep(name="test_step", function=lambda ctx: None)
        assert step.name == "test_step"
        assert step.params == {}
        assert step.depends_on == []
        assert step.status == "pending"
        assert step.result is None
        assert step.duration_seconds == 0.0
        assert step.error == ""

    def test_pipeline_step_custom(self) -> None:
        """PipelineStep stores custom values."""
        step = PipelineStep(
            name="filter",
            function=lambda ctx: [1, 2, 3],
            params={"min_length": 1000},
            depends_on=["load"],
            status="completed",
            result=[1, 2, 3],
            duration_seconds=1.5,
        )
        assert step.name == "filter"
        assert step.params["min_length"] == 1000
        assert step.depends_on == ["load"]
        assert step.status == "completed"
        assert step.result == [1, 2, 3]
        assert step.duration_seconds == 1.5


# ---------------------------------------------------------------------------
# PipelineResult dataclass
# ---------------------------------------------------------------------------


class TestPipelineResult:
    """Tests for PipelineResult dataclass."""

    def test_pipeline_result_defaults(self) -> None:
        """PipelineResult stores pipeline execution state."""
        result = PipelineResult(
            pipeline_name="qc",
            steps=[],
            success=True,
            total_duration=5.0,
            output_dir=Path("/tmp/output"),
        )
        assert result.pipeline_name == "qc"
        assert result.steps == []
        assert result.success is True
        assert result.total_duration == 5.0
        assert result.output_dir == Path("/tmp/output")
        assert result.summary == {}

    def test_pipeline_result_with_steps(self) -> None:
        """PipelineResult tracks steps with results."""
        steps = [
            PipelineStep(name="step1", function=lambda ctx: None, status="completed", duration_seconds=1.0),
            PipelineStep(name="step2", function=lambda ctx: None, status="failed", error="boom"),
        ]
        result = PipelineResult(
            pipeline_name="test",
            steps=steps,
            success=False,
            total_duration=3.0,
            output_dir=Path("/tmp"),
            summary={"step1": {"data": True}},
        )
        assert len(result.steps) == 2
        assert result.steps[0].status == "completed"
        assert result.steps[1].status == "failed"
        assert result.steps[1].error == "boom"
        assert result.summary["step1"]["data"] is True


# ---------------------------------------------------------------------------
# Pipeline configuration functions
# ---------------------------------------------------------------------------


class TestPipelineConfigs:
    """Tests for pipeline configuration factory functions."""

    def test_qc_pipeline_config_defaults(self) -> None:
        """get_qc_pipeline_config returns valid config with defaults."""
        config = get_qc_pipeline_config()
        assert config["pipeline_name"] == "qc"
        assert config["parameters"]["min_length"] == 1000
        assert config["parameters"]["min_quality"] == 7.0
        assert config["parameters"]["trim_adapters"] is True
        assert len(config["steps"]) > 0

    def test_qc_pipeline_config_custom(self) -> None:
        """get_qc_pipeline_config accepts custom parameters."""
        config = get_qc_pipeline_config(
            min_length=2000,
            min_quality=10.0,
            trim_adapters=False,
        )
        assert config["parameters"]["min_length"] == 2000
        assert config["parameters"]["min_quality"] == 10.0
        assert config["parameters"]["trim_adapters"] is False

    def test_assembly_pipeline_config(self) -> None:
        """get_assembly_pipeline_config returns valid config."""
        config = get_assembly_pipeline_config(
            min_overlap=3000,
            k=17,
            polish_iterations=3,
        )
        assert config["pipeline_name"] == "assembly"
        assert config["parameters"]["min_overlap"] == 3000
        assert config["parameters"]["k"] == 17
        assert config["parameters"]["polish_iterations"] == 3

    def test_methylation_pipeline_config(self) -> None:
        """get_methylation_pipeline_config returns valid config."""
        config = get_methylation_pipeline_config(
            modification_types=["5mC"],
            min_coverage=10,
        )
        assert config["pipeline_name"] == "methylation"
        assert config["parameters"]["modification_types"] == ["5mC"]
        assert config["parameters"]["min_coverage"] == 10

    def test_sv_pipeline_config(self) -> None:
        """get_sv_pipeline_config returns valid config."""
        config = get_sv_pipeline_config(
            min_sv_size=100,
            min_support=5,
            sv_types=["DEL", "INS"],
        )
        assert config["pipeline_name"] == "sv"
        assert config["parameters"]["min_sv_size"] == 100
        assert config["parameters"]["sv_types"] == ["DEL", "INS"]


# ---------------------------------------------------------------------------
# validate_pipeline_config
# ---------------------------------------------------------------------------


class TestValidatePipelineConfig:
    """Tests for validate_pipeline_config function."""

    def test_valid_qc_config(self) -> None:
        """Valid QC config passes validation."""
        config = get_qc_pipeline_config()
        errors = validate_pipeline_config(config, "qc")
        assert errors == []

    def test_valid_assembly_config(self) -> None:
        """Valid assembly config passes validation."""
        config = get_assembly_pipeline_config()
        errors = validate_pipeline_config(config, "assembly")
        assert errors == []

    def test_missing_pipeline_name(self) -> None:
        """Config without pipeline_name fails validation."""
        config: dict[str, Any] = {"parameters": {}}
        errors = validate_pipeline_config(config, "qc")
        assert any("pipeline_name" in e for e in errors)

    def test_pipeline_name_mismatch(self) -> None:
        """Config with wrong pipeline name fails validation."""
        config = get_qc_pipeline_config()
        errors = validate_pipeline_config(config, "assembly")
        assert any("mismatch" in e.lower() for e in errors)

    def test_parameter_out_of_range(self) -> None:
        """Config with out-of-range parameter fails validation."""
        config = get_qc_pipeline_config()
        config["parameters"]["min_quality"] = 100.0  # Max is 60.0
        errors = validate_pipeline_config(config, "qc")
        assert any("min_quality" in e for e in errors)


# ---------------------------------------------------------------------------
# LongReadOrchestrator
# ---------------------------------------------------------------------------


class TestLongReadOrchestrator:
    """Tests for LongReadOrchestrator class."""

    def test_orchestrator_init(self, tmp_path: Path) -> None:
        """LongReadOrchestrator initializes with config and output dir."""
        config = get_qc_pipeline_config()
        orchestrator = LongReadOrchestrator(config, tmp_path / "output")
        assert orchestrator.config == config
        assert orchestrator.output_dir == tmp_path / "output"
        assert (tmp_path / "output").exists()

    def test_orchestrator_run_pipeline_invalid(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_pipeline raises ValueError for unknown pipeline."""
        config = get_qc_pipeline_config()
        orchestrator = LongReadOrchestrator(config, tmp_path / "output")
        with pytest.raises(ValueError, match="Unknown pipeline"):
            orchestrator.run_pipeline("nonexistent", [])

    def test_orchestrator_qc_pipeline_basic(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_qc_pipeline processes reads end to end."""
        config = get_qc_pipeline_config(min_length=500, min_quality=0, trim_adapters=False)

        reads = [
            {
                "read_id": f"read_{i}",
                "sequence": _random_dna(2000 + i * 500, seed=i),
                "quality_string": _phred_string([15] * (2000 + i * 500)),
            }
            for i in range(5)
        ]

        orchestrator = LongReadOrchestrator(config, tmp_path / "qc_output")
        result = orchestrator.run_qc_pipeline(reads)

        assert isinstance(result, PipelineResult)
        assert result.pipeline_name == "qc"
        assert result.total_duration > 0
        # Check that steps completed
        completed_steps = [s for s in result.steps if s.status == "completed"]
        assert len(completed_steps) >= 3  # At minimum: read_input, filter_length, filter_quality

    def test_orchestrator_qc_pipeline_with_trim(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_qc_pipeline includes adapter trimming step."""
        config = get_qc_pipeline_config(min_length=100, min_quality=0, trim_adapters=True)

        reads = [
            {
                "read_id": f"read_{i}",
                "sequence": _random_dna(3000, seed=i + 100),
                "quality_string": _phred_string([20] * 3000),
            }
            for i in range(3)
        ]

        orchestrator = LongReadOrchestrator(config, tmp_path / "qc_trim")
        result = orchestrator.run_qc_pipeline(reads)

        assert isinstance(result, PipelineResult)
        step_names = [s.name for s in result.steps]
        assert "trim_adapters" in step_names

    def test_orchestrator_qc_pipeline_empty_reads(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_qc_pipeline handles empty read list."""
        config = get_qc_pipeline_config(min_length=100, trim_adapters=False)
        orchestrator = LongReadOrchestrator(config, tmp_path / "qc_empty")
        result = orchestrator.run_qc_pipeline([])

        assert isinstance(result, PipelineResult)
        # Pipeline should still complete (with empty data)

    def test_orchestrator_resolve_dependencies(self, tmp_path: Path) -> None:
        """LongReadOrchestrator resolves step dependencies in correct order."""
        config = get_qc_pipeline_config()
        orchestrator = LongReadOrchestrator(config, tmp_path / "dep_test")

        results_order: list[str] = []

        def step_a(ctx: dict[str, Any]) -> str:
            results_order.append("a")
            return "result_a"

        def step_b(ctx: dict[str, Any]) -> str:
            results_order.append("b")
            return "result_b"

        def step_c(ctx: dict[str, Any]) -> str:
            results_order.append("c")
            return "result_c"

        steps = [
            PipelineStep(name="c", function=step_c, depends_on=["a", "b"]),
            PipelineStep(name="a", function=step_a, depends_on=[]),
            PipelineStep(name="b", function=step_b, depends_on=["a"]),
        ]

        result = orchestrator._run_steps("test", steps)

        # a must come before b, and both before c
        assert results_order.index("a") < results_order.index("b")
        assert results_order.index("b") < results_order.index("c")

    def test_orchestrator_step_failure_skips_dependents(self, tmp_path: Path) -> None:
        """Failed step causes dependent steps to be skipped."""
        config = get_qc_pipeline_config()
        orchestrator = LongReadOrchestrator(config, tmp_path / "fail_test")

        def fail_step(ctx: dict[str, Any]) -> None:
            raise RuntimeError("Intentional failure")

        def dependent_step(ctx: dict[str, Any]) -> str:
            return "should_not_run"

        steps = [
            PipelineStep(name="failing", function=fail_step, depends_on=[]),
            PipelineStep(name="dependent", function=dependent_step, depends_on=["failing"]),
        ]

        result = orchestrator._run_steps("fail_test", steps)

        assert result.success is False
        failing_step = next(s for s in result.steps if s.name == "failing")
        assert failing_step.status == "failed"

        dependent = next(s for s in result.steps if s.name == "dependent")
        assert dependent.status == "skipped"

    def test_orchestrator_sv_pipeline_basic(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_sv_pipeline processes alignments."""
        config = get_sv_pipeline_config(min_sv_size=50, min_support=2)

        # Create alignments with a deletion visible in CIGAR
        alignments = [
            {
                "read_name": f"read_{i}",
                "reference_name": "chr1",
                "reference_start": 1000,
                "cigar_string": "2000M200D3000M",
                "cigar_tuples": [],
                "mapping_quality": 60,
                "is_unmapped": False,
                "is_reverse": False,
                "tags": {},
                "query_length": 5000,
                "query_sequence": "",
            }
            for i in range(3)
        ]

        orchestrator = LongReadOrchestrator(config, tmp_path / "sv_output")
        result = orchestrator.run_sv_pipeline(alignments)

        assert isinstance(result, PipelineResult)
        assert result.pipeline_name == "sv"

    def test_orchestrator_assembly_pipeline_basic(self, tmp_path: Path) -> None:
        """LongReadOrchestrator.run_assembly_pipeline processes reads."""
        config = get_assembly_pipeline_config(
            min_read_length=100,
            min_read_quality=0,
            min_overlap=50,
        )

        seq = _random_dna(500, seed=42)
        reads = [
            {
                "read_id": f"read_{i}",
                "sequence": seq,
                "quality_string": _phred_string([20] * 500),
            }
            for i in range(3)
        ]

        orchestrator = LongReadOrchestrator(config, tmp_path / "asm_output")
        result = orchestrator.run_assembly_pipeline(reads)

        assert isinstance(result, PipelineResult)
        assert result.pipeline_name == "assembly"
        assert result.total_duration > 0
