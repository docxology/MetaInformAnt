"""Tests for phenotype workflow pipeline: PipelineConfig, PipelineResult, PhenotypePipeline.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import json

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.phenotype.workflow.pipeline import (
    PhenotypePipeline,
    PipelineConfig,
    PipelineResult,
)


# ---------------------------------------------------------------------------
# PipelineConfig
# ---------------------------------------------------------------------------


class TestPipelineConfig:
    def test_defaults(self):
        cfg = PipelineConfig()
        assert cfg.name == "phenotype_pipeline"
        assert cfg.phenotype_types == ["morphological"]
        assert cfg.steps == ["load", "validate", "analyze", "summarize"]
        assert cfg.input_path is None
        assert cfg.output_path is None
        assert cfg.parameters == {}

    def test_custom_fields(self):
        cfg = PipelineConfig(
            name="custom",
            phenotype_types=["behavioral", "sonic"],
            input_path="/tmp/in.json",
            output_path="/tmp/out.json",
            parameters={"normalize": True},
            steps=["load", "analyze"],
        )
        assert cfg.name == "custom"
        assert cfg.phenotype_types == ["behavioral", "sonic"]
        assert cfg.input_path == "/tmp/in.json"
        assert cfg.parameters == {"normalize": True}
        assert len(cfg.steps) == 2

    def test_from_dict(self):
        data = {
            "name": "dict_pipeline",
            "phenotype_types": ["chemical"],
            "steps": ["load", "preprocess", "analyze", "export"],
            "parameters": {"threshold": 0.05},
        }
        cfg = PipelineConfig.from_dict(data)
        assert cfg.name == "dict_pipeline"
        assert cfg.phenotype_types == ["chemical"]
        assert cfg.steps == ["load", "preprocess", "analyze", "export"]
        assert cfg.parameters["threshold"] == 0.05

    def test_from_dict_defaults(self):
        cfg = PipelineConfig.from_dict({})
        assert cfg.name == "phenotype_pipeline"
        assert cfg.phenotype_types == ["morphological"]

    def test_validate_valid_config(self):
        cfg = PipelineConfig(
            phenotype_types=["morphological", "behavioral"],
            steps=["load", "validate", "analyze", "summarize"],
        )
        issues = cfg.validate()
        assert issues == []

    def test_validate_unknown_phenotype_type(self):
        cfg = PipelineConfig(phenotype_types=["morphological", "unknown_type"])
        issues = cfg.validate()
        assert any("Unknown phenotype type" in i for i in issues)

    def test_validate_unknown_step(self):
        cfg = PipelineConfig(steps=["load", "magic_step"])
        issues = cfg.validate()
        assert any("Unknown step" in i for i in issues)

    def test_validate_all_valid_types(self):
        cfg = PipelineConfig(
            phenotype_types=["morphological", "behavioral", "chemical", "electronic", "sonic"]
        )
        issues = cfg.validate()
        assert issues == []

    def test_validate_all_valid_steps(self):
        cfg = PipelineConfig(
            steps=["load", "validate", "preprocess", "analyze", "summarize", "export"]
        )
        issues = cfg.validate()
        assert issues == []


# ---------------------------------------------------------------------------
# PipelineResult
# ---------------------------------------------------------------------------


class TestPipelineResult:
    def test_to_dict(self):
        cfg = PipelineConfig(name="test_pipe", phenotype_types=["morphological"])
        result = PipelineResult(
            success=True,
            config=cfg,
            outputs={"load": {"status": "ok"}, "analyze": {"status": "ok"}},
            errors=[],
            metrics={"duration_seconds": 1.5},
        )
        d = result.to_dict()
        assert d["success"] is True
        assert d["config_name"] == "test_pipe"
        assert d["phenotype_types"] == ["morphological"]
        assert "load" in d["steps_executed"]
        assert "analyze" in d["steps_executed"]
        assert d["errors"] == []

    def test_to_dict_with_errors(self):
        cfg = PipelineConfig(name="fail_pipe")
        result = PipelineResult(
            success=False,
            config=cfg,
            errors=["Step 'analyze' failed: missing data"],
        )
        d = result.to_dict()
        assert d["success"] is False
        assert len(d["errors"]) == 1

    def test_save_json(self, tmp_path):
        cfg = PipelineConfig(name="save_test")
        result = PipelineResult(
            success=True,
            config=cfg,
            outputs={"load": {"status": "ok"}},
            metrics={"steps_completed": 1},
        )
        out = tmp_path / "result.json"
        result.save_json(out)
        assert out.exists()

        with open(out) as f:
            data = json.load(f)
        assert data["config_name"] == "save_test"
        assert data["success"] is True

    def test_save_json_creates_parents(self, tmp_path):
        cfg = PipelineConfig()
        result = PipelineResult(success=True, config=cfg)
        nested = tmp_path / "a" / "b" / "result.json"
        result.save_json(nested)
        assert nested.exists()

    def test_timestamp_present(self):
        cfg = PipelineConfig()
        result = PipelineResult(success=True, config=cfg)
        assert result.timestamp is not None
        assert len(result.timestamp) > 0


# ---------------------------------------------------------------------------
# PhenotypePipeline - basic run
# ---------------------------------------------------------------------------


class TestPipelineRun:
    def test_run_with_dict_data(self):
        cfg = PipelineConfig(
            name="basic",
            steps=["load", "validate", "analyze", "summarize"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data={"key": "value"})
        assert result.success is True
        assert "load" in result.outputs
        assert "validate" in result.outputs
        assert "analyze" in result.outputs
        assert "summarize" in result.outputs

    def test_run_with_list_data(self):
        cfg = PipelineConfig(
            name="list_test",
            steps=["load", "validate", "analyze", "summarize"],
        )
        pipeline = PhenotypePipeline(cfg)
        data = [{"specimen": 1}, {"specimen": 2}, {"specimen": 3}]
        result = pipeline.run(data=data)
        assert result.success is True
        assert result.outputs["load"]["n_items"] == 3

    def test_run_no_data(self):
        cfg = PipelineConfig(
            name="empty",
            steps=["load", "validate", "summarize"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run()
        assert result.success is True
        assert result.outputs["load"]["status"] == "no_data"

    def test_run_invalid_phenotype_type(self):
        cfg = PipelineConfig(
            name="bad_type",
            phenotype_types=["invalid_type"],
            steps=["load"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data={"key": "value"})
        assert result.success is False
        assert any("Unknown phenotype type" in e for e in result.errors)

    def test_run_metrics_present(self):
        cfg = PipelineConfig(steps=["load", "summarize"])
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data={"key": "value"})
        assert "duration_seconds" in result.metrics
        assert result.metrics["steps_completed"] == 2
        assert result.metrics["steps_total"] == 2

    def test_run_from_json_file(self, tmp_path):
        # Write JSON input
        input_file = tmp_path / "input.json"
        input_data = [{"name": "specimen_1"}, {"name": "specimen_2"}]
        with open(input_file, "w") as f:
            json.dump(input_data, f)

        cfg = PipelineConfig(
            name="file_load",
            input_path=str(input_file),
            steps=["load", "validate"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run()
        assert result.success is True
        assert result.outputs["load"]["status"] == "loaded"
        assert result.outputs["load"]["n_items"] == 2


# ---------------------------------------------------------------------------
# PhenotypePipeline - step sequences
# ---------------------------------------------------------------------------


class TestPipelineSteps:
    def test_preprocess_step(self):
        cfg = PipelineConfig(steps=["load", "validate", "preprocess", "summarize"])
        pipeline = PhenotypePipeline(cfg)
        data = [{"a": 1}, {}, {"b": 2}]
        result = pipeline.run(data=data)
        assert result.success is True
        assert "preprocess" in result.outputs

    def test_export_step(self, tmp_path):
        out_file = tmp_path / "export.json"
        cfg = PipelineConfig(
            output_path=str(out_file),
            steps=["load", "analyze", "summarize", "export"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data=[{"x": 1}])
        assert result.success is True
        assert result.outputs["export"]["status"] == "exported"
        assert out_file.exists()

    def test_export_no_output_path(self):
        cfg = PipelineConfig(steps=["load", "export"])
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data={"key": "value"})
        assert result.success is True
        assert result.outputs["export"]["status"] == "skipped"

    def test_all_six_steps(self, tmp_path):
        out_file = tmp_path / "full.json"
        cfg = PipelineConfig(
            output_path=str(out_file),
            steps=["load", "validate", "preprocess", "analyze", "summarize", "export"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data=[{"x": 1}, {"y": 2}])
        assert result.success is True
        assert result.metrics["steps_completed"] == 6


# ---------------------------------------------------------------------------
# PhenotypePipeline - custom steps
# ---------------------------------------------------------------------------


class TestCustomSteps:
    def test_register_and_run_custom_step(self):
        cfg = PipelineConfig(steps=["load", "custom_norm"])
        pipeline = PhenotypePipeline(cfg)

        def custom_norm():
            return {"status": "normalized", "method": "z-score"}

        pipeline.register_step("custom_norm", custom_norm)
        result = pipeline.run(data={"key": "value"})
        assert result.success is True
        assert result.outputs["custom_norm"]["status"] == "normalized"

    def test_custom_step_overrides_builtin(self):
        cfg = PipelineConfig(steps=["load", "validate"])
        pipeline = PhenotypePipeline(cfg)

        def custom_validate():
            return {"status": "custom_validated", "all_good": True}

        pipeline.register_step("validate", custom_validate)
        result = pipeline.run(data={"key": "value"})
        assert result.outputs["validate"]["status"] == "custom_validated"


# ---------------------------------------------------------------------------
# PhenotypePipeline - morphological analyzer
# ---------------------------------------------------------------------------


class TestMorphologicalAnalyzer:
    def test_raw_dict_data(self):
        cfg = PipelineConfig(
            phenotype_types=["morphological"],
            steps=["load", "analyze"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data=[{"specimen": 1}, {"specimen": 2}])
        assert result.success is True
        analysis = result.outputs["analyze"]["results"]["morphological"]
        assert analysis["n_specimens"] == 2

    def test_single_dict_data(self):
        cfg = PipelineConfig(
            phenotype_types=["morphological"],
            steps=["load", "analyze"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data={"specimen": 1})
        assert result.success is True
        analysis = result.outputs["analyze"]["results"]["morphological"]
        assert analysis["type"] == "single_specimen"


# ---------------------------------------------------------------------------
# PhenotypePipeline - behavioral analyzer
# ---------------------------------------------------------------------------


class TestBehavioralAnalyzer:
    def test_raw_data(self):
        cfg = PipelineConfig(
            phenotype_types=["behavioral"],
            steps=["load", "analyze"],
        )
        pipeline = PhenotypePipeline(cfg)
        result = pipeline.run(data=[{"event": "A"}, {"event": "B"}])
        assert result.success is True
        analysis = result.outputs["analyze"]["results"]["behavioral"]
        assert "n_items" in analysis
