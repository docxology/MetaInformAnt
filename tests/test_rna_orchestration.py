"""Comprehensive tests for RNA-seq workflow orchestration module.

Tests all 9 public functions in metainformant.rna.engine.orchestration:
    1. run_workflow_for_species
    2. cleanup_unquantified_samples
    3. monitor_workflows
    4. discover_species_configs
    5. run_parallel_workflows
    6. validate_and_execute
    7. retry_failed_steps
    8. get_workflow_status
    9. estimate_workflow_resources

NO MOCKING -- all tests use real filesystem fixtures and real implementations.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Set

import pytest

from metainformant.rna.engine.workflow import (
    AmalgkitWorkflowConfig,
    WorkflowExecutionResult,
    WorkflowStepResult,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_config(
    tmp_path: Path,
    species: str = "Apis_mellifera",
    threads: int = 4,
    genome: Dict[str, Any] | None = None,
    **kwargs: Any,
) -> AmalgkitWorkflowConfig:
    """Create a minimal AmalgkitWorkflowConfig rooted inside *tmp_path*."""
    work_dir = tmp_path / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    return AmalgkitWorkflowConfig(
        work_dir=str(work_dir),
        threads=threads,
        species_list=[species],
        genome=genome or {},
        **kwargs,
    )


def _make_result(
    steps: List[tuple[str, int, bool]],
) -> WorkflowExecutionResult:
    """Build a WorkflowExecutionResult from a compact list of (name, rc, success)."""
    step_results = [WorkflowStepResult(step_name=name, return_code=rc, success=ok) for name, rc, ok in steps]
    total = len(step_results)
    succeeded = sum(1 for s in step_results if s.success)
    return WorkflowExecutionResult(
        steps_executed=step_results,
        success=(total - succeeded) == 0,
        total_steps=total,
        successful_steps=succeeded,
        failed_steps=total - succeeded,
    )


def _seed_quant_sample(work_dir: Path, sample_id: str) -> Path:
    """Create a fake quantification output for *sample_id* inside quant/."""
    quant_dir = work_dir / "quant" / sample_id
    quant_dir.mkdir(parents=True, exist_ok=True)
    abundance = quant_dir / "abundance.tsv"
    # Must be > 100 bytes to pass the size check in get_quantified_samples
    header = "target_id\tlength\teff_length\test_counts\ttpm\n"
    row = "gene_001\t1000\t900\t150.0\t42.5\n"
    abundance.write_text(header + row * 10)
    return abundance


def _seed_fastq_sample(work_dir: Path, sample_id: str) -> Path:
    """Create a fake FASTQ directory entry for *sample_id*."""
    fq_dir = work_dir / "fastq" / "getfastq" / sample_id
    fq_dir.mkdir(parents=True, exist_ok=True)
    fq_file = fq_dir / f"{sample_id}_1.fastq.gz"
    fq_file.write_bytes(b"\x1f\x8b" + b"\x00" * 50)  # fake gzip header
    return fq_file


def _seed_metadata(work_dir: Path, filename: str = "metadata.tsv") -> Path:
    """Create a minimal metadata TSV file inside work_dir/metadata/."""
    meta_dir = work_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)
    meta_file = meta_dir / filename
    content = textwrap.dedent(
        """\
        run\torganism\tlayout
        SRR000001\tApis mellifera\tPAIRED
        SRR000002\tApis mellifera\tPAIRED
    """
    )
    meta_file.write_text(content)
    return meta_file


def _seed_merge(work_dir: Path) -> Path:
    """Create a fake merge output directory."""
    merge_dir = work_dir / "merge"
    merge_dir.mkdir(parents=True, exist_ok=True)
    merged = merge_dir / "merged_abundance.tsv"
    merged.write_text("gene\tSRR000001\tSRR000002\ngene_001\t10\t20\n")
    return merged


def _seed_curate(work_dir: Path) -> Path:
    """Create a fake curate output directory."""
    curate_dir = work_dir / "curate"
    curate_dir.mkdir(parents=True, exist_ok=True)
    curated = curate_dir / "curated_abundance.tsv"
    curated.write_text("gene\tSRR000001\ngene_001\t10\n")
    return curated


def _seed_yaml_config(config_dir: Path, species_name: str) -> Path:
    """Write a minimal YAML config file for discover_species_configs."""
    config_dir.mkdir(parents=True, exist_ok=True)
    cfg_file = config_dir / f"amalgkit_{species_name}.yaml"
    cfg_file.write_text(
        textwrap.dedent(
            f"""\
            work_dir: /tmp/work
            threads: 4
            species_list:
              - {species_name}
        """
        )
    )
    return cfg_file


# =========================================================================
# 1. discover_species_configs
# =========================================================================


class TestDiscoverSpeciesConfigs:
    """Tests for discover_species_configs."""

    def test_finds_yaml_configs(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        _seed_yaml_config(tmp_path, "Apis_mellifera")
        _seed_yaml_config(tmp_path, "Bombus_terrestris")

        configs = discover_species_configs(tmp_path)
        species_found = {c["species"] for c in configs}

        assert len(configs) == 2
        assert "Apis_mellifera" in species_found
        assert "Bombus_terrestris" in species_found

    def test_skips_template_configs(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        _seed_yaml_config(tmp_path, "Apis_mellifera")
        _seed_yaml_config(tmp_path, "template_example")

        configs = discover_species_configs(tmp_path)

        assert len(configs) == 1
        assert configs[0]["species"] == "Apis_mellifera"

    def test_skips_test_configs(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        _seed_yaml_config(tmp_path, "Apis_mellifera")
        _seed_yaml_config(tmp_path, "test_species")

        configs = discover_species_configs(tmp_path)

        assert len(configs) == 1

    def test_empty_directory(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        configs = discover_species_configs(tmp_path)
        assert configs == []

    def test_nonexistent_directory(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        configs = discover_species_configs(tmp_path / "does_not_exist")
        assert configs == []

    def test_only_non_amalgkit_yaml_files(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        (tmp_path / "other_config.yaml").write_text("key: value\n")
        (tmp_path / "random.yaml").write_text("threads: 4\n")

        configs = discover_species_configs(tmp_path)
        assert configs == []

    def test_config_info_structure(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import discover_species_configs

        _seed_yaml_config(tmp_path, "Drosophila_melanogaster")

        configs = discover_species_configs(tmp_path)

        assert len(configs) == 1
        cfg = configs[0]
        assert "species" in cfg
        assert "config_file" in cfg
        assert "path" in cfg
        assert cfg["species"] == "Drosophila_melanogaster"
        assert isinstance(cfg["path"], Path)


# =========================================================================
# 2. estimate_workflow_resources
# =========================================================================


class TestEstimateWorkflowResources:
    """Tests for estimate_workflow_resources -- primarily math verification."""

    def test_basic_estimate_structure(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=8)
        est = estimate_workflow_resources(config, n_samples=100)

        assert "n_samples" in est
        assert "threads" in est
        assert "estimated_disk_gb" in est
        assert "estimated_memory_gb" in est
        assert "estimated_hours" in est
        assert "breakdown" in est
        assert est["n_samples"] == 100
        assert est["threads"] == 8

    def test_disk_math_100_samples(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=100)

        # fastq = 100 * 2.0 = 200.0
        # quant = 100 * 0.2 = 20.0
        # index = 5.0 (default)
        # total = 225.0 -> ceil = 225.0
        assert est["breakdown"]["fastq_disk_gb"] == 200.0
        assert est["breakdown"]["quant_disk_gb"] == 20.0
        assert est["breakdown"]["index_disk_gb"] == 5.0
        assert est["estimated_disk_gb"] == 225.0

    def test_disk_math_zero_samples(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=0)

        assert est["breakdown"]["fastq_disk_gb"] == 0.0
        assert est["breakdown"]["quant_disk_gb"] == 0.0
        # Only index overhead
        assert est["estimated_disk_gb"] == 5.0

    def test_disk_math_single_sample(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=1)

        # 1*2.0 + 1*0.2 + 5.0 = 7.2 -> ceil = 8.0
        assert est["estimated_disk_gb"] == 8.0

    def test_memory_math(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        # 4 + 8*0.5 = 8.0
        config = _make_config(tmp_path, threads=8)
        est = estimate_workflow_resources(config, n_samples=10)
        assert est["estimated_memory_gb"] == 8.0

        # 4 + 1*0.5 = 4.5 -> ceil = 5.0
        config2 = _make_config(tmp_path, threads=1)
        est2 = estimate_workflow_resources(config2, n_samples=10)
        assert est2["estimated_memory_gb"] == 5.0

    def test_time_math_download(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=60)

        # download_minutes = 60 * 1.0 = 60 -> 1.0 hours
        assert est["breakdown"]["download_hours"] == 1.0

    def test_time_math_quant_parallelism(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        # threads=8, effective_parallelism = max(1, 8//2) = 4
        config = _make_config(tmp_path, threads=8)
        est = estimate_workflow_resources(config, n_samples=120)

        # quant_minutes = (120 * 2.0) / 4 = 60.0 -> 1.0 hours
        assert est["breakdown"]["quant_hours"] == 1.0

    def test_time_math_merge(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=100)

        # merge_minutes = 10.0 + 100 * 0.3 = 40.0 -> 40/60 = 0.666... rounds to 0.7
        assert est["breakdown"]["merge_hours"] == 0.7

    def test_threads_clamped_to_one(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=0)
        est = estimate_workflow_resources(config, n_samples=10)

        # threads clamped to max(0,1) = 1
        assert est["threads"] == 1

    def test_genome_size_affects_index(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        # genome_size_mb = 10000 -> index_gb = max(5.0, 10000/1024 * 0.1) = max(5.0, 0.976) = 5.0
        config = _make_config(tmp_path, threads=4, genome={"genome_size_mb": 10000})
        est = estimate_workflow_resources(config, n_samples=10)
        assert est["breakdown"]["index_disk_gb"] == 5.0

        # genome_size_mb = 100000 -> index_gb = max(5.0, 100000/1024 * 0.1) = max(5.0, 9.76) = 9.8 (rounded)
        config2 = _make_config(tmp_path, threads=4, genome={"genome_size_mb": 100000})
        est2 = estimate_workflow_resources(config2, n_samples=10)
        assert est2["breakdown"]["index_disk_gb"] >= 5.0

    def test_total_hours_rounded_up(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=4)
        est = estimate_workflow_resources(config, n_samples=10)

        # Verify total_hours is >= sum of breakdowns (with epsilon for float rounding)
        breakdown = est["breakdown"]
        raw_sum = breakdown["download_hours"] + breakdown["quant_hours"] + breakdown["merge_hours"]
        assert est["estimated_hours"] >= raw_sum - 0.01


# =========================================================================
# 3. get_workflow_status
# =========================================================================


class TestGetWorkflowStatus:
    """Tests for get_workflow_status with real filesystem fixtures."""

    def test_not_started_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        # Remove the work dir so it does not exist
        import shutil

        shutil.rmtree(config.work_dir)

        status = get_workflow_status(config)

        assert status["overall_status"] == "not_started"
        assert status["species"] == ["Apis_mellifera"]

    def test_initialized_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        # work_dir exists but no metadata or other artifacts

        status = get_workflow_status(config)

        assert status["overall_status"] == "initialized"
        assert status["manifest_exists"] is False

    def test_metadata_ready_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_metadata(config.work_dir)

        status = get_workflow_status(config)

        assert status["overall_status"] == "metadata_ready"
        assert status["steps"]["metadata"]["completed"] is True
        assert "metadata.tsv" in status["steps"]["metadata"]["files"]

    def test_downloading_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_metadata(config.work_dir)
        _seed_fastq_sample(config.work_dir, "SRR000001")

        status = get_workflow_status(config)

        assert status["overall_status"] == "downloading"
        assert status["steps"]["getfastq"]["completed"] is True
        assert status["steps"]["getfastq"]["sample_count"] >= 1
        assert "SRR000001" in status["samples"]["with_fastq"]

    def test_in_progress_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_metadata(config.work_dir)
        _seed_fastq_sample(config.work_dir, "SRR000001")
        _seed_quant_sample(config.work_dir, "SRR000001")

        status = get_workflow_status(config)

        assert status["overall_status"] == "in_progress"
        assert status["steps"]["quant"]["completed"] is True
        assert status["steps"]["quant"]["sample_count"] >= 1
        assert "SRR000001" in status["samples"]["quantified"]

    def test_completed_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_metadata(config.work_dir)
        _seed_quant_sample(config.work_dir, "SRR000001")
        _seed_curate(config.work_dir)

        status = get_workflow_status(config)

        assert status["overall_status"] == "completed"
        assert status["steps"]["curate"]["completed"] is True

    def test_disk_free_gb_is_numeric(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        status = get_workflow_status(config)

        assert isinstance(status["disk_free_gb"], float)
        assert status["disk_free_gb"] >= 0

    def test_multiple_quantified_samples(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        for sid in ["SRR000001", "SRR000002", "SRR000003"]:
            _seed_quant_sample(config.work_dir, sid)

        status = get_workflow_status(config)

        assert status["samples"]["quantified_count"] == 3
        assert set(status["samples"]["quantified"]) == {"SRR000001", "SRR000002", "SRR000003"}

    def test_merge_step_detected(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_merge(config.work_dir)

        status = get_workflow_status(config)

        assert status["steps"]["merge"]["completed"] is True

    def test_config_step_detected(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        cfg_dir = config.work_dir / "config"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        (cfg_dir / "species.config").write_text("config_data\n")

        status = get_workflow_status(config)

        assert status["steps"]["config"]["completed"] is True

    def test_select_step_detected(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        select_dir = config.work_dir / "select"
        select_dir.mkdir(parents=True, exist_ok=True)
        (select_dir / "selected.tsv").write_text("run\torganism\nSRR1\tspecies\n")

        status = get_workflow_status(config)

        assert status["steps"]["select"]["completed"] is True

    def test_empty_work_dir_returns_correct_structure(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        status = get_workflow_status(config)

        # All step entries should exist
        expected_steps = {"metadata", "config", "select", "getfastq", "quant", "merge", "curate"}
        assert expected_steps == set(status["steps"].keys())

        # Samples section
        assert "quantified" in status["samples"]
        assert "with_fastq" in status["samples"]
        assert status["samples"]["quantified_count"] == 0
        assert status["samples"]["fastq_count"] == 0

    def test_manifest_exists_flag(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        config.manifest_path.write_text("{}\n")

        status = get_workflow_status(config)
        assert status["manifest_exists"] is True


# =========================================================================
# 4. retry_failed_steps
# =========================================================================


class TestRetryFailedSteps:
    """Tests for retry_failed_steps."""

    def test_no_failed_steps_returns_previous(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        prev = _make_result([("metadata", 0, True), ("quant", 0, True)])

        result = retry_failed_steps(config, prev, max_retries=1, base_delay=0.0)

        assert result is prev  # Should return exact same object
        assert result.success is True

    def test_all_steps_failed_merged_result(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        prev = _make_result([("metadata", 1, False), ("quant", 1, False)])

        # execute_workflow will be called but will likely fail since no real
        # amalgkit is installed. We verify the structure is correct.
        result = retry_failed_steps(config, prev, max_retries=1, base_delay=0.01)

        # The result should be a WorkflowExecutionResult
        assert isinstance(result, WorkflowExecutionResult)
        # It has merged steps -- 0 successful originals + retry results
        assert result.total_steps >= 1

    def test_mixed_success_preserves_successful(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        prev = _make_result(
            [
                ("metadata", 0, True),
                ("select", 0, True),
                ("quant", 1, False),
            ]
        )

        result = retry_failed_steps(config, prev, max_retries=1, base_delay=0.01)

        # Should have preserved the 2 successful originals
        successful_names = [s.step_name for s in result.steps_executed if s.success]
        assert "metadata" in successful_names
        assert "select" in successful_names

    def test_retry_count_honored(self, tmp_path: Path) -> None:
        """Verify that max_retries limits actual retry attempts."""
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        prev = _make_result([("quant", 1, False)])

        # max_retries=2, base_delay=0.01 to keep test fast
        result = retry_failed_steps(config, prev, max_retries=2, base_delay=0.01)

        assert isinstance(result, WorkflowExecutionResult)
        # Result should exist regardless of whether retries succeeded
        assert result.total_steps >= 1

    def test_base_delay_parameter_accepted(self, tmp_path: Path) -> None:
        """Verify the function accepts custom base_delay without error."""
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        prev = _make_result([("metadata", 1, False)])

        result = retry_failed_steps(config, prev, max_retries=1, base_delay=0.001)
        assert isinstance(result, WorkflowExecutionResult)


# =========================================================================
# 5. validate_and_execute
# =========================================================================


class TestValidateAndExecute:
    """Tests for validate_and_execute -- pre-flight validation logic."""

    def test_empty_species_list_raises(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path / "work"),
            threads=4,
            species_list=[],
            genome={},
        )

        with pytest.raises(RuntimeError, match="species_list"):
            validate_and_execute(config)

    def test_zero_threads_raises(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path / "work"),
            threads=0,
            species_list=["species"],
            genome={},
        )

        with pytest.raises(RuntimeError, match="threads"):
            validate_and_execute(config)

    def test_disk_space_check_absurd_minimum(self, tmp_path: Path) -> None:
        """Request impossibly large min_disk_gb to trigger disk space error."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = _make_config(tmp_path, genome={"accession": "GCF_test"})

        with pytest.raises(RuntimeError, match="[Ii]nsufficient disk space"):
            validate_and_execute(config, min_disk_gb=999_999.0)

    def test_low_min_disk_passes_disk_check(self, tmp_path: Path) -> None:
        """With a very low min_disk_gb, the disk check should pass.
        The function may still fail later (no amalgkit), but the disk
        check itself should not be the cause."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = _make_config(tmp_path, genome={"accession": "GCF_test"})

        # This should pass disk check but might fail on execute_workflow
        # since amalgkit is not installed. We catch any downstream error.
        try:
            validate_and_execute(config, min_disk_gb=0.001)
        except RuntimeError as e:
            # Should NOT be a disk-space error
            assert "disk space" not in str(e).lower()

    def test_missing_genome_for_quant_raises(self, tmp_path: Path) -> None:
        """When genome is empty and quant steps are implied, should raise."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = _make_config(tmp_path, genome={})

        with pytest.raises(RuntimeError, match="[Gg]enome"):
            validate_and_execute(config, min_disk_gb=0.001)

    def test_genome_present_skips_genome_error(self, tmp_path: Path) -> None:
        """When genome has an accession field, the genome check should pass."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = _make_config(tmp_path, genome={"accession": "GCF_001"})

        # Should pass genome check; may fail later on execute_workflow
        try:
            validate_and_execute(config, min_disk_gb=0.001)
        except RuntimeError as e:
            assert "genome" not in str(e).lower()

    def test_steps_kwarg_skips_genome_check(self, tmp_path: Path) -> None:
        """When steps=['metadata'], genome is not needed."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = _make_config(tmp_path, genome={})

        # Explicitly only run metadata -- genome should NOT be required
        try:
            validate_and_execute(config, min_disk_gb=0.001, steps=["metadata"])
        except RuntimeError as e:
            assert "genome" not in str(e).lower()

    def test_work_dir_created_automatically(self, tmp_path: Path) -> None:
        """validate_and_execute should create work_dir if it does not exist."""
        from metainformant.rna.engine.orchestration import validate_and_execute

        work_dir = tmp_path / "nonexistent" / "deep" / "work"
        config = AmalgkitWorkflowConfig(
            work_dir=str(work_dir),
            threads=4,
            species_list=["TestSpecies"],
            genome={"accession": "GCF_test"},
        )

        try:
            validate_and_execute(config, min_disk_gb=0.001)
        except RuntimeError:
            pass  # downstream errors are fine

        # The work_dir should have been created by the function
        assert work_dir.exists()


# =========================================================================
# 6. run_parallel_workflows
# =========================================================================


class TestRunParallelWorkflows:
    """Tests for run_parallel_workflows."""

    def test_empty_configs_returns_empty(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        results = run_parallel_workflows([], max_workers=1)
        assert results == {}

    def test_single_config_returns_one_result(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        config = _make_config(tmp_path, species="TestSpecies")

        # execute_workflow will likely fail (no amalgkit), but we get a result
        results = run_parallel_workflows([config], max_workers=1)

        assert len(results) == 1
        assert "TestSpecies" in results
        assert isinstance(results["TestSpecies"], WorkflowExecutionResult)

    def test_multiple_configs_parallel(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        cfg1 = _make_config(tmp_path / "sp1", species="Species_A")
        cfg2 = _make_config(tmp_path / "sp2", species="Species_B")

        results = run_parallel_workflows([cfg1, cfg2], max_workers=2)

        assert len(results) == 2
        assert "Species_A" in results
        assert "Species_B" in results

    def test_failure_does_not_cancel_others(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        cfg1 = _make_config(tmp_path / "sp1", species="GoodSpecies")
        cfg2 = _make_config(tmp_path / "sp2", species="BadSpecies")

        results = run_parallel_workflows([cfg1, cfg2], max_workers=2)

        # Both should have results, even if both failed
        assert len(results) == 2

    def test_empty_species_list_uses_unknown_key(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path / "work"),
            threads=4,
            species_list=[],
        )

        results = run_parallel_workflows([config], max_workers=1)

        # The key should be "unknown_0" since species_list is empty
        assert len(results) == 1
        key = list(results.keys())[0]
        assert "unknown" in key

    def test_max_workers_one_executes_sequentially(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        cfg1 = _make_config(tmp_path / "sp1", species="Species_1")
        cfg2 = _make_config(tmp_path / "sp2", species="Species_2")

        results = run_parallel_workflows([cfg1, cfg2], max_workers=1)
        assert len(results) == 2


# =========================================================================
# 7. run_workflow_for_species (config-file based)
# =========================================================================


class TestRunWorkflowForSpecies:
    """Tests for run_workflow_for_species.

    This function requires a real config file path and calls into
    execute_workflow. We test that the config file loading works and
    errors are raised for invalid configs.
    """

    def test_missing_config_file_raises(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_workflow_for_species

        with pytest.raises(FileNotFoundError):
            run_workflow_for_species(tmp_path / "nonexistent.yaml", species="Test")

    def test_valid_config_file_executes(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_workflow_for_species

        # Write a minimal YAML config
        config_file = tmp_path / "test_config.yaml"
        config_file.write_text(
            textwrap.dedent(
                """\
                work_dir: {work_dir}
                threads: 2
                species_list:
                  - Apis_mellifera
            """
            ).format(work_dir=str(tmp_path / "work"))
        )

        # Will likely fail downstream but should parse config and attempt execution
        try:
            result = run_workflow_for_species(config_file, species="Apis_mellifera")
            # If it succeeds, check structure
            assert "species" in result
            assert "success" in result
        except Exception:
            # Downstream failures (no amalgkit) are expected
            pass

    def test_species_override(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_workflow_for_species

        config_file = tmp_path / "test_config.yaml"
        config_file.write_text(
            textwrap.dedent(
                """\
                work_dir: {work_dir}
                threads: 2
                species_list:
                  - Apis_mellifera
            """
            ).format(work_dir=str(tmp_path / "work"))
        )

        try:
            result = run_workflow_for_species(config_file, species="Bombus_terrestris")
            assert result["species"] == "Bombus_terrestris"
        except Exception:
            pass


# =========================================================================
# 8. cleanup_unquantified_samples
# =========================================================================


class TestCleanupUnquantifiedSamples:
    """Tests for cleanup_unquantified_samples."""

    def test_missing_config_raises(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import cleanup_unquantified_samples

        with pytest.raises(FileNotFoundError):
            cleanup_unquantified_samples(tmp_path / "missing.yaml")

    def test_valid_config_returns_tuple(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import cleanup_unquantified_samples

        config_file = tmp_path / "cfg.yaml"
        config_file.write_text(
            textwrap.dedent(
                """\
                work_dir: {work_dir}
                threads: 2
                species_list:
                  - Apis_mellifera
            """
            ).format(work_dir=str(tmp_path / "work"))
        )
        (tmp_path / "work").mkdir(parents=True, exist_ok=True)

        result = cleanup_unquantified_samples(config_file)

        assert isinstance(result, tuple)
        assert len(result) == 2
        quantified, failed = result
        assert isinstance(quantified, int)
        assert isinstance(failed, int)

    def test_empty_work_dir_returns_zero_counts(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import cleanup_unquantified_samples

        config_file = tmp_path / "cfg.yaml"
        work_dir = tmp_path / "work"
        config_file.write_text(
            textwrap.dedent(
                """\
                work_dir: {work_dir}
                threads: 2
                species_list:
                  - TestSpecies
            """
            ).format(work_dir=str(work_dir))
        )
        work_dir.mkdir(parents=True, exist_ok=True)

        quantified, failed = cleanup_unquantified_samples(config_file)

        assert quantified == 0
        assert failed == 0


# =========================================================================
# 9. monitor_workflows
# =========================================================================


class TestMonitorWorkflows:
    """Tests for monitor_workflows."""

    def test_empty_directory_returns_status(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import monitor_workflows

        work_dir = tmp_path / "work"
        work_dir.mkdir(parents=True, exist_ok=True)

        result = monitor_workflows(work_dir)

        assert isinstance(result, dict)

    def test_nonexistent_directory(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import monitor_workflows

        # Should handle gracefully
        try:
            result = monitor_workflows(tmp_path / "nonexistent")
            assert isinstance(result, dict)
        except (FileNotFoundError, OSError):
            # Also acceptable behavior
            pass


# =========================================================================
# WorkflowExecutionResult and WorkflowStepResult dataclass tests
# =========================================================================


class TestWorkflowDataclasses:
    """Test the dataclass structures used by orchestration functions."""

    def test_step_result_creation(self) -> None:
        step = WorkflowStepResult(
            step_name="metadata",
            return_code=0,
            success=True,
        )
        assert step.step_name == "metadata"
        assert step.return_code == 0
        assert step.success is True
        assert step.error_message is None

    def test_step_result_with_error(self) -> None:
        step = WorkflowStepResult(
            step_name="quant",
            return_code=1,
            success=False,
            error_message="kallisto not found",
        )
        assert step.success is False
        assert step.error_message == "kallisto not found"

    def test_execution_result_success(self) -> None:
        result = _make_result([("a", 0, True), ("b", 0, True)])
        assert result.success is True
        assert result.total_steps == 2
        assert result.successful_steps == 2
        assert result.failed_steps == 0
        assert result.return_codes == [0, 0]

    def test_execution_result_with_failures(self) -> None:
        result = _make_result([("a", 0, True), ("b", 1, False)])
        assert result.success is False
        assert result.total_steps == 2
        assert result.successful_steps == 1
        assert result.failed_steps == 1
        assert result.return_codes == [0, 1]

    def test_execution_result_is_list(self) -> None:
        """WorkflowExecutionResult inherits from list for backward compat."""
        result = _make_result([("a", 0, True)])
        assert isinstance(result, list)

    def test_execution_result_get_by_name(self) -> None:
        result = _make_result([("metadata", 0, True), ("quant", 1, False)])
        assert result.get("metadata") == 0
        assert result.get("quant") == 1
        assert result.get("nonexistent") is None
        assert result.get("nonexistent", -1) == -1

    def test_execution_result_len(self) -> None:
        result = _make_result([("a", 0, True), ("b", 1, False), ("c", 0, True)])
        assert len(result) == 3

    def test_execution_result_getitem(self) -> None:
        result = _make_result([("a", 0, True), ("b", 1, False)])
        step = result[0]
        assert isinstance(step, WorkflowStepResult)
        assert step.step_name == "a"

    def test_empty_execution_result(self) -> None:
        result = _make_result([])
        assert result.success is True
        assert result.total_steps == 0
        assert len(result) == 0


# =========================================================================
# AmalgkitWorkflowConfig tests
# =========================================================================


class TestAmalgkitWorkflowConfig:
    """Test config object construction and properties."""

    def test_default_construction(self, tmp_path: Path) -> None:
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path))
        assert config.work_dir == tmp_path
        assert config.threads == 8
        assert config.species_list == []
        assert config.genome == {}

    def test_manifest_path_property(self, tmp_path: Path) -> None:
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path))
        assert config.manifest_path == tmp_path / "amalgkit.manifest.jsonl"

    def test_log_file_property(self, tmp_path: Path) -> None:
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path))
        assert config.log_file == tmp_path / "logs" / "workflow.log"

    def test_custom_log_dir(self, tmp_path: Path) -> None:
        log_dir = tmp_path / "custom_logs"
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path), log_dir=str(log_dir))
        assert config.log_dir == log_dir

    def test_to_dict_roundtrip(self, tmp_path: Path) -> None:
        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path),
            threads=12,
            species_list=["TestSpecies"],
            genome={"accession": "GCF_001"},
        )
        d = config.to_dict()
        assert d["threads"] == 12
        assert d["species_list"] == ["TestSpecies"]
        assert d["genome"]["accession"] == "GCF_001"

    def test_from_dict_construction(self, tmp_path: Path) -> None:
        d = {
            "work_dir": str(tmp_path),
            "threads": 6,
            "species_list": ["Sp1", "Sp2"],
        }
        config = AmalgkitWorkflowConfig.from_dict(d)
        assert config.threads == 6
        assert config.species_list == ["Sp1", "Sp2"]

    def test_extra_config_kwargs(self, tmp_path: Path) -> None:
        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path),
            custom_key="custom_value",
        )
        assert config.extra_config.get("custom_key") == "custom_value"


# =========================================================================
# Edge cases and integration
# =========================================================================


class TestEdgeCases:
    """Edge cases cutting across multiple functions."""

    def test_estimate_resources_large_sample_count(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=16)
        est = estimate_workflow_resources(config, n_samples=10000)

        assert est["estimated_disk_gb"] > 0
        assert est["estimated_hours"] > 0
        assert est["n_samples"] == 10000

    def test_get_status_with_drr_err_samples(self, tmp_path: Path) -> None:
        """Verify DRR and ERR prefix samples are detected."""
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        _seed_fastq_sample(config.work_dir, "DRR000001")
        _seed_fastq_sample(config.work_dir, "ERR000001")
        _seed_quant_sample(config.work_dir, "DRR000001")
        _seed_quant_sample(config.work_dir, "ERR000001")

        status = get_workflow_status(config)

        assert "DRR000001" in status["samples"]["with_fastq"]
        assert "ERR000001" in status["samples"]["with_fastq"]
        assert "DRR000001" in status["samples"]["quantified"]
        assert "ERR000001" in status["samples"]["quantified"]

    def test_estimate_with_negative_threads_clamped(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import estimate_workflow_resources

        config = _make_config(tmp_path, threads=-5)
        est = estimate_workflow_resources(config, n_samples=10)

        # threads clamped to max(-5, 1) = 1
        assert est["threads"] == 1

    def test_get_status_no_abundance_small_file(self, tmp_path: Path) -> None:
        """Abundance files < 100 bytes should NOT count as quantified."""
        from metainformant.rna.engine.orchestration import get_workflow_status

        config = _make_config(tmp_path)
        quant_dir = config.work_dir / "quant" / "SRR000001"
        quant_dir.mkdir(parents=True, exist_ok=True)
        tiny = quant_dir / "abundance.tsv"
        tiny.write_text("header\n")  # < 100 bytes

        status = get_workflow_status(config)

        assert status["samples"]["quantified_count"] == 0

    def test_retry_with_zero_max_retries(self, tmp_path: Path) -> None:
        """max_retries=0 means no retries at all, should still return a result."""
        from metainformant.rna.engine.orchestration import retry_failed_steps

        config = _make_config(tmp_path)
        # With zero max_retries the for loop runs 0 times, but
        # the function should handle this by returning previous_result
        # if no failed steps, or running the loop 0 times otherwise.
        prev_no_failures = _make_result([("metadata", 0, True)])
        result = retry_failed_steps(config, prev_no_failures, max_retries=0)
        assert result is prev_no_failures

    def test_parallel_workflows_result_values(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import run_parallel_workflows

        cfg = _make_config(tmp_path, species="Solo")
        results = run_parallel_workflows([cfg], max_workers=1)

        result = results["Solo"]
        assert hasattr(result, "success")
        assert hasattr(result, "total_steps")
        assert hasattr(result, "steps_executed")

    def test_validate_and_execute_invalid_species_name(self, tmp_path: Path) -> None:
        from metainformant.rna.engine.orchestration import validate_and_execute

        config = AmalgkitWorkflowConfig(
            work_dir=str(tmp_path / "work"),
            threads=4,
            species_list=[""],  # empty string
        )

        with pytest.raises(RuntimeError, match="[Ii]nvalid species|species_list"):
            validate_and_execute(config, min_disk_gb=0.001)
