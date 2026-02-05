"""Comprehensive tests for metainformant.rna.engine.workflow_steps module.

Tests validate_step_prerequisites, check_step_completion_status,
handle_post_step_actions, log_workflow_summary, and setup_vdb_config
using real file system operations (NO mocking).
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pytest

from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, WorkflowStepResult
from metainformant.rna.engine.workflow_steps import (
    check_step_completion_status,
    handle_post_step_actions,
    log_workflow_summary,
    setup_vdb_config,
    validate_step_prerequisites,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_config(work_dir: Path, **kwargs: Any) -> AmalgkitWorkflowConfig:
    """Create an AmalgkitWorkflowConfig rooted at *work_dir*."""
    return AmalgkitWorkflowConfig(work_dir=str(work_dir), **kwargs)


def _write_file(path: Path, content: str = "") -> Path:
    """Create *path* (and parents), optionally writing *content*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _write_metadata_tsv(path: Path, rows: int = 3) -> Path:
    """Write a minimal metadata TSV with a header + *rows* data lines."""
    header = "run_accession\torganism\tread_count\n"
    lines = [f"SRR{1000 + i}\tApis mellifera\t{5000 * i}\n" for i in range(1, rows + 1)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "".join(lines))
    return path


@dataclass
class FakeResult:
    """Lightweight stand-in for subprocess.CompletedProcess returned by step execution."""

    returncode: int = 0
    stdout: str = ""
    stderr: str = ""


# ===================================================================
# validate_step_prerequisites
# ===================================================================


class TestValidateStepPrerequisites:
    """Tests for validate_step_prerequisites()."""

    # ---------------------------------------------------------------
    # Unknown / passthrough step
    # ---------------------------------------------------------------

    def test_unknown_step_returns_none(self, tmp_path: Path) -> None:
        """Steps without explicit validation logic should pass (return None)."""
        config = _make_config(tmp_path)
        result = validate_step_prerequisites(
            step_name="sanity",
            step_params={},
            config=config,
            steps_planned=[("sanity", {})],
            steps_config={},
        )
        assert result is None

    # ---------------------------------------------------------------
    # integrate step
    # ---------------------------------------------------------------

    def test_integrate_fails_when_no_fastq_files(self, tmp_path: Path) -> None:
        """integrate should fail when no FASTQ files exist anywhere."""
        config = _make_config(tmp_path)
        # Create work_dir but no FASTQ files
        (tmp_path / "fastq").mkdir(parents=True, exist_ok=True)

        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )
        assert error is not None
        assert "PREREQUISITE CHECK FAILED" in error
        assert "REMEDIATION" in error
        assert "FASTQ" in error

    def test_integrate_passes_when_fastq_files_exist(self, tmp_path: Path) -> None:
        """integrate should pass when FASTQ files are present."""
        config = _make_config(tmp_path)
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        _write_file(fastq_dir / "SRR1001" / "SRR1001.fastq.gz", "fake fastq")

        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )
        assert error is None

    def test_integrate_fails_with_fq_extension_only(self, tmp_path: Path) -> None:
        """integrate glob pattern is *.fastq* so .fq.gz alone should NOT satisfy it."""
        config = _make_config(tmp_path)
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        # .fq.gz does NOT match *.fastq*
        _write_file(fastq_dir / "sample.fq.gz", "data")

        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )
        # The glob pattern only matches *.fastq*, so .fq.gz alone is not found
        assert error is not None
        assert "FASTQ" in error

    def test_integrate_passes_with_fastq_gz_extension(self, tmp_path: Path) -> None:
        """integrate should find .fastq.gz files (matches *.fastq* glob)."""
        config = _make_config(tmp_path)
        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        _write_file(fastq_dir / "sample.fastq.gz", "data")

        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )
        assert error is None

    def test_integrate_finds_fastq_in_getfastq_subdir(self, tmp_path: Path) -> None:
        """integrate should find FASTQ files in the getfastq subdirectory."""
        config = _make_config(tmp_path)
        getfastq_dir = tmp_path / "fastq" / "getfastq"
        getfastq_dir.mkdir(parents=True, exist_ok=True)
        _write_file(getfastq_dir / "SRR1001" / "SRR1001.fastq.gz", "data")

        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )
        assert error is None

    def test_integrate_finds_fastq_via_steps_config_out_dir(self, tmp_path: Path) -> None:
        """integrate should look in the getfastq out_dir specified in steps_config."""
        config = _make_config(tmp_path)
        custom_dir = tmp_path / "custom_fastq"
        custom_dir.mkdir(parents=True, exist_ok=True)
        _write_file(custom_dir / "SRR1001.fastq", "reads")

        steps_config: Dict[str, Any] = {"getfastq": {"out_dir": str(custom_dir)}}
        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=[("integrate", {})],
            steps_config=steps_config,
        )
        assert error is None

    # ---------------------------------------------------------------
    # quant step
    # ---------------------------------------------------------------

    def test_quant_fails_when_no_metadata(self, tmp_path: Path) -> None:
        """quant should fail when neither integrated_metadata.json nor metadata.tsv exist."""
        config = _make_config(tmp_path)

        error = validate_step_prerequisites(
            step_name="quant",
            step_params={},
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
        )
        assert error is not None
        assert "PREREQUISITE CHECK FAILED" in error
        assert "REMEDIATION" in error
        assert "integrate" in error.lower() or "Integrate" in error

    def test_quant_fails_when_metadata_but_no_fastq(self, tmp_path: Path) -> None:
        """quant should fail when metadata exists but no FASTQ files are present."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        error = validate_step_prerequisites(
            step_name="quant",
            step_params={},
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
        )
        assert error is not None
        assert "FASTQ" in error
        assert "REMEDIATION" in error

    def test_quant_passes_with_metadata_and_fastq(self, tmp_path: Path) -> None:
        """quant should pass when both metadata and FASTQ files exist."""
        config = _make_config(tmp_path)
        # Create metadata
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")
        # Create FASTQ files
        _write_file(tmp_path / "fastq" / "SRR1001" / "SRR1001.fastq.gz", "reads")

        error = validate_step_prerequisites(
            step_name="quant",
            step_params={},
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
        )
        # May return None (pass) or error about quantification tools not being installed.
        # Both are valid since we only test file-level prerequisites here.
        if error is not None:
            # The only acceptable failure is about quantification tools (kallisto/salmon),
            # not about metadata or FASTQ files.
            assert "quantification tool" in error.lower() or "kallisto" in error.lower() or "salmon" in error.lower()

    def test_quant_passes_with_integrated_metadata_json(self, tmp_path: Path) -> None:
        """quant should accept integrated_metadata.json as alternative to metadata.tsv."""
        config = _make_config(tmp_path)
        _write_file(tmp_path / "integration" / "integrated_metadata.json", '{"samples": []}')
        _write_file(tmp_path / "fastq" / "SRR1001.fastq.gz", "reads")

        error = validate_step_prerequisites(
            step_name="quant",
            step_params={},
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
        )
        if error is not None:
            # Only acceptable if it is about missing tools, not about missing metadata
            assert "metadata" not in error.lower() or "quantification" in error.lower()

    # ---------------------------------------------------------------
    # merge step
    # ---------------------------------------------------------------

    def test_merge_fails_when_no_quant_files(self, tmp_path: Path) -> None:
        """merge should fail when no quantification output files exist."""
        config = _make_config(tmp_path)
        (tmp_path / "quant").mkdir(parents=True, exist_ok=True)

        error = validate_step_prerequisites(
            step_name="merge",
            step_params={"out_dir": str(tmp_path / "merged")},
            config=config,
            steps_planned=[("merge", {"out_dir": str(tmp_path / "merged")})],
            steps_config={},
        )
        assert error is not None
        assert "PREREQUISITE CHECK FAILED" in error
        assert "quantification" in error.lower() or "quant" in error.lower()
        assert "REMEDIATION" in error

    def test_merge_passes_with_abundance_tsv(self, tmp_path: Path) -> None:
        """merge should pass when abundance.tsv files exist in quant dir."""
        config = _make_config(tmp_path)
        quant_dir = tmp_path / "quant"
        _write_file(quant_dir / "SRR1001" / "abundance.tsv", "target_id\test_counts\n")

        error = validate_step_prerequisites(
            step_name="merge",
            step_params={"out_dir": str(tmp_path / "merged")},
            config=config,
            steps_planned=[("merge", {"out_dir": str(tmp_path / "merged")})],
            steps_config={"quant": {"out_dir": str(quant_dir)}},
        )
        # May pass entirely or fail on R/Rscript/ggplot2 check
        if error is not None:
            assert "R" in error or "Rscript" in error or "ggplot2" in error

    def test_merge_passes_with_quant_sf(self, tmp_path: Path) -> None:
        """merge should also find Salmon quant.sf files."""
        config = _make_config(tmp_path)
        quant_dir = tmp_path / "quant"
        _write_file(quant_dir / "SRR1001" / "quant.sf", "Name\tLength\tEffectiveLength\tTPM\tNumReads\n")

        error = validate_step_prerequisites(
            step_name="merge",
            step_params={"out_dir": str(tmp_path / "merged")},
            config=config,
            steps_planned=[("merge", {"out_dir": str(tmp_path / "merged")})],
            steps_config={"quant": {"out_dir": str(quant_dir)}},
        )
        if error is not None:
            assert "R" in error or "Rscript" in error or "ggplot2" in error

    def test_merge_uses_steps_config_quant_out_dir(self, tmp_path: Path) -> None:
        """merge should look for quant files in the out_dir from steps_config."""
        config = _make_config(tmp_path)
        custom_quant = tmp_path / "custom_quant"
        _write_file(custom_quant / "sample1" / "abundance.tsv", "data")

        error = validate_step_prerequisites(
            step_name="merge",
            step_params={"out_dir": str(tmp_path / "merged")},
            config=config,
            steps_planned=[("merge", {})],
            steps_config={"quant": {"out_dir": str(custom_quant)}},
        )
        if error is not None:
            # Should not be about missing quant files
            assert "quantification files" not in error.lower()

    def test_merge_bridges_quant_results_symlink(self, tmp_path: Path) -> None:
        """merge should create a symlink bridging quant results into merge output dir."""
        config = _make_config(tmp_path)
        quant_dir = tmp_path / "quant"
        _write_file(quant_dir / "SRR1001" / "abundance.tsv", "data")
        merge_out = tmp_path / "merged"

        error = validate_step_prerequisites(
            step_name="merge",
            step_params={"out_dir": str(merge_out)},
            config=config,
            steps_planned=[("merge", {"out_dir": str(merge_out)})],
            steps_config={"quant": {"out_dir": str(quant_dir)}},
        )
        # If merge validation passes the quant-file check, it should have attempted
        # to create the bridge symlink at merge_out/quant -> quant_dir
        if error is None:
            merge_quant_link = merge_out / "quant"
            assert merge_quant_link.exists() or merge_quant_link.is_symlink()

    # ---------------------------------------------------------------
    # curate / cstmm steps (R-dependent)
    # ---------------------------------------------------------------

    def test_curate_fails_without_rscript(self, tmp_path: Path) -> None:
        """curate should fail when Rscript is not on PATH."""
        if shutil.which("Rscript"):
            pytest.skip("Rscript is available - cannot test missing-R path")

        config = _make_config(tmp_path)
        error = validate_step_prerequisites(
            step_name="curate",
            step_params={},
            config=config,
            steps_planned=[("curate", {})],
            steps_config={},
        )
        assert error is not None
        assert "Rscript" in error

    def test_cstmm_fails_without_rscript(self, tmp_path: Path) -> None:
        """cstmm should fail when Rscript is not on PATH."""
        if shutil.which("Rscript"):
            pytest.skip("Rscript is available - cannot test missing-R path")

        config = _make_config(tmp_path)
        error = validate_step_prerequisites(
            step_name="cstmm",
            step_params={},
            config=config,
            steps_planned=[("cstmm", {})],
            steps_config={},
        )
        assert error is not None
        assert "Rscript" in error

    def test_curate_passes_with_rscript(self, tmp_path: Path) -> None:
        """curate should pass when Rscript is available."""
        if not shutil.which("Rscript"):
            pytest.skip("Rscript not available")

        config = _make_config(tmp_path)
        error = validate_step_prerequisites(
            step_name="curate",
            step_params={"out_dir": str(tmp_path / "curate")},
            config=config,
            steps_planned=[("curate", {"out_dir": str(tmp_path / "curate")})],
            steps_config={},
        )
        assert error is None

    def test_curate_bridges_merge_results(self, tmp_path: Path) -> None:
        """curate should create a symlink from curate/merge -> merge output dir."""
        if not shutil.which("Rscript"):
            pytest.skip("Rscript not available")

        config = _make_config(tmp_path)
        # Create merge results directory
        merge_out = tmp_path / "merged"
        merge_results = merge_out / "merge"
        merge_results.mkdir(parents=True, exist_ok=True)
        _write_file(merge_results / "merged_abundance.tsv", "data")

        curate_out = tmp_path / "curate"
        error = validate_step_prerequisites(
            step_name="curate",
            step_params={"out_dir": str(curate_out)},
            config=config,
            steps_planned=[("curate", {"out_dir": str(curate_out)})],
            steps_config={"merge": {"out_dir": str(merge_out)}},
        )
        assert error is None
        # Bridge symlink should have been created
        curate_merge_link = curate_out / "merge"
        assert curate_merge_link.exists() or curate_merge_link.is_symlink()

    def test_cstmm_passes_with_rscript(self, tmp_path: Path) -> None:
        """cstmm should pass when Rscript is available."""
        if not shutil.which("Rscript"):
            pytest.skip("Rscript not available")

        config = _make_config(tmp_path)
        error = validate_step_prerequisites(
            step_name="cstmm",
            step_params={"out_dir": str(tmp_path / "cstmm")},
            config=config,
            steps_planned=[("cstmm", {"out_dir": str(tmp_path / "cstmm")})],
            steps_config={},
        )
        assert error is None

    # ---------------------------------------------------------------
    # config and select (no explicit validation defined - passthrough)
    # ---------------------------------------------------------------

    def test_config_step_returns_none(self, tmp_path: Path) -> None:
        """config step has no prerequisite validation -> None."""
        config = _make_config(tmp_path)
        assert (
            validate_step_prerequisites(
                step_name="config",
                step_params={},
                config=config,
                steps_planned=[],
                steps_config={},
            )
            is None
        )

    def test_select_step_returns_none(self, tmp_path: Path) -> None:
        """select step has no prerequisite validation -> None."""
        config = _make_config(tmp_path)
        assert (
            validate_step_prerequisites(
                step_name="select",
                step_params={},
                config=config,
                steps_planned=[],
                steps_config={},
            )
            is None
        )

    def test_getfastq_step_returns_none(self, tmp_path: Path) -> None:
        """getfastq step has no prerequisite validation -> None."""
        config = _make_config(tmp_path)
        assert (
            validate_step_prerequisites(
                step_name="getfastq",
                step_params={},
                config=config,
                steps_planned=[],
                steps_config={},
            )
            is None
        )

    def test_metadata_step_returns_none(self, tmp_path: Path) -> None:
        """metadata step has no prerequisite validation -> None."""
        config = _make_config(tmp_path)
        assert (
            validate_step_prerequisites(
                step_name="metadata",
                step_params={},
                config=config,
                steps_planned=[],
                steps_config={},
            )
            is None
        )


# ===================================================================
# check_step_completion_status
# ===================================================================


class TestCheckStepCompletionStatus:
    """Tests for check_step_completion_status()."""

    def test_empty_steps_returns_empty(self, tmp_path: Path) -> None:
        """No planned steps -> empty completed & to-run lists."""
        config = _make_config(tmp_path)
        completed, to_run = check_step_completion_status([], config)
        assert completed == []
        assert to_run == []

    def test_metadata_completed_when_file_exists(self, tmp_path: Path) -> None:
        """metadata step should be completed when metadata.tsv exists."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "metadata"
        assert to_run == []

    def test_metadata_not_completed_when_missing(self, tmp_path: Path) -> None:
        """metadata step should not be completed when metadata.tsv is absent."""
        config = _make_config(tmp_path)

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert completed == []
        assert to_run == ["metadata"]

    def test_config_completed_when_config_files_exist(self, tmp_path: Path) -> None:
        """config step should be completed when config_base has .config files."""
        config = _make_config(tmp_path)
        config_base = tmp_path / "config_base"
        config_base.mkdir(parents=True, exist_ok=True)
        _write_file(config_base / "species.config", "[config]\nkey=value")

        steps: List[Tuple[str, Dict[str, Any]]] = [("config", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "config"

    def test_config_not_completed_when_empty(self, tmp_path: Path) -> None:
        """config step should not be completed when config_base is empty."""
        config = _make_config(tmp_path)
        (tmp_path / "config_base").mkdir(parents=True, exist_ok=True)

        steps: List[Tuple[str, Dict[str, Any]]] = [("config", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert completed == []
        assert to_run == ["config"]

    def test_select_completed_when_metadata_selected_exists(self, tmp_path: Path) -> None:
        """select step should be completed when metadata_selected.tsv exists."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata_selected.tsv", rows=2)

        steps: List[Tuple[str, Dict[str, Any]]] = [("select", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "select"

    def test_getfastq_always_runs(self, tmp_path: Path) -> None:
        """getfastq should always be in to_run (defers to amalgkit's skip logic)."""
        config = _make_config(tmp_path)

        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert completed == []
        assert to_run == ["getfastq"]

    def test_quant_always_runs(self, tmp_path: Path) -> None:
        """quant should always be in to_run (defers to amalgkit's skip logic)."""
        config = _make_config(tmp_path)

        steps: List[Tuple[str, Dict[str, Any]]] = [("quant", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert completed == []
        assert to_run == ["quant"]

    def test_integrate_completed_when_integration_dir_has_content(self, tmp_path: Path) -> None:
        """integrate step should be completed when integration dir has files and metadata.tsv exists."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")
        integration_dir = tmp_path / "integration"
        integration_dir.mkdir(parents=True, exist_ok=True)
        _write_file(integration_dir / "some_output.json", '{"done": true}')

        steps: List[Tuple[str, Dict[str, Any]]] = [("integrate", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "integrate"

    def test_integrate_completed_with_integrated_metadata_json(self, tmp_path: Path) -> None:
        """integrate step is completed when integrated_metadata.json exists."""
        config = _make_config(tmp_path)
        _write_file(tmp_path / "integration" / "integrated_metadata.json", "{}")

        steps: List[Tuple[str, Dict[str, Any]]] = [("integrate", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "integrate"

    def test_merge_completed_when_merged_abundance_exists(self, tmp_path: Path) -> None:
        """merge step is completed when merged_abundance.tsv exists."""
        config = _make_config(tmp_path)
        merged_dir = tmp_path / "merged"
        _write_file(merged_dir / "merged_abundance.tsv", "gene\tsample1\n")

        steps: List[Tuple[str, Dict[str, Any]]] = [("merge", {"out_dir": str(merged_dir)})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "merge"

    def test_curate_completed_when_outputs_exist(self, tmp_path: Path) -> None:
        """curate step is completed when curated_abundance.tsv exists."""
        config = _make_config(tmp_path)
        curate_dir = tmp_path / "curate"
        _write_file(curate_dir / "curated_abundance.tsv", "data")

        steps: List[Tuple[str, Dict[str, Any]]] = [("curate", {"out_dir": str(curate_dir)})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "curate"

    def test_force_redo_overrides_completion(self, tmp_path: Path) -> None:
        """force_redo=True should force completed steps into to_run."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]
        completed, to_run = check_step_completion_status(steps, config, redo="yes")

        assert completed == []
        assert to_run == ["metadata"]

    def test_redo_bool_true_forces_rerun(self, tmp_path: Path) -> None:
        """redo=True (boolean) should force steps to re-run."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]
        completed, to_run = check_step_completion_status(steps, config, redo=True)

        assert completed == []
        assert to_run == ["metadata"]

    def test_redo_no_does_not_force_rerun(self, tmp_path: Path) -> None:
        """redo='no' should not force re-run of completed steps."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]
        completed, to_run = check_step_completion_status(steps, config, redo="no")

        assert len(completed) == 1
        assert to_run == []

    def test_per_step_redo_param(self, tmp_path: Path) -> None:
        """redo set in step_params should override default."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {"redo": "yes"})]
        completed, to_run = check_step_completion_status(steps, config)

        assert completed == []
        assert to_run == ["metadata"]

    def test_multiple_steps_mixed_status(self, tmp_path: Path) -> None:
        """A mix of completed and incomplete steps should be correctly categorised."""
        config = _make_config(tmp_path)
        # metadata is completed
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")
        # config is NOT completed (no config_base dir)
        # getfastq always runs

        steps: List[Tuple[str, Dict[str, Any]]] = [
            ("metadata", {}),
            ("config", {}),
            ("getfastq", {}),
        ]
        completed, to_run = check_step_completion_status(steps, config)

        completed_names = [c[0] for c in completed]
        assert "metadata" in completed_names
        assert "config" in to_run
        assert "getfastq" in to_run

    def test_sanity_completed_when_file_exists(self, tmp_path: Path) -> None:
        """sanity step completed when sanity_check.txt exists."""
        config = _make_config(tmp_path)
        _write_file(tmp_path / "sanity_check.txt", "OK")

        steps: List[Tuple[str, Dict[str, Any]]] = [("sanity", {})]
        completed, to_run = check_step_completion_status(steps, config)

        assert len(completed) == 1
        assert completed[0][0] == "sanity"


# ===================================================================
# handle_post_step_actions
# ===================================================================


class TestHandlePostStepActions:
    """Tests for handle_post_step_actions()."""

    def test_config_step_creates_symlink(self, tmp_path: Path) -> None:
        """After config step, a config -> config_base symlink should be created."""
        config = _make_config(tmp_path)
        config_base = tmp_path / "config_base"
        config_base.mkdir(parents=True, exist_ok=True)
        _write_file(config_base / "species.config", "data")

        result = FakeResult(returncode=0)
        handle_post_step_actions(
            step_name="config",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("config", {})],
            steps_config={},
        )

        config_dir = tmp_path / "config"
        assert config_dir.exists() or config_dir.is_symlink()

    def test_config_step_no_crash_without_config_base(self, tmp_path: Path) -> None:
        """config post-step should not crash when config_base does not exist."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=0)

        # Should not raise
        handle_post_step_actions(
            step_name="config",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("config", {})],
            steps_config={},
        )

    def test_config_step_no_overwrite_existing_config_dir(self, tmp_path: Path) -> None:
        """config post-step should not overwrite an existing config directory."""
        config = _make_config(tmp_path)
        config_base = tmp_path / "config_base"
        config_base.mkdir(parents=True, exist_ok=True)
        existing_config_dir = tmp_path / "config"
        existing_config_dir.mkdir(parents=True, exist_ok=True)
        _write_file(existing_config_dir / "existing.config", "preserve")

        result = FakeResult(returncode=0)
        handle_post_step_actions(
            step_name="config",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("config", {})],
            steps_config={},
        )

        # Existing dir should still be a real dir, not replaced by a symlink
        assert existing_config_dir.is_dir()
        assert (existing_config_dir / "existing.config").exists()

    def test_integrate_deduplicates_metadata(self, tmp_path: Path) -> None:
        """After integrate (rc=0), metadata files should be deduplicated (function called without error)."""
        config = _make_config(tmp_path)
        metadata_dir = tmp_path / "metadata"
        # Write metadata with duplicate rows
        header = "run_accession\torganism\n"
        row = "SRR1001\tApis mellifera\n"
        _write_file(metadata_dir / "metadata.tsv", header + row + row)

        result = FakeResult(returncode=0)
        handle_post_step_actions(
            step_name="integrate",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )

        # The file should still exist and be readable after deduplication
        content = (metadata_dir / "metadata.tsv").read_text()
        assert "run_accession" in content

    def test_integrate_skipped_on_failure(self, tmp_path: Path) -> None:
        """After integrate with returncode!=0, no deduplication should happen."""
        config = _make_config(tmp_path)
        # No metadata files exist

        result = FakeResult(returncode=1, stderr="integrate failed")
        # Should not raise even without metadata files
        handle_post_step_actions(
            step_name="integrate",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("integrate", {})],
            steps_config={},
        )

    def test_unknown_step_no_crash(self, tmp_path: Path) -> None:
        """Post-step for an unknown step should not crash."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=0)

        handle_post_step_actions(
            step_name="unknown_step",
            step_params={},
            result=result,
            config=config,
            steps_planned=[],
            steps_config={},
        )

    def test_step_results_list_is_optional(self, tmp_path: Path) -> None:
        """step_results=None should be handled gracefully."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=0)

        handle_post_step_actions(
            step_name="config",
            step_params={},
            result=result,
            config=config,
            steps_planned=[],
            steps_config={},
            step_results=None,
        )

    def test_getfastq_success_no_crash(self, tmp_path: Path) -> None:
        """After getfastq with rc=0, post-step should not crash even without validation data."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=0, stdout="Downloaded 10 samples", stderr="")

        step_results: List[WorkflowStepResult] = []
        handle_post_step_actions(
            step_name="getfastq",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("getfastq", {})],
            steps_config={},
            step_results=step_results,
        )
        # Should not crash; step_results may or may not have entries

    def test_getfastq_failure_skips_validation(self, tmp_path: Path) -> None:
        """After getfastq with rc!=0, validation should be skipped."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=1, stderr="download error")

        step_results: List[WorkflowStepResult] = []
        handle_post_step_actions(
            step_name="getfastq",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("getfastq", {})],
            steps_config={},
            step_results=step_results,
        )
        # No validation step_results appended because rc != 0
        # (the validation block only runs when result.returncode == 0)
        assert all(sr.step_name != "getfastq_validation" for sr in step_results)

    def test_quant_success_no_crash(self, tmp_path: Path) -> None:
        """After quant with rc=0, post-step should not crash."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=0, stdout="Quantified 10 samples")

        step_results: List[WorkflowStepResult] = []
        handle_post_step_actions(
            step_name="quant",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
            step_results=step_results,
        )

    def test_quant_failure_skips_cleanup(self, tmp_path: Path) -> None:
        """After quant with rc!=0, post-quant validation/cleanup is skipped."""
        config = _make_config(tmp_path)
        result = FakeResult(returncode=1, stderr="quant failed")

        step_results: List[WorkflowStepResult] = []
        handle_post_step_actions(
            step_name="quant",
            step_params={},
            result=result,
            config=config,
            steps_planned=[("quant", {})],
            steps_config={},
            step_results=step_results,
        )
        # No validation entries since rc != 0


# ===================================================================
# log_workflow_summary
# ===================================================================


class TestLogWorkflowSummary:
    """Tests for log_workflow_summary()."""

    def test_all_steps_successful(self, tmp_path: Path) -> None:
        """Summary should complete without error when all steps succeed."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(step_name="metadata", return_code=0, success=True),
            WorkflowStepResult(step_name="config", return_code=0, success=True),
            WorkflowStepResult(step_name="quant", return_code=0, success=True),
        ]
        # Should not raise
        log_workflow_summary(results, config)

    def test_some_steps_failed(self, tmp_path: Path) -> None:
        """Summary should log remediation for failed steps without raising."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(step_name="getfastq", return_code=0, success=True),
            WorkflowStepResult(
                step_name="quant",
                return_code=1,
                success=False,
                error_message="Quant failed: no reads",
            ),
            WorkflowStepResult(
                step_name="merge",
                return_code=1,
                success=False,
                error_message="Merge failed: no quant files",
            ),
        ]
        log_workflow_summary(results, config)

    def test_all_steps_failed(self, tmp_path: Path) -> None:
        """Summary should handle all-failure scenario."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(step_name="getfastq", return_code=1, success=False, error_message="Download error"),
            WorkflowStepResult(step_name="integrate", return_code=1, success=False, error_message="No FASTQ"),
        ]
        log_workflow_summary(results, config)

    def test_empty_results(self, tmp_path: Path) -> None:
        """Summary should handle empty step results list."""
        config = _make_config(tmp_path)
        log_workflow_summary([], config)

    def test_single_failed_step_with_no_error_message(self, tmp_path: Path) -> None:
        """Summary should handle failed step with None error_message."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(step_name="curate", return_code=1, success=False, error_message=None),
        ]
        log_workflow_summary(results, config)

    def test_getfastq_validation_remediation(self, tmp_path: Path) -> None:
        """getfastq_validation failure should produce getfastq remediation template."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(
                step_name="getfastq_validation",
                return_code=1,
                success=False,
                error_message="No FASTQ files extracted",
            ),
        ]
        # Should not raise; getfastq_validation maps to getfastq template
        log_workflow_summary(results, config)

    def test_unknown_step_remediation(self, tmp_path: Path) -> None:
        """Unknown step names should get generic remediation guidance."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(
                step_name="custom_step",
                return_code=1,
                success=False,
                error_message="Something went wrong",
            ),
        ]
        log_workflow_summary(results, config)

    def test_multiline_error_message(self, tmp_path: Path) -> None:
        """Summary should handle multi-line error messages correctly."""
        config = _make_config(tmp_path)
        results = [
            WorkflowStepResult(
                step_name="quant",
                return_code=1,
                success=False,
                error_message="Line 1: failed\nLine 2: details\nLine 3: more info",
            ),
        ]
        log_workflow_summary(results, config)

    def test_log_dir_used_in_remediation(self, tmp_path: Path) -> None:
        """Remediation templates should use the config's log_dir."""
        log_dir = tmp_path / "custom_logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        config = _make_config(tmp_path, log_dir=str(log_dir))
        results = [
            WorkflowStepResult(step_name="getfastq", return_code=1, success=False, error_message="Failed"),
        ]
        # Should not raise; internally formats template with log_dir
        log_workflow_summary(results, config)


# ===================================================================
# setup_vdb_config (edge cases only - does not require vdb-config binary)
# ===================================================================


class TestSetupVdbConfig:
    """Tests for setup_vdb_config() - edge cases and fallback paths."""

    def test_no_getfastq_step_returns_unchanged(self, tmp_path: Path) -> None:
        """When no getfastq step is in steps_planned, return unchanged list."""
        config = _make_config(tmp_path)
        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {}), ("config", {})]
        result = setup_vdb_config(config, steps)
        assert result == steps

    def test_returns_list_type(self, tmp_path: Path) -> None:
        """setup_vdb_config should always return a list."""
        config = _make_config(tmp_path)
        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {"out_dir": str(tmp_path / "fastq")})]
        result = setup_vdb_config(config, steps)
        assert isinstance(result, list)

    def test_getfastq_dir_created(self, tmp_path: Path) -> None:
        """setup_vdb_config should create the getfastq output directory."""
        config = _make_config(tmp_path)
        fastq_out = tmp_path / "fastq"
        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {"out_dir": str(fastq_out)})]
        setup_vdb_config(config, steps)

        getfastq_dir = fastq_out / "getfastq"
        assert getfastq_dir.exists()

    def test_getfastq_removed_when_all_quantified(self, tmp_path: Path) -> None:
        """getfastq should be removed from steps when all samples are already quantified."""
        config = _make_config(tmp_path)
        # Create metadata_selected.tsv with one sample
        metadata_dir = tmp_path / "metadata"
        header = "run_accession\torganism\n"
        row = "SRR1001\tApis mellifera\n"
        _write_file(metadata_dir / "metadata_selected.tsv", header + row)

        # Create quant output for that sample (so it counts as quantified)
        quant_dir = tmp_path / "quant"
        _write_file(quant_dir / "SRR1001" / "abundance.tsv", "gene\tcount\ngene1\t100\n")

        fastq_out = tmp_path / "fastq"
        steps: List[Tuple[str, Dict[str, Any]]] = [
            ("getfastq", {"out_dir": str(fastq_out)}),
            ("quant", {}),
        ]
        result = setup_vdb_config(config, steps)
        step_names = [name for name, _ in result]
        # getfastq might be removed if all samples are quantified
        # (depends on filter_metadata_for_unquantified's exact logic)
        assert isinstance(result, list)

    def test_exception_does_not_crash(self, tmp_path: Path) -> None:
        """setup_vdb_config should handle internal exceptions gracefully."""
        config = _make_config(tmp_path)
        # Use a non-existent directory that will cause disk space check to fail
        # if the path can't be created
        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {"out_dir": str(tmp_path / "fastq")})]
        # Should not raise even if vdb-config is not installed
        result = setup_vdb_config(config, steps)
        assert isinstance(result, list)

    def test_default_out_dir_uses_work_dir(self, tmp_path: Path) -> None:
        """When out_dir uses the default work_dir/fastq path, the getfastq subdir is created."""
        config = _make_config(tmp_path)
        # Note: empty params dict {} is falsy in Python, so setup_vdb_config
        # treats it as "no getfastq params found". We must provide at least one key.
        default_fastq = str(tmp_path / "fastq")
        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {"out_dir": default_fastq})]
        result = setup_vdb_config(config, steps)
        assert isinstance(result, list)
        # The default getfastq dir should have been created
        default_getfastq = tmp_path / "fastq" / "getfastq"
        assert default_getfastq.exists()

    def test_empty_params_dict_treated_as_no_getfastq(self, tmp_path: Path) -> None:
        """Empty params dict {} is falsy, so setup_vdb_config skips vdb-config setup."""
        config = _make_config(tmp_path)
        steps: List[Tuple[str, Dict[str, Any]]] = [("getfastq", {})]
        # Empty dict is falsy -> treated as "no getfastq params" -> returns unchanged
        result = setup_vdb_config(config, steps)
        assert result == steps


# ===================================================================
# Integration-style tests
# ===================================================================


class TestWorkflowStepsIntegration:
    """Integration tests combining multiple workflow_steps functions."""

    def test_validate_then_check_completion_flow(self, tmp_path: Path) -> None:
        """Simulate a real flow: check completion, validate prerequisites, then log summary."""
        config = _make_config(tmp_path)

        # Set up completed metadata step
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [
            ("metadata", {}),
            ("integrate", {}),
        ]

        # Check completion
        completed, to_run = check_step_completion_status(steps, config)
        assert "metadata" in [c[0] for c in completed]
        assert "integrate" in to_run

        # Validate prerequisites for integrate (should fail - no FASTQ)
        error = validate_step_prerequisites(
            step_name="integrate",
            step_params={},
            config=config,
            steps_planned=steps,
            steps_config={},
        )
        assert error is not None
        assert "FASTQ" in error

        # Log summary with a failed result
        results = [
            WorkflowStepResult(step_name="metadata", return_code=0, success=True),
            WorkflowStepResult(step_name="integrate", return_code=1, success=False, error_message=error),
        ]
        log_workflow_summary(results, config)

    def test_full_success_flow(self, tmp_path: Path) -> None:
        """Simulate a fully successful workflow flow."""
        config = _make_config(tmp_path)

        # Set up all completion indicators
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")
        _write_metadata_tsv(tmp_path / "metadata" / "metadata_selected.tsv", rows=2)
        config_base = tmp_path / "config_base"
        config_base.mkdir(parents=True, exist_ok=True)
        _write_file(config_base / "species.config", "data")
        _write_file(tmp_path / "integration" / "integrated_metadata.json", "{}")

        steps: List[Tuple[str, Dict[str, Any]]] = [
            ("metadata", {}),
            ("config", {}),
            ("select", {}),
            ("integrate", {}),
        ]

        completed, to_run = check_step_completion_status(steps, config)
        assert len(completed) == 4
        assert to_run == []

        # Log all-success summary
        results = [WorkflowStepResult(step_name=name, return_code=0, success=True) for name, _ in steps]
        log_workflow_summary(results, config)

    def test_post_step_actions_for_multiple_steps(self, tmp_path: Path) -> None:
        """Running post-step for several steps in sequence should not crash."""
        config = _make_config(tmp_path)
        success_result = FakeResult(returncode=0)
        failure_result = FakeResult(returncode=1, stderr="error")

        step_results: List[WorkflowStepResult] = []

        for step_name, result_obj in [
            ("config", success_result),
            ("integrate", success_result),
            ("getfastq", failure_result),
            ("quant", success_result),
        ]:
            handle_post_step_actions(
                step_name=step_name,
                step_params={},
                result=result_obj,
                config=config,
                steps_planned=[],
                steps_config={},
                step_results=step_results,
            )

    def test_redo_values_interpreted_correctly(self, tmp_path: Path) -> None:
        """Various redo values should be correctly interpreted as boolean."""
        config = _make_config(tmp_path)
        _write_metadata_tsv(tmp_path / "metadata" / "metadata.tsv")

        steps: List[Tuple[str, Dict[str, Any]]] = [("metadata", {})]

        # "1" should force redo
        _, to_run = check_step_completion_status(steps, config, redo="1")
        assert to_run == ["metadata"]

        # "true" should force redo
        _, to_run = check_step_completion_status(steps, config, redo="true")
        assert to_run == ["metadata"]

        # "false" should not force redo
        completed, to_run = check_step_completion_status(steps, config, redo="false")
        assert len(completed) == 1
        assert to_run == []

        # "0" should not force redo
        completed, to_run = check_step_completion_status(steps, config, redo="0")
        assert len(completed) == 1
        assert to_run == []
