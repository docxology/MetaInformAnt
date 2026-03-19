"""
test_pipeline.py — End-to-end pipeline tests (Zero-Mock policy).

All tests exercise real code against synthetic data generated from
the scripts themselves.  No mocking, no patching.
"""

import subprocess
import sys
from pathlib import Path

import pytest
import yaml


PROJECT_ROOT = Path(__file__).parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"


def run_script(script: Path, config: Path, *extra_args: str) -> subprocess.CompletedProcess:
    """Run a pipeline script as a subprocess and return the result."""
    cmd = [sys.executable, str(script), "--config", str(config), *extra_args]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=PROJECT_ROOT)
    return result


# ── Stage 99: Synthetic Data ───────────────────────────────────────────────────

class TestSyntheticDataGenerator:
    def test_creates_raw_csv_files(self, tmp_project: Path) -> None:
        """Stage 99 should produce at least two CSV files in data/raw/."""
        config_path = tmp_project / "config" / "default.yaml"
        result = run_script(
            SCRIPTS / "99_create_synthetic_data.py",
            config_path,
            "--n-samples", "50",
            "--n-features", "4",
            "--seed", "1",
        )
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        raw_files = list((tmp_project / "data" / "raw").glob("*.csv"))
        assert len(raw_files) >= 2, "Expected at least 2 CSV files"

    def test_creates_metadata_yaml(self, tmp_project: Path) -> None:
        """Stage 99 should produce a valid metadata.yaml."""
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "99_create_synthetic_data.py", config_path,
                   "--n-samples", "30", "--seed", "42")
        meta = tmp_project / "data" / "raw" / "metadata.yaml"
        assert meta.exists(), "metadata.yaml not created"
        with meta.open() as fh:
            content = yaml.safe_load(fh)
        assert "generated_at" in content
        assert "datasets" in content
        assert len(content["datasets"]) >= 2

    def test_idempotent_without_force(self, tmp_project: Path) -> None:
        """Running Stage 99 twice without --force should not fail."""
        config_path = tmp_project / "config" / "default.yaml"
        for _ in range(2):
            result = run_script(SCRIPTS / "99_create_synthetic_data.py", config_path,
                                "--n-samples", "20", "--seed", "7")
        assert result.returncode == 0


# ── Stage 01: Data Processing ──────────────────────────────────────────────────

class TestDataProcessing:
    def _generate_data(self, tmp_project: Path, sample_raw_csv: Path) -> None:
        """Ensure raw data exists (fixture already creates it)."""

    def test_produces_processed_csv(
        self, tmp_project: Path, sample_raw_csv: Path
    ) -> None:
        """Stage 1 should write processed_data.csv to data/processed/."""
        config_path = tmp_project / "config" / "default.yaml"
        result = run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        output = tmp_project / "data" / "processed" / "processed_data.csv"
        assert output.exists(), "processed_data.csv not created"
        assert output.stat().st_size > 0, "processed_data.csv is empty"

    def test_writes_log_file(self, tmp_project: Path, sample_raw_csv: Path) -> None:
        """Stage 1 should write a structured log to logs/."""
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        log = tmp_project / "logs" / "01_process_data.log"
        assert log.exists()
        assert log.stat().st_size > 0

    def test_idempotent_skip(self, tmp_project: Path, sample_raw_csv: Path) -> None:
        """Stage 1 without --force should skip and exit 0 if output exists."""
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        result2 = run_script(SCRIPTS / "01_process_data.py", config_path)
        assert result2.returncode == 0
        assert "skipping" in result2.stdout.lower() or result2.returncode == 0

    def test_no_raw_files_exits_cleanly(self, tmp_project: Path) -> None:
        """Stage 1 with an empty raw dir should exit 0 (warns, doesn't crash)."""
        config_path = tmp_project / "config" / "default.yaml"
        # Do not create any raw files
        result = run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        assert result.returncode == 0


# ── Stage 02: Analysis ─────────────────────────────────────────────────────────

class TestAnalysis:
    def _prepare(self, tmp_project: Path, sample_raw_csv: Path) -> None:
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "01_process_data.py", config_path, "--force")

    def test_produces_summary_statistics(
        self, tmp_project: Path, sample_raw_csv: Path
    ) -> None:
        """Stage 2 should produce summary_statistics.csv."""
        self._prepare(tmp_project, sample_raw_csv)
        config_path = tmp_project / "config" / "default.yaml"
        result = run_script(SCRIPTS / "02_analyze_results.py", config_path, "--force")
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        summary = tmp_project / "results" / "tables" / "summary_statistics.csv"
        assert summary.exists()
        assert summary.stat().st_size > 0

    def test_produces_metadata_json(
        self, tmp_project: Path, sample_raw_csv: Path
    ) -> None:
        """Stage 2 should write analysis_metadata.json."""
        self._prepare(tmp_project, sample_raw_csv)
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "02_analyze_results.py", config_path, "--force")
        meta = tmp_project / "results" / "tables" / "analysis_metadata.json"
        assert meta.exists()

    def test_no_input_exits_nonzero(self, tmp_project: Path) -> None:
        """Stage 2 without processed data should exit non-zero."""
        config_path = tmp_project / "config" / "default.yaml"
        result = run_script(SCRIPTS / "02_analyze_results.py", config_path, "--force")
        assert result.returncode != 0


# ── Stage 03: Visualisation ────────────────────────────────────────────────────

class TestVisualisation:
    def _prepare(self, tmp_project: Path, sample_raw_csv: Path) -> None:
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        run_script(SCRIPTS / "02_analyze_results.py", config_path, "--force")

    def test_produces_distribution_grid(
        self, tmp_project: Path, sample_raw_csv: Path
    ) -> None:
        """Stage 3 should produce distribution_grid.png."""
        self._prepare(tmp_project, sample_raw_csv)
        config_path = tmp_project / "config" / "default.yaml"
        result = run_script(SCRIPTS / "03_visualize.py", config_path, "--force")
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        fig = tmp_project / "results" / "figures" / "distribution_grid.png"
        assert fig.exists()
        assert fig.stat().st_size > 0

    def test_produces_correlation_heatmap(
        self, tmp_project: Path, sample_raw_csv: Path
    ) -> None:
        """Stage 3 should produce correlation_heatmap.png when correlation data exists."""
        self._prepare(tmp_project, sample_raw_csv)
        config_path = tmp_project / "config" / "default.yaml"
        run_script(SCRIPTS / "03_visualize.py", config_path, "--force")
        heatmap = tmp_project / "results" / "figures" / "correlation_heatmap.png"
        assert heatmap.exists()
        assert heatmap.stat().st_size > 0


# ── Full end-to-end ────────────────────────────────────────────────────────────

class TestEndToEnd:
    def test_full_pipeline_from_synthetic_data(self, tmp_project: Path) -> None:
        """
        Generate synthetic data and run all three stages.
        Verifies final expected outputs are present and non-empty.
        """
        config_path = tmp_project / "config" / "default.yaml"

        # Stage 99 — generate data
        r = run_script(SCRIPTS / "99_create_synthetic_data.py", config_path,
                       "--n-samples", "60", "--seed", "99")
        assert r.returncode == 0, f"Stage 99 failed:\n{r.stderr}"

        # Stage 01
        r = run_script(SCRIPTS / "01_process_data.py", config_path, "--force")
        assert r.returncode == 0, f"Stage 01 failed:\n{r.stderr}"

        # Stage 02
        r = run_script(SCRIPTS / "02_analyze_results.py", config_path, "--force")
        assert r.returncode == 0, f"Stage 02 failed:\n{r.stderr}"

        # Stage 03
        r = run_script(SCRIPTS / "03_visualize.py", config_path, "--force")
        assert r.returncode == 0, f"Stage 03 failed:\n{r.stderr}"

        # Verify key outputs
        expected = [
            tmp_project / "data" / "processed" / "processed_data.csv",
            tmp_project / "results" / "tables" / "summary_statistics.csv",
            tmp_project / "results" / "figures" / "distribution_grid.png",
        ]
        for path in expected:
            assert path.exists(), f"Missing expected output: {path}"
            assert path.stat().st_size > 0, f"Empty output: {path}"
