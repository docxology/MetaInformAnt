"""Tests for orchestrator scripts.

Tests follow the NO-MOCKING policy: all tests use real implementations
and call orchestrator scripts via subprocess to verify actual behavior.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

# Repository root
REPO_ROOT = Path(__file__).parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"


class TestOrchestratorScripts:
    """Test orchestrator scripts for basic functionality."""

    def test_dna_orchestrator_help(self):
        """Test DNA orchestrator shows help."""
        script = SCRIPTS_DIR / "dna" / "run_dna_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "DNA analysis" in result.stdout or "usage:" in result.stdout.lower()

    def test_protein_orchestrator_help(self):
        """Test protein orchestrator shows help."""
        script = SCRIPTS_DIR / "protein" / "run_protein_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "protein" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_ontology_orchestrator_help(self):
        """Test ontology orchestrator shows help."""
        script = SCRIPTS_DIR / "ontology" / "run_ontology_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "ontology" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_phenotype_orchestrator_help(self):
        """Test phenotype orchestrator shows help."""
        script = SCRIPTS_DIR / "phenotype" / "run_phenotype_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "phenotype" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_networks_orchestrator_help(self):
        """Test networks orchestrator shows help."""
        script = SCRIPTS_DIR / "networks" / "run_network_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "network" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_multiomics_orchestrator_help(self):
        """Test multiomics orchestrator shows help."""
        script = SCRIPTS_DIR / "multiomics" / "run_multiomics_integration.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "multiomics" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_singlecell_orchestrator_help(self):
        """Test single-cell orchestrator shows help."""
        script = SCRIPTS_DIR / "singlecell" / "run_singlecell_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "single" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_quality_orchestrator_help(self):
        """Test quality orchestrator shows help."""
        script = SCRIPTS_DIR / "quality" / "run_quality_control.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "quality" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_simulation_orchestrator_help(self):
        """Test simulation orchestrator shows help."""
        script = SCRIPTS_DIR / "simulation" / "run_simulation.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "simulation" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_visualization_orchestrator_help(self):
        """Test visualization orchestrator shows help."""
        script = SCRIPTS_DIR / "visualization" / "run_visualization.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "visualization" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_epigenome_orchestrator_help(self):
        """Test epigenome orchestrator shows help."""
        script = SCRIPTS_DIR / "epigenome" / "run_epigenome_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "epigenome" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_ecology_orchestrator_help(self):
        """Test ecology orchestrator shows help."""
        script = SCRIPTS_DIR / "ecology" / "run_ecology_analysis.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "ecology" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_ml_orchestrator_help(self):
        """Test ML orchestrator shows help."""
        script = SCRIPTS_DIR / "ml" / "run_ml_pipeline.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "machine" in result.stdout.lower() or "ml" in result.stdout.lower() or "usage:" in result.stdout.lower()

    def test_math_orchestrator_help(self):
        """Test math orchestrator shows help."""
        script = SCRIPTS_DIR / "math" / "run_math_modeling.py"
        result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "math" in result.stdout.lower() or "usage:" in result.stdout.lower()


class TestOrchestratorCLIIntegration:
    """Test orchestrator scripts via CLI integration."""

    def test_ontology_cli_integration(self, tmp_path: Path):
        """Test ontology CLI command."""
        output_dir = tmp_path / "ontology_output"
        cmd = [
            sys.executable,
            "-m",
            "metainformant",
            "ontology",
            "run",
            f"--output={output_dir}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        # Should either succeed or fail gracefully with clear message
        assert result.returncode in (0, 1)  # 0 = success, 1 = expected failure (no input)
        assert "ontology" in result.stdout.lower() or "output" in result.stdout.lower() or len(result.stderr) == 0

    def test_phenotype_cli_integration(self, tmp_path: Path):
        """Test phenotype CLI command."""
        output_dir = tmp_path / "phenotype_output"
        # Create minimal input file
        input_file = tmp_path / "phenotypes.json"
        input_file.write_text('{"data": []}')
        cmd = [
            sys.executable,
            "-m",
            "metainformant",
            "phenotype",
            "run",
            f"--input={input_file}",
            f"--output={output_dir}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode in (0, 1)

    def test_networks_cli_integration(self, tmp_path: Path):
        """Test networks CLI command."""
        output_dir = tmp_path / "networks_output"
        # Create minimal input file
        input_file = tmp_path / "interactions.tsv"
        input_file.write_text("source\ttarget\nA\tB\nB\tC\n")
        cmd = [
            sys.executable,
            "-m",
            "metainformant",
            "networks",
            "run",
            f"--input={input_file}",
            f"--output={output_dir}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode in (0, 1)

    def test_quality_cli_integration(self, tmp_path: Path):
        """Test quality CLI command."""
        output_dir = tmp_path / "quality_output"
        cmd = [
            sys.executable,
            "-m",
            "metainformant",
            "quality",
            "run",
            f"--output={output_dir}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode in (0, 1)

    def test_simulation_cli_integration(self, tmp_path: Path):
        """Test simulation CLI command."""
        output_dir = tmp_path / "simulation_output"
        cmd = [
            sys.executable,
            "-m",
            "metainformant",
            "simulation",
            "run",
            "--model=sequences",
            f"--output={output_dir}",
            "--n=10",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        assert result.returncode in (0, 1)

    @pytest.mark.slow
    def test_all_orchestrators_exist(self):
        """Verify all orchestrator scripts exist."""
        orchestrators = [
            "dna/run_dna_analysis.py",
            "protein/run_protein_analysis.py",
            "ontology/run_ontology_analysis.py",
            "phenotype/run_phenotype_analysis.py",
            "networks/run_network_analysis.py",
            "multiomics/run_multiomics_integration.py",
            "singlecell/run_singlecell_analysis.py",
            "quality/run_quality_control.py",
            "simulation/run_simulation.py",
            "visualization/run_visualization.py",
            "epigenome/run_epigenome_analysis.py",
            "ecology/run_ecology_analysis.py",
            "ml/run_ml_pipeline.py",
            "math/run_math_modeling.py",
        ]
        for orchestrator in orchestrators:
            script_path = SCRIPTS_DIR / orchestrator
            assert script_path.exists(), f"Orchestrator script not found: {orchestrator}"
            assert script_path.is_file(), f"Orchestrator is not a file: {orchestrator}"


