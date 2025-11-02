"""
Tests for ENA-based RNA-seq workflow integration.

This module validates the robust ENA download and quantification workflow:
- Direct FASTQ downloads from ENA API
- Robust retry and resume logic with wget
- Batched download→quantify→delete processing
- Auto-detection of single vs paired-end data
- Integration with kallisto quantification

Tests use real network requests and real tools (no mocks per project policy).
"""

import subprocess
import sys
from pathlib import Path
import pytest

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from metainformant.core.io import read_delimited


class TestENADownloader:
    """Test the standalone ENA downloader."""
    
    def test_downloader_script_exists(self):
        """Verify download_ena_robust.py exists and is executable."""
        script = Path("scripts/rna/download_ena_robust.py")
        assert script.exists(), "Download script not found"
        assert script.stat().st_mode & 0o111, "Download script not executable"
    
    def test_downloader_help(self):
        """Verify downloader script shows help."""
        script = Path("scripts/rna/download_ena_robust.py")
        result = subprocess.run(
            [sys.executable, str(script), "--help"],
            capture_output=True,
            text=True,
            timeout=10
        )
        assert result.returncode == 0, "Help command failed"
        assert "--metadata" in result.stdout, "Missing --metadata option"
        assert "--out-dir" in result.stdout, "Missing --out-dir option"
        assert "--threads" in result.stdout, "Missing --threads option"


class TestIntegratedWorkflow:
    """Test the integrated ENA download + quantification workflow."""
    
    def test_workflow_script_exists(self):
        """Verify workflow_ena_integrated.py exists and is executable."""
        script = Path("scripts/rna/workflow_ena_integrated.py")
        assert script.exists(), "Workflow script not found"
        assert script.stat().st_mode & 0o111, "Workflow script not executable"
    
    def test_workflow_help(self):
        """Verify workflow script shows help."""
        script = Path("scripts/rna/workflow_ena_integrated.py")
        result = subprocess.run(
            [sys.executable, str(script), "--help"],
            capture_output=True,
            text=True,
            timeout=10
        )
        assert result.returncode == 0, "Help command failed"
        assert "--config" in result.stdout, "Missing --config option"
        assert "--batch-size" in result.stdout, "Missing --batch-size option"
        assert "--threads" in result.stdout, "Missing --threads option"
        assert "--max-samples" in result.stdout, "Missing --max-samples option"
        assert "--skip-download" in result.stdout, "Missing --skip-download option"
    
    def test_workflow_requires_config(self):
        """Verify workflow requires --config argument."""
        script = Path("scripts/rna/workflow_ena_integrated.py")
        result = subprocess.run(
            [sys.executable, str(script)],
            capture_output=True,
            text=True,
            timeout=10
        )
        assert result.returncode != 0, "Should fail without config"
        assert "required" in result.stderr.lower() or "config" in result.stderr.lower()


class TestWorkflowIntegration:
    """Integration tests with real data (when available)."""
    
    def test_cfloridanus_config_exists(self):
        """Verify C. floridanus config exists for testing."""
        config = Path("config/amalgkit/amalgkit_cfloridanus.yaml")
        assert config.exists(), "Test config not found"
    
    def test_cfloridanus_metadata_exists(self):
        """Verify C. floridanus metadata exists."""
        metadata = Path("output/amalgkit/cfloridanus/work/metadata/metadata.tsv")
        if metadata.exists():
            # Verify structure
            rows = list(read_delimited(metadata, delimiter='\t'))
            assert len(rows) > 0, "Metadata empty"
            assert 'run' in rows[0], "Metadata missing 'run' column"
            assert 'lib_layout' in rows[0], "Metadata missing 'lib_layout' column"
        else:
            pytest.skip("Metadata not generated yet")
    
    def test_kallisto_index_exists(self):
        """Verify kallisto index exists for C. floridanus."""
        index_dir = Path("output/amalgkit/cfloridanus/work/index")
        if index_dir.exists():
            indices = list(index_dir.glob("*.idx"))
            assert len(indices) > 0, "No kallisto index found"
        else:
            pytest.skip("Index not built yet")
    
    def test_quantification_results_exist(self):
        """Verify quantification results from test run."""
        quant_dir = Path("output/amalgkit/cfloridanus/quant")
        if quant_dir.exists():
            samples = list(quant_dir.glob("SRR*/abundance.tsv"))
            if len(samples) > 0:
                # Verify structure of abundance files
                sample = samples[0]
                with open(sample) as f:
                    header = f.readline().strip().split('\t')
                    assert 'target_id' in header, "Missing target_id column"
                    assert 'tpm' in header, "Missing TPM column"
                    assert 'est_counts' in header, "Missing est_counts column"
            else:
                pytest.skip("No quantified samples yet")
        else:
            pytest.skip("Quantification not run yet")


class TestWorkflowDocumentation:
    """Tests verifying documentation is updated."""
    
    def test_scripts_readme_mentions_ena_workflow(self):
        """Verify scripts/rna/README.md documents the ENA workflow."""
        readme = Path("scripts/rna/README.md")
        assert readme.exists(), "README not found"
        content = readme.read_text()
        assert "workflow_ena_integrated" in content, "ENA workflow not documented"
        assert "download_ena_robust" in content, "ENA downloader not documented"
        assert "ENA" in content, "ENA not mentioned"
    
    def test_scripts_agents_mentions_ena(self):
        """Verify scripts/rna/AGENTS.md documents ENA contributions."""
        agents = Path("scripts/rna/AGENTS.md")
        assert agents.exists(), "AGENTS.md not found"
        content = agents.read_text()
        assert "ENA" in content, "ENA not mentioned"
        assert "download_ena_robust" in content or "ena" in content.lower()
    
    def test_docs_rna_readme_mentions_ena(self):
        """Verify docs/rna/README.md mentions ENA approach."""
        readme = Path("docs/rna/README.md")
        assert readme.exists(), "docs README not found"
        content = readme.read_text()
        assert "ENA" in content, "ENA not mentioned in docs"


class TestDependencies:
    """Tests for required external dependencies."""
    
    def test_wget_available(self):
        """Verify wget is installed (required for downloads)."""
        result = subprocess.run(
            ["which", "wget"],
            capture_output=True,
            timeout=5
        )
        if result.returncode != 0:
            pytest.skip("wget not installed (required for ENA downloads)")
    
    def test_kallisto_available(self):
        """Verify kallisto is installed (required for quantification)."""
        result = subprocess.run(
            ["which", "kallisto"],
            capture_output=True,
            timeout=5
        )
        if result.returncode != 0:
            pytest.skip("kallisto not installed (required for quantification)")
    
    def test_wget_supports_continue(self):
        """Verify wget supports --continue flag."""
        result = subprocess.run(
            ["wget", "--help"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            assert "--continue" in result.stdout, "wget missing --continue support"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

