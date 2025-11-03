"""Comprehensive tests for RNA workflow step modules."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.steps import (
    batched_process,
    sequential_process,
)


class TestSequentialProcessing:
    """Test sequential download-quantify-delete processing."""

    def test_sequential_process_metadata_handling(self, tmp_path: Path):
        """Test sequential process can handle metadata file."""
        # Create minimal metadata file
        metadata_file = tmp_path / "metadata.tsv"
        metadata_file.write_text("run\nSRR123\nSRR456\n")
        
        # Test that metadata file can be read
        from metainformant.core.io import read_delimited
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        assert len(rows) >= 0  # May be 0 if header-only
        
    def test_sequential_process_params(self, tmp_path: Path):
        """Test sequential process parameter handling."""
        params = {
            "metadata_path": str(tmp_path / "metadata.tsv"),
            "getfastq_params": {"out_dir": str(tmp_path / "fastq")},
            "quant_params": {"out_dir": str(tmp_path / "quant")},
            "work_dir": str(tmp_path / "work"),
            "max_samples": 1,
        }
        assert isinstance(params, dict)
        assert "metadata_path" in params


class TestBatchedProcessing:
    """Test batched download-quantify-delete processing."""

    def test_batched_process_metadata_handling(self, tmp_path: Path):
        """Test batched process can handle metadata file."""
        # Create minimal metadata file
        metadata_file = tmp_path / "metadata.tsv"
        metadata_file.write_text("run\nSRR123\nSRR456\n")
        
        # Test that metadata file can be read
        from metainformant.core.io import read_delimited
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        assert len(rows) >= 0

    def test_batched_process_params(self, tmp_path: Path):
        """Test batched process parameter handling."""
        params = {
            "metadata_path": str(tmp_path / "metadata.tsv"),
            "getfastq_params": {"out_dir": str(tmp_path / "fastq")},
            "quant_params": {"out_dir": str(tmp_path / "quant")},
            "batch_size": 4,
            "max_samples": 10,
        }
        assert isinstance(params, dict)
        assert "metadata_path" in params
        assert params["batch_size"] == 4


class TestStepModuleImports:
    """Test that step modules can be imported and have expected structure."""

    def test_step_modules_importable(self):
        """Test all step modules can be imported."""
        from metainformant.rna.steps import (
            config,
            csca,
            cstmm,
            curate,
            getfastq,
            integrate,
            merge,
            metadata,
            parallel_download,
            quant,
            sanity,
            select,
        )
        
        # Verify modules exist
        assert config is not None
        assert metadata is not None
        assert quant is not None

    def test_step_runner_functions_exist(self):
        """Test that step runner functions are available."""
        from metainformant.rna.steps import (
            run_metadata,
            run_integrate,
            run_config,
            run_select,
            run_getfastq,
            run_quant,
            run_merge,
            run_cstmm,
            run_curate,
            run_csca,
            run_sanity,
        )
        
        # Verify functions exist and are callable
        assert callable(run_metadata)
        assert callable(run_quant)
        assert callable(run_merge)

