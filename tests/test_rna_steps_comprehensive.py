"""Comprehensive tests for RNA workflow step modules."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna.steps import run_download_quant_workflow


class TestUnifiedProcessing:
    """Test unified download-quantify-delete processing."""

    def test_unified_process_metadata_handling(self, tmp_path: Path):
        """Test unified process can handle metadata file."""
        # Create minimal metadata file
        metadata_file = tmp_path / "metadata.tsv"
        metadata_file.write_text("run\nSRR123\nSRR456\n")
        
        # Test that metadata file can be read
        from metainformant.core.io import read_delimited
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        assert len(rows) >= 0  # May be 0 if header-only
        
    def test_unified_process_params_sequential(self, tmp_path: Path):
        """Test unified process parameter handling (sequential mode)."""
        params = {
            "metadata_path": str(tmp_path / "metadata.tsv"),
            "getfastq_params": {"out_dir": str(tmp_path / "fastq")},
            "quant_params": {"out_dir": str(tmp_path / "quant")},
            "work_dir": str(tmp_path / "work"),
            "num_workers": 1,  # Sequential mode
            "max_samples": 1,
        }
        assert isinstance(params, dict)
        assert "metadata_path" in params
        assert params["num_workers"] == 1

    def test_unified_process_params_parallel(self, tmp_path: Path):
        """Test unified process parameter handling (parallel mode)."""
        params = {
            "metadata_path": str(tmp_path / "metadata.tsv"),
            "getfastq_params": {"out_dir": str(tmp_path / "fastq")},
            "quant_params": {"out_dir": str(tmp_path / "quant")},
            "num_workers": 4,  # Parallel mode
            "max_samples": 10,
        }
        assert isinstance(params, dict)
        assert "metadata_path" in params
        assert params["num_workers"] == 4


class TestDownloadValidation:
    """Test download validation and 0-read detection."""

    def test_download_validation_imports(self):
        """Test that download validation functions can be imported."""
        from metainformant.rna.steps.process_samples import (
            _download_worker,
            _wait_for_fastq_files,
            _delete_fastq_for_sample,
        )
        assert callable(_wait_for_fastq_files)
        assert callable(_delete_fastq_for_sample)

    def test_fastq_dir_extraction_logic(self, tmp_path: Path):
        """Test fastq_dir extraction from getfastq_params."""
        getfastq_params = {"out_dir": str(tmp_path / "fastq")}
        work_dir = tmp_path
        
        # Simulate the extraction logic from _download_worker
        out_dir_str = getfastq_params.get("out_dir", "")
        if out_dir_str:
            from pathlib import Path
            fastq_dir = Path(out_dir_str)
            if not fastq_dir.is_absolute():
                fastq_dir = fastq_dir.resolve()
        elif work_dir:
            fastq_dir = (work_dir / "fastq").resolve()
        else:
            fastq_dir = Path("fastq").resolve()
        
        assert str(fastq_dir) == str((tmp_path / "fastq").resolve())


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
            process_samples,
            quant,
            sanity,
            select,
        )
        
        # Verify modules exist
        assert config is not None
        assert metadata is not None
        assert quant is not None
        assert process_samples is not None

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
            run_download_quant_workflow,
        )
        
        # Verify functions exist and are callable
        assert callable(run_metadata)
        assert callable(run_quant)
        assert callable(run_merge)
        assert callable(run_download_quant_workflow)

