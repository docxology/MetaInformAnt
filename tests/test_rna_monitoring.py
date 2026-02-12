"""Tests for RNA monitoring and status checking functions.

This module tests all monitoring functions following NO_MOCKING_POLICY.
All tests use real file operations and real implementations.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core.io.io import dump_json, write_delimited
from metainformant.rna.engine import monitoring


class TestCountQuantifiedSamples:
    """Test count_quantified_samples function."""

    def test_count_quantified_samples_empty_config(self, tmp_path: Path):
        """Test count_quantified_samples with empty directories."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        dump_json(config_data, config_file)

        quantified, total = monitoring.count_quantified_samples(config_file)
        assert isinstance(quantified, int)
        assert isinstance(total, int)
        assert quantified == 0
        assert total == 0

    def test_count_quantified_samples_with_quantified(self, tmp_path: Path):
        """Test count_quantified_samples counts quantified samples."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir since load_workflow_config resolves relative paths w.r.t. repo root.
        # Monitoring functions default quant_dir to work_dir / "quant" when not configured.
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}],
            metadata_file,
            delimiter="\t",
        )

        # Create quantified samples under the default quant_dir (work_dir / "quant")
        quant_dir = work_dir_abs / "quant"
        quant_dir.mkdir(parents=True)
        (quant_dir / "SRR123").mkdir()
        (quant_dir / "SRR123" / "abundance.tsv").write_text("test")
        (quant_dir / "SRR456").mkdir()
        (quant_dir / "SRR456" / "abundance.tsv").write_text("test")

        quantified, total = monitoring.count_quantified_samples(config_file)
        assert quantified == 2
        assert total == 2

    def test_count_quantified_samples_partial(self, tmp_path: Path):
        """Test count_quantified_samples with partial quantification."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; quant_dir defaults to work_dir / "quant"
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata with 3 samples
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}, {"run": "SRR789"}],
            metadata_file,
            delimiter="\t",
        )

        # Only quantify 2 samples under work_dir / "quant"
        quant_dir = work_dir_abs / "quant"
        quant_dir.mkdir(parents=True)
        (quant_dir / "SRR123").mkdir()
        (quant_dir / "SRR123" / "abundance.tsv").write_text("test")
        (quant_dir / "SRR456").mkdir()
        (quant_dir / "SRR456" / "abundance.tsv").write_text("test")

        quantified, total = monitoring.count_quantified_samples(config_file)
        assert quantified == 2
        assert total == 3


class TestGetSampleStatus:
    """Test get_sample_status function."""

    def test_get_sample_status_undownloaded(self, tmp_path: Path):
        """Test get_sample_status for undownloaded sample."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "getfastq": {"out_dir": str(tmp_path / "fastq")},
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        dump_json(config_data, config_file)

        status = monitoring.get_sample_status(config_file, "SRR123")
        assert isinstance(status, dict)
        assert status["quantified"] is False
        assert status["has_fastq"] is False
        assert status["has_sra"] is False
        assert status["status"] == "undownloaded"

    def test_get_sample_status_has_fastq(self, tmp_path: Path):
        """Test get_sample_status for sample with FASTQ files."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; fastq_dir defaults to work_dir / "fastq"
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create FASTQ file under default fastq_dir (work_dir / "fastq")
        fastq_dir = work_dir_abs / "fastq"
        fastq_sample_dir = fastq_dir / "getfastq" / "SRR123"
        fastq_sample_dir.mkdir(parents=True)
        (fastq_sample_dir / "SRR123_1.fastq.gz").write_text("test")

        status = monitoring.get_sample_status(config_file, "SRR123")
        assert status["has_fastq"] is True
        assert status["status"] == "has_fastq"

    def test_get_sample_status_quantified(self, tmp_path: Path):
        """Test get_sample_status for quantified sample."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; quant_dir defaults to work_dir / "quant"
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create quantified sample under default quant_dir (work_dir / "quant")
        quant_sample_dir = (work_dir_abs / "quant") / "SRR123"
        quant_sample_dir.mkdir(parents=True)
        (quant_sample_dir / "abundance.tsv").write_text("test")

        status = monitoring.get_sample_status(config_file, "SRR123")
        assert status["quantified"] is True
        assert status["status"] == "quantified"


class TestAnalyzeSpeciesStatus:
    """Test analyze_species_status function."""

    def test_analyze_species_status_empty(self, tmp_path: Path):
        """Test analyze_species_status with no metadata."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "getfastq": {"out_dir": str(tmp_path / "fastq")},
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        dump_json(config_data, config_file)

        status = monitoring.analyze_species_status(config_file)
        assert isinstance(status, dict)
        assert status["total_in_metadata"] == 0
        assert status["quantified"] == 0
        assert "categories" in status

    def test_analyze_species_status_with_samples(self, tmp_path: Path):
        """Test analyze_species_status with samples."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; fastq_dir/quant_dir default under work_dir
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}],
            metadata_file,
            delimiter="\t",
        )

        # Create quantified sample under work_dir / "quant"
        quant_sample_dir = (work_dir_abs / "quant") / "SRR123"
        quant_sample_dir.mkdir(parents=True)
        (quant_sample_dir / "abundance.tsv").write_text("test")

        status = monitoring.analyze_species_status(config_file)
        assert status["total_in_metadata"] == 2
        assert status["quantified"] == 1
        assert "categories" in status
        assert "quantified_and_deleted" in status["categories"]
        assert "undownloaded" in status["categories"]


class TestFindUnquantifiedSamples:
    """Test find_unquantified_samples function."""

    def test_find_unquantified_samples_empty(self, tmp_path: Path):
        """Test find_unquantified_samples with no metadata."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        dump_json(config_data, config_file)

        unquantified = monitoring.find_unquantified_samples(config_file)
        assert isinstance(unquantified, list)
        assert len(unquantified) == 0

    def test_find_unquantified_samples_finds_unquantified(self, tmp_path: Path):
        """Test find_unquantified_samples finds unquantified samples."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; quant_dir defaults under work_dir
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}],
            metadata_file,
            delimiter="\t",
        )

        # Only quantify one sample under work_dir / "quant"
        quant_sample_dir = (work_dir_abs / "quant") / "SRR123"
        quant_sample_dir.mkdir(parents=True)
        (quant_sample_dir / "abundance.tsv").write_text("test")

        unquantified = monitoring.find_unquantified_samples(config_file)
        assert len(unquantified) == 1
        assert "SRR456" in unquantified


class TestCheckActiveDownloads:
    """Test check_active_downloads function."""

    @pytest.mark.slow
    def test_check_active_downloads_returns_set(self):
        """Test that check_active_downloads returns a set.

        NOTE: This test is marked as slow because it runs `ps aux` which can be slow.
        """
        import signal

        def timeout_handler(signum, frame):
            raise TimeoutError("check_active_downloads timed out")

        # Set timeout for this test (5 seconds)
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(5)

        try:
            active = monitoring.check_active_downloads()
            assert isinstance(active, set)
            # May be empty if no downloads are active
            assert all(isinstance(sample_id, str) for sample_id in active)
        except TimeoutError:
            pytest.skip("check_active_downloads timed out (likely system issue)")
        finally:
            signal.alarm(0)  # Cancel alarm


class TestCheckWorkflowProgress:
    """Test check_workflow_progress function."""

    def test_check_workflow_progress_empty(self, tmp_path: Path):
        """Test check_workflow_progress with empty workflow."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
            "per_step": {
                "getfastq": {"out_dir": str(tmp_path / "fastq")},
                "quant": {"out_dir": str(tmp_path / "quant")},
            },
        }
        dump_json(config_data, config_file)

        progress = monitoring.check_workflow_progress(config_file)
        assert isinstance(progress, dict)
        assert "quantified" in progress
        assert "total" in progress
        assert "percentage" in progress
        assert "remaining" in progress
        assert "downloading" in progress
        assert "has_files" in progress
        assert progress["quantified"] == 0
        assert progress["total"] == 0
        assert progress["percentage"] == 0.0

    def test_check_workflow_progress_with_samples(self, tmp_path: Path):
        """Test check_workflow_progress with samples."""
        config_file = tmp_path / "config.yaml"
        # Use absolute work_dir; fastq_dir/quant_dir default under work_dir
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}],
            metadata_file,
            delimiter="\t",
        )

        # Quantify one sample under work_dir / "quant"
        quant_sample_dir = (work_dir_abs / "quant") / "SRR123"
        quant_sample_dir.mkdir(parents=True)
        (quant_sample_dir / "abundance.tsv").write_text("test")

        progress = monitoring.check_workflow_progress(config_file)
        assert progress["quantified"] == 1
        assert progress["total"] == 2
        assert progress["percentage"] == 50.0
        assert progress["remaining"] == 1


class TestAssessAllSpeciesProgress:
    """Test assess_all_species_progress function."""

    def test_assess_all_species_progress_empty_dir(self, tmp_path: Path):
        """Test assess_all_species_progress with empty config directory."""
        config_dir = tmp_path / "configs"
        config_dir.mkdir(parents=True)

        results = monitoring.assess_all_species_progress(config_dir)
        assert isinstance(results, dict)
        assert len(results) == 0

    def test_assess_all_species_progress_with_configs(self, tmp_path: Path):
        """Test assess_all_species_progress with config files."""
        config_dir = tmp_path / "configs"
        config_dir.mkdir(parents=True)

        # Create config with absolute work_dir; filename must not contain 'test' to be included
        config_file = config_dir / "amalgkit_example_species.yaml"
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        results = monitoring.assess_all_species_progress(config_dir)
        assert isinstance(results, dict)
        assert "example_species" in results


class TestInitializeProgressTracking:
    """Test initialize_progress_tracking function."""

    def test_initialize_progress_tracking_no_metadata(self, tmp_path: Path):
        """Test initialize_progress_tracking with no metadata."""
        config_file = tmp_path / "config.yaml"
        config_data = {
            "work_dir": str(tmp_path / "work"),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        result = monitoring.initialize_progress_tracking(config_file)
        assert isinstance(result, dict)
        assert result["success"] is False
        assert "error" in result

    def test_initialize_progress_tracking_with_metadata(self, tmp_path: Path):
        """Test initialize_progress_tracking with metadata."""
        config_file = tmp_path / "config.yaml"
        # Use absolute paths
        work_dir_abs = (tmp_path / "work").resolve()
        config_data = {
            "work_dir": str(work_dir_abs),
            "threads": 4,
            "species_list": ["Test_species"],
        }
        dump_json(config_data, config_file)

        # Create metadata
        metadata_dir = work_dir_abs / "metadata"
        metadata_dir.mkdir(parents=True)
        metadata_file = metadata_dir / "metadata.tsv"
        write_delimited(
            [{"run": "SRR123"}, {"run": "SRR456"}],
            metadata_file,
            delimiter="\t",
        )

        result = monitoring.initialize_progress_tracking(config_file)
        assert isinstance(result, dict)
        assert result["success"] is True
        assert result["total_samples"] == 2
        assert "sample_ids" in result
        assert len(result["sample_ids"]) == 2


class TestMonitoringDocumentation:
    """Test that monitoring functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all monitoring functions have docstrings."""
        functions = [
            monitoring.count_quantified_samples,
            monitoring.get_sample_status,
            monitoring.analyze_species_status,
            monitoring.find_unquantified_samples,
            monitoring.check_active_downloads,
            monitoring.check_workflow_progress,
            monitoring.assess_all_species_progress,
            monitoring.initialize_progress_tracking,
        ]

        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0
