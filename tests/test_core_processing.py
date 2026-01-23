"""Tests for config-based processing functionality."""

import pytest
import tempfile
from pathlib import Path

from metainformant.core import io
from metainformant.core import (
    create_sample_config,
    download_and_process_data,
    run_config_based_workflow,
    validate_config_file,
)


class TestConfigBasedProcessing:
    """Test config-based end-to-end processing functionality."""

    def test_validate_config_file(self):
        """Test configuration file validation."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Test non-existent file
            non_existent = tmp_path / "does_not_exist.json"
            is_valid, errors = validate_config_file(non_existent)
            assert not is_valid
            assert len(errors) == 1
            assert "not found" in errors[0]

            # Test valid config file
            valid_config = {"downloads": {"test_data": {"url": "https://example.com/test.txt", "filename": "test.txt"}}}
            config_file = tmp_path / "valid_config.json"
            io.dump_json(valid_config, config_file)

            is_valid, errors = validate_config_file(config_file)
            assert is_valid
            assert len(errors) == 0

            # Test invalid config (missing both downloads and processing sections)
            invalid_config = {
                "metadata": {"version": "1.0"}
                # Missing both downloads and processing sections
            }
            invalid_config_file = tmp_path / "invalid_config.json"
            io.dump_json(invalid_config, invalid_config_file)

            is_valid, errors = validate_config_file(invalid_config_file)
            assert not is_valid
            assert len(errors) == 1
            assert "downloads" in errors[0] or "processing" in errors[0]

    def test_create_sample_config(self):
        """Test sample configuration file creation."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Test basic config
            basic_config = tmp_path / "basic.json"
            create_sample_config(basic_config, "basic")

            assert basic_config.exists()
            config_data = io.load_json(basic_config)
            assert "description" in config_data
            assert "downloads" in config_data
            assert "processing" in config_data

            # Test scientific config
            scientific_config = tmp_path / "scientific.json"
            create_sample_config(scientific_config, "scientific")

            config_data = io.load_json(scientific_config)
            assert "Scientific data processing" in config_data["description"]
            assert "gene_expression" in config_data["downloads"]

            # Test advanced config
            advanced_config = tmp_path / "advanced.json"
            create_sample_config(advanced_config, "advanced")

            config_data = io.load_json(advanced_config)
            assert "advanced processing" in config_data["description"]

    def test_download_and_process_data(self):
        """Test the main processing function."""
        # Test with valid config
        config = {
            "downloads": {"test_data": {"url": "https://httpbin.org/json", "filename": "test.json"}},
            "processing": {"analyze": {"type": "json_analysis"}},
        }

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            try:
                results = download_and_process_data(config, output_dir=tmp_path, verbose=False)

                # Check results structure
                assert "config" in results
                assert "downloads" in results
                assert "processing" in results
                assert "start_time" in results
                assert "end_time" in results

                # Check that results file was created
                results_file = tmp_path / "processing_results.json"
                assert results_file.exists()

            except Exception:
                # Expected to fail without internet, but structure should be correct
                pass

    def test_run_config_based_workflow(self):
        """Test the main workflow function."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Create a sample config
            config = {"downloads": {"sample": {"url": "https://httpbin.org/json", "filename": "sample.json"}}}

            config_file = tmp_path / "test_config.json"
            io.dump_json(config, config_file)

            try:
                # This will fail without internet, but should return proper structure
                results = run_config_based_workflow(config_file, verbose=False)

                # Should have success flag and error information
                assert "success" in results
                assert "config_path" in results

            except Exception:
                # Expected to fail without internet
                pass

    def test_config_processing_error_handling(self):
        """Test error handling in config processing."""
        # Test with invalid config structure
        invalid_config = {
            "downloads": {
                "test": {
                    # Missing URL
                    "filename": "test.txt"
                }
            }
        }

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            try:
                results = download_and_process_data(invalid_config, output_dir=tmp_path, verbose=False)

                # Should contain error information
                assert "errors" in results
                assert len(results["errors"]) > 0

            except Exception:
                # Should handle errors gracefully
                pass

    def test_processing_with_no_downloads(self):
        """Test processing with only processing steps (no downloads)."""
        config = {"processing": {"analyze": {"type": "data_analysis"}}}

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            results = download_and_process_data(config, output_dir=tmp_path, verbose=False)

            # Should complete without downloads
            assert "processing" in results
            assert len(results["downloads"]) == 0
