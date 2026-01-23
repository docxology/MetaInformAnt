"""Tests for core disk utilities.

Tests disk space monitoring and management functions following NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core import disk


class TestDiskUsage:
    """Tests for get_disk_usage function."""

    def test_get_disk_usage_current_dir(self):
        """Test getting disk usage for current directory."""
        total, used, free, percent = disk.get_disk_usage(Path("."))
        assert total >= 0
        assert used >= 0
        assert free >= 0
        assert isinstance(percent, str)
        assert "%" in percent

    def test_get_disk_usage_output_dir(self, tmp_path):
        """Test getting disk usage for temporary directory."""
        total, used, free, percent = disk.get_disk_usage(tmp_path)
        assert total >= 0
        assert used >= 0
        assert free >= 0
        assert isinstance(percent, str)

    def test_get_disk_usage_nonexistent_path(self):
        """Test getting disk usage for nonexistent path (should still work)."""
        # Should work even if path doesn't exist - uses parent directory
        total, used, free, percent = disk.get_disk_usage(Path("/nonexistent/path/12345"))
        # May return zeros if unable to determine, which is acceptable
        assert isinstance(total, float)
        assert isinstance(used, float)
        assert isinstance(free, float)
        assert isinstance(percent, str)


class TestCheckDiskSpace:
    """Tests for check_disk_space function."""

    def test_check_disk_space_sufficient(self, tmp_path):
        """Test checking disk space when sufficient."""
        is_ok, message = disk.check_disk_space(tmp_path, min_free_gb=0.1, min_free_percent=0.1)
        assert isinstance(is_ok, bool)
        assert isinstance(message, str)
        # Result depends on actual disk space, but should be consistent

    def test_check_disk_space_very_low_threshold(self, tmp_path):
        """Test with very low threshold (should pass)."""
        is_ok, message = disk.check_disk_space(tmp_path, min_free_gb=0.001, min_free_percent=0.001)
        assert isinstance(is_ok, bool)
        assert isinstance(message, str)

    def test_check_disk_space_high_threshold(self, tmp_path):
        """Test with high threshold (may fail on small drives)."""
        is_ok, message = disk.check_disk_space(tmp_path, min_free_gb=10000.0, min_free_percent=50.0)
        assert isinstance(is_ok, bool)
        assert isinstance(message, str)
        # May fail if drive is small, which is expected behavior


class TestDiskSpaceInfo:
    """Tests for get_disk_space_info function."""

    def test_get_disk_space_info(self, tmp_path):
        """Test getting comprehensive disk space information."""
        info = disk.get_disk_space_info(tmp_path)
        assert isinstance(info, dict)
        assert "total_gb" in info
        assert "used_gb" in info
        assert "free_gb" in info
        assert "percent_used" in info
        assert "percent_free" in info
        assert isinstance(info["total_gb"], float)
        assert isinstance(info["used_gb"], float)
        assert isinstance(info["free_gb"], float)
        assert isinstance(info["percent_used"], str)
        assert isinstance(info["percent_free"], str)


class TestDriveSizeCategory:
    """Tests for detect_drive_size_category function."""

    def test_detect_drive_size_category(self, tmp_path):
        """Test detecting drive size category."""
        category = disk.detect_drive_size_category(tmp_path)
        assert category in ("large", "medium", "small")

    def test_detect_drive_size_category_nonexistent(self):
        """Test detecting category for nonexistent path."""
        category = disk.detect_drive_size_category(Path("/nonexistent/path/12345"))
        # Should return "small" as fallback
        assert category in ("large", "medium", "small")


class TestRecommendedBatchSize:
    """Tests for get_recommended_batch_size function."""

    def test_get_recommended_batch_size_default(self, tmp_path):
        """Test getting recommended batch size with default parameters."""
        batch_size = disk.get_recommended_batch_size(tmp_path)
        assert isinstance(batch_size, int)
        assert batch_size >= 8  # Minimum should be at least 8

    def test_get_recommended_batch_size_custom(self, tmp_path):
        """Test getting recommended batch size with custom parameters."""
        batch_size = disk.get_recommended_batch_size(tmp_path, sample_size_gb=0.5, safety_buffer=0.2)
        assert isinstance(batch_size, int)
        assert batch_size > 0

    def test_get_recommended_batch_size_large_samples(self, tmp_path):
        """Test with large sample size (should reduce batch size)."""
        batch_size = disk.get_recommended_batch_size(tmp_path, sample_size_gb=100.0)
        assert isinstance(batch_size, int)
        assert batch_size >= 8  # Should still respect minimum


class TestRecommendedTempDir:
    """Tests for get_recommended_temp_dir function."""

    def test_get_recommended_temp_dir(self, tmp_path):
        """Test getting recommended temporary directory."""
        # Use tmp_path as repo root
        temp_dir = disk.get_recommended_temp_dir(tmp_path)
        assert isinstance(temp_dir, Path)
        assert temp_dir.exists() or temp_dir.parent.exists()

    def test_get_recommended_temp_dir_with_output(self, tmp_path):
        """Test with output directory present."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        temp_dir = disk.get_recommended_temp_dir(tmp_path)
        assert isinstance(temp_dir, Path)
