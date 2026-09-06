"""Tests for core disk utilities.

Tests disk space monitoring and management functions following real-implementation policy.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

from metainformant.core.io import disk


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


class TestFreeSpaceAndSizes:
    """Tests for get_free_space, get_directory_size, and get_largest_files."""

    def test_get_free_space_returns_nonnegative_int(self, tmp_path):
        free = disk.get_free_space(tmp_path)
        assert isinstance(free, int)
        assert free >= 0

    def test_get_directory_size_sums_files_recursively(self, tmp_path):
        (tmp_path / "a.txt").write_bytes(b"x" * 100)
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "b.bin").write_bytes(b"y" * 50)
        assert disk.get_directory_size(tmp_path) == 150

    def test_get_directory_size_missing_dir_is_zero(self, tmp_path):
        assert disk.get_directory_size(tmp_path / "missing") == 0

    def test_get_largest_files_sorted_and_capped(self, tmp_path):
        (tmp_path / "small.txt").write_bytes(b"x")
        (tmp_path / "large.txt").write_bytes(b"x" * 100)
        (tmp_path / "medium.txt").write_bytes(b"x" * 10)
        top = disk.get_largest_files(tmp_path, n=2)
        assert [path.name for path, _ in top] == ["large.txt", "medium.txt"]
        assert top[0][1] == 100

    def test_get_largest_files_missing_dir_is_empty(self, tmp_path):
        assert disk.get_largest_files(tmp_path / "missing") == []


class TestCleanupFunctions:
    """Tests for cleanup_temp_files and cleanup_old_files."""

    def test_cleanup_temp_files_removes_only_old_files(self, tmp_path):
        old = tmp_path / "old.tmp"
        old.write_text("stale")
        old_time = time.time() - 48 * 3600
        os.utime(old, (old_time, old_time))
        fresh = tmp_path / "fresh.tmp"
        fresh.write_text("current")

        removed = disk.cleanup_temp_files(tmp_path, max_age_hours=24)

        assert removed == 1
        assert fresh.exists()

    def test_cleanup_temp_files_missing_dir_returns_zero(self, tmp_path):
        assert disk.cleanup_temp_files(tmp_path / "missing") == 0

    def test_cleanup_old_files_respects_exclude_patterns(self, tmp_path):
        old_time = time.time() - 60 * 24 * 3600
        keep = tmp_path / "important.keep"
        keep.write_text("keep me")
        drop = tmp_path / "scratch.csv"
        drop.write_text("drop me")
        for path in (keep, drop):
            os.utime(path, (old_time, old_time))

        removed = disk.cleanup_old_files(tmp_path, max_age_days=30, exclude_patterns=["*.keep"])

        assert removed == 1
        assert keep.exists()
        assert not drop.exists()


class TestMonitorAndSafeRemove:
    """Tests for monitor_disk_space, ensure_disk_space, and safe_remove_directory."""

    def test_monitor_disk_space_reports_ok_when_below_thresholds(self, tmp_path):
        result = disk.monitor_disk_space(tmp_path, warning_threshold=1.1, critical_threshold=1.2)
        assert result["status"] == "ok"
        assert result["usage"]["path"] == str(tmp_path)

    def test_monitor_disk_space_flags_warning_when_over_threshold(self, tmp_path):
        result = disk.monitor_disk_space(tmp_path, warning_threshold=-1.0, critical_threshold=2.0)
        assert result["status"] == "warning"

    def test_ensure_disk_space_sufficient_returns_true(self, tmp_path):
        assert disk.ensure_disk_space(tmp_path, 1) is True

    def test_ensure_disk_space_insufficient_raises(self, tmp_path):
        with pytest.raises(RuntimeError, match="Insufficient disk space"):
            disk.ensure_disk_space(tmp_path, 10**30)

    def test_safe_remove_directory_removes_tree(self, tmp_path):
        victim = tmp_path / "victim"
        (victim / "nested").mkdir(parents=True)
        (victim / "nested" / "file.txt").write_text("data")

        assert disk.safe_remove_directory(victim) is True
        assert not victim.exists()

    def test_safe_remove_directory_missing_dir_is_success(self, tmp_path):
        assert disk.safe_remove_directory(tmp_path / "ghost") is True

    def test_safe_remove_directory_refuses_oversized(self, tmp_path):
        victim = tmp_path / "big"
        victim.mkdir()
        (victim / "payload.bin").write_bytes(b"\0" * (2 * 1024 * 1024))

        assert disk.safe_remove_directory(victim, confirm_size_mb=1) is False
        assert victim.exists()
