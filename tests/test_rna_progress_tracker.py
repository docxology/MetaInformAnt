"""Tests for RNA progress tracker functionality.

This module tests progress tracking following NO_MOCKING_POLICY.
All tests use real file operations.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.rna.progress_tracker import ProgressTracker, get_tracker


class TestProgressTracker:
    """Test ProgressTracker class."""

    def test_progress_tracker_initialization(self, tmp_path: Path):
        """Test ProgressTracker initialization."""
        state_file = tmp_path / "progress_state.json"
        dashboard_file = tmp_path / "progress_dashboard.txt"
        
        tracker = ProgressTracker(
            state_file=state_file,
            dashboard_file=dashboard_file,
        )
        
        assert tracker.state_file == state_file
        assert tracker.dashboard_file == dashboard_file
        assert isinstance(tracker.state, dict)

    def test_initialize_species(self, tmp_path: Path):
        """Test initialize_species method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        sample_ids = ["SRR123", "SRR456", "SRR789"]
        tracker.initialize_species("test_species", 3, sample_ids)
        
        assert tracker.state["test_species"]["total"] == 3
        assert tracker.state["test_species"]["need_download"] == set(sample_ids)
        assert len(tracker.state["test_species"]["need_download"]) == 3

    def test_on_download_start(self, tmp_path: Path):
        """Test on_download_start method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        tracker.initialize_species("test_species", 1, ["SRR123"])
        tracker.on_download_start("test_species", "SRR123")
        
        assert "SRR123" in tracker.state["test_species"]["ongoing_download"]
        assert "SRR123" not in tracker.state["test_species"]["need_download"]

    def test_on_download_complete(self, tmp_path: Path):
        """Test on_download_complete method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        tracker.initialize_species("test_species", 1, ["SRR123"])
        tracker.on_download_start("test_species", "SRR123")
        tracker.on_download_complete("test_species", "SRR123")
        
        assert "SRR123" in tracker.state["test_species"]["needs_quant"]
        assert "SRR123" not in tracker.state["test_species"]["ongoing_download"]

    def test_on_quantification_complete(self, tmp_path: Path):
        """Test on_quantification_complete method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        tracker.initialize_species("test_species", 1, ["SRR123"])
        tracker.on_download_start("test_species", "SRR123")
        tracker.on_download_complete("test_species", "SRR123")
        tracker.on_quantification_complete("test_species", "SRR123")
        
        assert "SRR123" in tracker.state["test_species"]["needs_delete"]
        assert "SRR123" not in tracker.state["test_species"]["needs_quant"]

    def test_on_deletion_complete(self, tmp_path: Path):
        """Test on_deletion_complete method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        tracker.initialize_species("test_species", 1, ["SRR123"])
        tracker.on_download_start("test_species", "SRR123")
        tracker.on_download_complete("test_species", "SRR123")
        tracker.on_quantification_complete("test_species", "SRR123")
        tracker.on_deletion_complete("test_species", "SRR123")
        
        assert "SRR123" in tracker.state["test_species"]["completed"]
        assert "SRR123" not in tracker.state["test_species"]["needs_delete"]

    def test_state_persistence(self, tmp_path: Path):
        """Test that state is persisted to disk."""
        state_file = tmp_path / "progress_state.json"
        tracker1 = ProgressTracker(state_file=state_file)
        
        tracker1.initialize_species("test_species", 2, ["SRR123", "SRR456"])
        tracker1.on_download_start("test_species", "SRR123")
        
        # Create new tracker and load state
        tracker2 = ProgressTracker(state_file=state_file)
        
        assert tracker2.state["test_species"]["total"] == 2
        assert "SRR123" in tracker2.state["test_species"]["ongoing_download"]
        assert "SRR456" in tracker2.state["test_species"]["need_download"]

    def test_get_species_summary(self, tmp_path: Path):
        """Test get_species_summary method."""
        state_file = tmp_path / "progress_state.json"
        tracker = ProgressTracker(state_file=state_file)
        
        tracker.initialize_species("test_species", 3, ["SRR123", "SRR456", "SRR789"])
        tracker.on_download_start("test_species", "SRR123")
        tracker.on_download_complete("test_species", "SRR123")
        tracker.on_quantification_complete("test_species", "SRR123")
        tracker.on_deletion_complete("test_species", "SRR123")
        
        summary = tracker.get_species_summary("test_species")
        assert isinstance(summary, dict)
        assert summary["total"] == 3
        assert summary["completed"] == 1
        assert summary["need_download"] == 2


class TestGetTracker:
    """Test get_tracker function."""

    def test_get_tracker_default(self):
        """Test get_tracker with default parameters."""
        tracker = get_tracker()
        assert isinstance(tracker, ProgressTracker)

    def test_get_tracker_custom_paths(self, tmp_path: Path):
        """Test get_tracker with custom paths."""
        state_file = tmp_path / "custom_state.json"
        dashboard_file = tmp_path / "custom_dashboard.txt"
        
        tracker = get_tracker(state_file=state_file, dashboard_file=dashboard_file)
        assert tracker.state_file == state_file
        assert tracker.dashboard_file == dashboard_file


class TestProgressTrackerDocumentation:
    """Test that progress tracker has proper documentation."""

    def test_class_has_docstring(self):
        """Verify ProgressTracker class has docstring."""
        assert ProgressTracker.__doc__ is not None
        assert len(ProgressTracker.__doc__.strip()) > 0

    def test_get_tracker_has_docstring(self):
        """Verify get_tracker function has docstring."""
        assert get_tracker.__doc__ is not None
        assert len(get_tracker.__doc__.strip()) > 0

