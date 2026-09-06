"""Tests for core.utils.optional_deps and core.utils.progress.

Real implementations only: exercises the actual warning-state machine and the
real progress logging helpers (no tqdm mocked; fallback shim used when tqdm is
absent).
"""

from __future__ import annotations

import logging

import pytest

from metainformant.core.utils import optional_deps, progress


@pytest.fixture(autouse=True)
def _clean_warning_state():
    """Isolate the process-global warning state for each test."""

    optional_deps.reset_warning_state()
    yield
    optional_deps.reset_warning_state()


class TestOptionalDependencyWarnings:
    """Tests for the optional dependency warning state machine."""

    def test_warning_issued_once_per_module(self, caplog) -> None:
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("seaborn", "enhanced plots")
            optional_deps.warn_optional_dependency("seaborn", "enhanced plots")

        assert len(caplog.records) == 1
        assert "seaborn not available" in caplog.text

    def test_distinct_functionality_warns_independently(self, caplog) -> None:
        with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
            optional_deps.warn_optional_dependency("anndata", "single-cell IO")
            optional_deps.warn_optional_dependency("anndata", "single-cell plots")

        assert len(caplog.records) == 2

    def test_suppressed_warnings_emit_nothing(self, caplog) -> None:
        optional_deps.suppress_optional_warnings()
        try:
            with caplog.at_level(logging.WARNING, logger="metainformant.core.optional_deps"):
                optional_deps.warn_optional_dependency("torch", "deep learning")
        finally:
            optional_deps.enable_optional_warnings()

        assert caplog.records == []
        assert optional_deps.get_warning_state()["warnings_issued"] == set()

    def test_state_reflects_issued_warnings(self) -> None:
        optional_deps.warn_optional_dependency("pod5", "read processing")
        state = optional_deps.get_warning_state()
        assert state["suppress_warnings"] is False
        assert "pod5:read processing" in state["warnings_issued"]


class TestProgressHelpers:
    """Tests for log_progress and task_context."""

    def test_log_progress_with_total(self, caplog) -> None:
        with caplog.at_level(logging.INFO, logger="metainformant.core.utils.progress"):
            progress.log_progress(5, 10, "Processing")

        assert "5/10 (50.0%)" in caplog.text

    def test_log_progress_zero_total_is_percent_zero(self, caplog) -> None:
        """total=0 is a known total, not an indeterminate one."""

        with caplog.at_level(logging.INFO, logger="metainformant.core.utils.progress"):
            progress.log_progress(0, 0, "Empty batch")

        assert "0/0 (0.0%)" in caplog.text

    def test_log_progress_without_total(self, caplog) -> None:
        with caplog.at_level(logging.INFO, logger="metainformant.core.utils.progress"):
            progress.log_progress(3, None, "Streaming")

        assert "3 items" in caplog.text

    def test_task_context_logs_completion(self, caplog) -> None:
        with caplog.at_level(logging.INFO, logger="metainformant.core.utils.progress"):
            with progress.task_context("indexing", total_steps=2) as tracker:
                tracker.update()
                tracker.update(1)

        assert "Starting task: indexing" in caplog.text
        assert "Completed task: indexing" in caplog.text

    def test_task_context_logs_failure_and_reraises(self, caplog) -> None:
        with pytest.raises(RuntimeError, match="boom"):
            with caplog.at_level(logging.INFO, logger="metainformant.core.utils.progress"):
                with progress.task_context("failing-task"):
                    raise RuntimeError("boom")

        assert "Task failed: failing-task" in caplog.text
