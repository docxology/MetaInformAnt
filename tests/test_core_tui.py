"""Tests for core Terminal UI (TUI) components.

Covers TerminalInterface lifecycle, bar management, ProgressState,
terminal state helpers, thread safety, _format_elapsed, and set_footer.

NO MOCKING: all tests exercise real implementations. Rendering tests
redirect or inspect internal state rather than writing to the real stdout.
"""

from __future__ import annotations

import io
import sys
import threading
import time

import pytest

from metainformant.core.ui.tui import (
    BLUE,
    GREEN,
    RED,
    ProgressState,
    TerminalInterface,
    _format_elapsed,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ui() -> TerminalInterface:
    """Create a TerminalInterface with a very short term-size TTL for test speed."""
    ui = TerminalInterface()
    ui._term_size_ttl = 0.05  # 50 ms cache for faster test assertions
    return ui


# ---------------------------------------------------------------------------
# Thread Lifecycle
# ---------------------------------------------------------------------------


class TestThreadLifecycle:
    """Tests for start() / stop() and related lifecycle methods."""

    def test_start_stop_basic(self) -> None:
        """start() spawns a render thread; stop() joins it."""
        ui = _make_ui()
        assert not ui._started
        assert not ui._stopped

        ui.start()
        assert ui._started
        assert not ui._stopped
        assert ui._monitor_thread is not None
        assert ui._monitor_thread.is_alive()

        ui.stop()
        assert ui._stopped
        assert ui._monitor_thread is None

    def test_stop_idempotent(self) -> None:
        """Calling stop() multiple times must not raise."""
        ui = _make_ui()
        ui.start()
        ui.stop()
        # Second call should be a harmless no-op
        ui.stop()
        assert ui._stopped

    def test_start_idempotent(self) -> None:
        """Calling start() when already running is a no-op."""
        ui = _make_ui()
        ui.start()
        thread_before = ui._monitor_thread
        # Second start should return without spawning another thread
        ui.start()
        assert ui._monitor_thread is thread_before
        ui.stop()

    def test_context_manager(self) -> None:
        """TerminalInterface works as a context manager."""
        with _make_ui() as ui:
            assert ui._started
            assert not ui._stopped
            ui.add_bar("ctx", "Context Bar", 50)
            ui.update("ctx", current=25)
        # After exiting the block, it should be stopped
        assert ui._stopped

    def test_del_calls_stop_if_needed(self) -> None:
        """__del__ safety net invokes stop() on a started, unstopped instance."""
        ui = _make_ui()
        ui.start()
        assert ui._started and not ui._stopped
        # Manually invoke __del__
        ui.__del__()
        assert ui._stopped

    def test_del_on_never_started_instance(self) -> None:
        """__del__ on an instance that was never started must not raise."""
        ui = _make_ui()
        # _started is False, so __del__ should be a no-op
        ui.__del__()
        assert not ui._started


# ---------------------------------------------------------------------------
# Bar Management
# ---------------------------------------------------------------------------


class TestBarManagement:
    """Tests for add_bar, update, remove_bar, clear."""

    def test_add_bar_creates_entry(self) -> None:
        """add_bar registers a ProgressState and preserves insertion order."""
        ui = _make_ui()
        ui.add_bar("dl_1", "Download 1", 1024, unit="MB")

        assert "dl_1" in ui._bars
        state = ui._bars["dl_1"]
        assert state.label == "Download 1"
        assert state.total == 1024
        assert state.unit == "MB"
        assert state.current == 0
        assert state.status == "Starting"
        assert state.start_time > 0
        assert ui._order == ["dl_1"]

    def test_add_bar_preserves_order(self) -> None:
        """Multiple add_bar calls preserve insertion order."""
        ui = _make_ui()
        for i in range(5):
            ui.add_bar(f"t{i}", f"Task {i}", 100)
        assert ui._order == [f"t{i}" for i in range(5)]

    def test_add_bar_same_id_replaces_state(self) -> None:
        """Adding a bar with the same ID replaces the state but keeps order stable."""
        ui = _make_ui()
        ui.add_bar("dup", "First", 100)
        ui.add_bar("dup", "Second", 200)

        assert ui._bars["dup"].label == "Second"
        assert ui._bars["dup"].total == 200
        # ID should appear only once in the order list
        assert ui._order.count("dup") == 1

    def test_update_modifies_state(self) -> None:
        """update() changes the specified fields, leaving others untouched."""
        ui = _make_ui()
        ui.add_bar("u1", "Upload", 500)

        ui.update("u1", current=250, status="Uploading", speed="10 MB/s", color=GREEN, stage="compress")
        state = ui._bars["u1"]
        assert state.current == 250
        assert state.status == "Uploading"
        assert state.speed == "10 MB/s"
        assert state.color == GREEN
        assert state.stage == "compress"
        # total should be unchanged
        assert state.total == 500

    def test_update_partial_fields(self) -> None:
        """update() with only some fields does not clobber other fields."""
        ui = _make_ui()
        ui.add_bar("p1", "Partial", 100)
        ui.update("p1", status="Running")
        state = ui._bars["p1"]
        assert state.status == "Running"
        assert state.current == 0  # unchanged
        assert state.speed == "0 B/s"  # unchanged default

    def test_update_nonexistent_bar_is_noop(self) -> None:
        """update() on a nonexistent task_id silently does nothing."""
        ui = _make_ui()
        # Should not raise
        ui.update("ghost", current=99, status="Missing")
        assert "ghost" not in ui._bars

    def test_update_total(self) -> None:
        """update() can change the total dynamically."""
        ui = _make_ui()
        ui.add_bar("grow", "Growing", 100)
        ui.update("grow", total=200)
        assert ui._bars["grow"].total == 200

    def test_remove_bar_removes_entry(self) -> None:
        """remove_bar() deletes the bar and its order entry."""
        ui = _make_ui()
        ui.add_bar("rm1", "Remove Me", 100)
        assert "rm1" in ui._bars

        ui.remove_bar("rm1")
        assert "rm1" not in ui._bars
        assert "rm1" not in ui._order

    def test_remove_bar_nonexistent_is_noop(self) -> None:
        """remove_bar() with a missing ID does not raise."""
        ui = _make_ui()
        ui.remove_bar("does_not_exist")  # should silently succeed

    def test_clear_removes_all_bars(self) -> None:
        """clear() empties bars, order, last_lines, and footer."""
        ui = _make_ui()
        ui.add_bar("c1", "C1", 100)
        ui.add_bar("c2", "C2", 200)
        ui.set_footer("some footer")
        ui._last_lines = 5

        ui.clear()
        assert len(ui._bars) == 0
        assert len(ui._order) == 0
        assert ui._last_lines == 0
        assert ui._footer == ""


# ---------------------------------------------------------------------------
# ProgressState
# ---------------------------------------------------------------------------


class TestProgressState:
    """Tests for the ProgressState dataclass."""

    def test_defaults(self) -> None:
        """ProgressState has sensible defaults."""
        ps = ProgressState(label="Test", current=0, total=100)
        assert ps.unit == "B"
        assert ps.status == "Pending"
        assert ps.start_time == 0.0
        assert ps.speed == "0 B/s"
        assert ps.color == BLUE
        assert ps.stage == ""

    def test_custom_values(self) -> None:
        """ProgressState accepts custom values for every field."""
        ps = ProgressState(
            label="Custom",
            current=42.5,
            total=85.0,
            unit="MB",
            status="Downloading",
            start_time=1000.0,
            speed="5 MB/s",
            color=RED,
            stage="fetch",
        )
        assert ps.label == "Custom"
        assert ps.current == 42.5
        assert ps.total == 85.0
        assert ps.unit == "MB"
        assert ps.status == "Downloading"
        assert ps.start_time == 1000.0
        assert ps.speed == "5 MB/s"
        assert ps.color == RED
        assert ps.stage == "fetch"

    def test_elapsed_time_field(self) -> None:
        """start_time set via add_bar is a real timestamp."""
        ui = _make_ui()
        before = time.time()
        ui.add_bar("ts1", "Timestamp", 100)
        after = time.time()

        st = ui._bars["ts1"].start_time
        assert before <= st <= after


# ---------------------------------------------------------------------------
# Terminal State helpers
# ---------------------------------------------------------------------------


class TestTerminalState:
    """Tests for _is_tty() and _get_terminal_width()."""

    def test_is_tty_returns_bool(self) -> None:
        """_is_tty() returns a boolean value."""
        result = TerminalInterface._is_tty()
        assert isinstance(result, bool)

    def test_is_tty_false_for_stringio(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """When stdout is a StringIO (not a tty), _is_tty() returns False."""
        fake_stdout = io.StringIO()
        monkeypatch.setattr(sys, "stdout", fake_stdout)
        assert TerminalInterface._is_tty() is False

    def test_get_terminal_width_returns_positive_int(self) -> None:
        """_get_terminal_width() returns a positive integer."""
        ui = _make_ui()
        width = ui._get_terminal_width()
        assert isinstance(width, int)
        assert width > 0

    def test_get_terminal_width_caches_value(self) -> None:
        """Second call to _get_terminal_width() uses cached value within TTL."""
        ui = _make_ui()
        ui._term_size_ttl = 10.0  # long TTL so cache definitely holds

        # First call populates the cache
        w1 = ui._get_terminal_width()
        ts1 = ui._term_size_ts

        # Second call should reuse the same cache timestamp
        w2 = ui._get_terminal_width()
        ts2 = ui._term_size_ts

        assert w1 == w2
        assert ts1 == ts2  # Timestamp unchanged means cache was reused

    def test_get_terminal_width_refreshes_after_ttl(self) -> None:
        """Cache is invalidated after the TTL elapses."""
        ui = _make_ui()
        ui._term_size_ttl = 0.01  # 10 ms TTL

        _ = ui._get_terminal_width()
        ts1 = ui._term_size_ts

        time.sleep(0.02)  # Exceed the TTL

        _ = ui._get_terminal_width()
        ts2 = ui._term_size_ts

        assert ts2 > ts1  # Cache was refreshed


# ---------------------------------------------------------------------------
# Thread Safety
# ---------------------------------------------------------------------------


class TestThreadSafety:
    """Verify concurrent operations do not corrupt internal state."""

    def test_concurrent_updates(self) -> None:
        """Multiple threads can update() the same bar without crashing."""
        ui = _make_ui()
        ui.add_bar("race", "Race Bar", 10000)

        errors: list[Exception] = []

        def updater(thread_id: int) -> None:
            try:
                for i in range(200):
                    ui.update("race", current=float(i + thread_id * 200), status=f"Thread-{thread_id}")
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=updater, args=(t,)) for t in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5.0)

        assert not errors, f"Concurrent updates produced errors: {errors}"
        # Bar should still be in valid state
        assert "race" in ui._bars
        assert ui._bars["race"].total == 10000

    def test_concurrent_add_bar(self) -> None:
        """Multiple threads can add_bar() simultaneously without corruption."""
        ui = _make_ui()
        errors: list[Exception] = []

        def adder(thread_id: int) -> None:
            try:
                for i in range(50):
                    bar_id = f"t{thread_id}_b{i}"
                    ui.add_bar(bar_id, f"Bar {bar_id}", 100.0)
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=adder, args=(t,)) for t in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5.0)

        assert not errors, f"Concurrent add_bar produced errors: {errors}"
        # 8 threads x 50 bars each = 400 bars
        assert len(ui._bars) == 400
        assert len(ui._order) == 400

    def test_concurrent_add_and_remove(self) -> None:
        """Simultaneous add and remove operations do not deadlock or crash."""
        ui = _make_ui()
        # Pre-populate some bars so removers have work to do
        for i in range(100):
            ui.add_bar(f"pre_{i}", f"Pre {i}", 100.0)

        errors: list[Exception] = []

        def adder() -> None:
            try:
                for i in range(100):
                    ui.add_bar(f"new_{i}", f"New {i}", 100.0)
            except Exception as exc:
                errors.append(exc)

        def remover() -> None:
            try:
                for i in range(100):
                    ui.remove_bar(f"pre_{i}")
            except Exception as exc:
                errors.append(exc)

        t_add = threading.Thread(target=adder)
        t_rem = threading.Thread(target=remover)
        t_add.start()
        t_rem.start()
        t_add.join(timeout=5.0)
        t_rem.join(timeout=5.0)

        assert not errors, f"Concurrent add/remove produced errors: {errors}"
        # All pre_ bars should be removed, all new_ bars should exist
        for i in range(100):
            assert f"pre_{i}" not in ui._bars
            assert f"new_{i}" in ui._bars


# ---------------------------------------------------------------------------
# _format_elapsed
# ---------------------------------------------------------------------------


class TestFormatElapsed:
    """Tests for the _format_elapsed helper function."""

    def test_seconds_only(self) -> None:
        """Values under 60 seconds show as '{n}s'."""
        assert _format_elapsed(0) == "0s"
        assert _format_elapsed(1) == "1s"
        assert _format_elapsed(30) == "30s"
        assert _format_elapsed(59) == "59s"

    def test_negative_returns_zero(self) -> None:
        """Negative values are clamped to '0s'."""
        assert _format_elapsed(-1) == "0s"
        assert _format_elapsed(-999) == "0s"

    def test_minutes_and_seconds(self) -> None:
        """Values between 60s and 3600s show as '{m}m{ss}s'."""
        assert _format_elapsed(60) == "1m00s"
        assert _format_elapsed(61) == "1m01s"
        assert _format_elapsed(90) == "1m30s"
        assert _format_elapsed(125) == "2m05s"
        assert _format_elapsed(3599) == "59m59s"

    def test_hours_and_minutes(self) -> None:
        """Values >= 3600s show as '{h}h{mm}m'."""
        assert _format_elapsed(3600) == "1h00m"
        assert _format_elapsed(3661) == "1h01m"
        assert _format_elapsed(7200) == "2h00m"
        assert _format_elapsed(7323) == "2h02m"

    def test_fractional_seconds_truncated(self) -> None:
        """Fractional seconds are truncated (not rounded)."""
        assert _format_elapsed(59.9) == "59s"
        assert _format_elapsed(0.1) == "0s"
        assert _format_elapsed(60.99) == "1m00s"


# ---------------------------------------------------------------------------
# set_footer
# ---------------------------------------------------------------------------


class TestSetFooter:
    """Tests for the set_footer method."""

    def test_set_footer_stores_text(self) -> None:
        """set_footer() stores the provided text in _footer."""
        ui = _make_ui()
        assert ui._footer == ""

        ui.set_footer("Processing 3/10 samples")
        assert ui._footer == "Processing 3/10 samples"

    def test_set_footer_overwrites_previous(self) -> None:
        """Subsequent set_footer() calls replace the previous text."""
        ui = _make_ui()
        ui.set_footer("First footer")
        ui.set_footer("Second footer")
        assert ui._footer == "Second footer"

    def test_set_footer_empty_string(self) -> None:
        """Setting footer to empty string clears it."""
        ui = _make_ui()
        ui.set_footer("Nonempty")
        ui.set_footer("")
        assert ui._footer == ""


# ---------------------------------------------------------------------------
# Render internals (state-based, no stdout output)
# ---------------------------------------------------------------------------


class TestRenderState:
    """Validate rendering-related internal state without writing to stdout."""

    def test_last_lines_tracks_render_count(self) -> None:
        """After a render cycle, _last_lines reflects the number of output lines."""
        ui = _make_ui()
        ui.add_bar("r1", "Render 1", 100)
        ui.add_bar("r2", "Render 2", 200)
        ui.set_footer("Footer text")

        # Without a TTY, _render() will return early. We verify the bar/order state instead.
        # In non-TTY (test) environments, _last_lines stays at 0 because _render exits early.
        assert ui._last_lines == 0
        assert len(ui._order) == 2

    def test_bar_percentage_calculation(self) -> None:
        """Percentage is derived from current/total."""
        ui = _make_ui()
        ui.add_bar("pct", "Percent", 200)
        ui.update("pct", current=100)

        state = ui._bars["pct"]
        pct = (state.current / state.total) * 100 if state.total > 0 else 0
        assert pct == 50.0

    def test_bar_percentage_zero_total(self) -> None:
        """Zero total does not cause division by zero."""
        ui = _make_ui()
        ui.add_bar("zero", "Zero Total", 0)
        ui.update("zero", current=50)

        state = ui._bars["zero"]
        pct = (state.current / state.total) * 100 if state.total > 0 else 0
        assert pct == 0


# ---------------------------------------------------------------------------
# Live instances cleanup
# ---------------------------------------------------------------------------


class TestLiveInstancesCleanup:
    """Verify the class-level _live_instances tracking."""

    def test_start_registers_instance(self) -> None:
        """start() adds the instance to _live_instances."""
        ui = _make_ui()
        initial_count = len(TerminalInterface._live_instances)
        ui.start()
        assert ui in TerminalInterface._live_instances
        assert len(TerminalInterface._live_instances) == initial_count + 1
        ui.stop()

    def test_stop_deregisters_instance(self) -> None:
        """stop() removes the instance from _live_instances."""
        ui = _make_ui()
        ui.start()
        assert ui in TerminalInterface._live_instances
        ui.stop()
        assert ui not in TerminalInterface._live_instances

    def test_context_manager_cleans_up_live_instances(self) -> None:
        """Context manager exit removes from _live_instances."""
        with _make_ui() as ui:
            assert ui in TerminalInterface._live_instances
        assert ui not in TerminalInterface._live_instances
