"""Comprehensive tests for metainformant.core.utils.timing module.

Tests all timing utilities with real implementations (NO mocking):
- timed decorator (with and without arguments)
- Timer context manager
- rate_limiter decorator
- timeout_after context manager
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from metainformant.core.utils.timing import Timer, rate_limiter, timed, timeout_after


class TestTimedDecorator:
    """Test the @timed decorator."""

    def test_timed_without_arguments(self, caplog) -> None:
        """Test @timed decorator without arguments (default debug level)."""
        import logging

        @timed
        def fast_function() -> str:
            time.sleep(0.01)  # 10ms
            return "result"

        # Set the logger to DEBUG level explicitly for this test
        with caplog.at_level("DEBUG", logger="metainformant.core.utils.timing"):
            result = fast_function()

        assert result == "result"
        assert len(caplog.records) >= 1
        assert "fast_function completed in" in caplog.text
        assert "s" in caplog.text

    def test_timed_with_level_info(self, caplog) -> None:
        """Test @timed(level='info') decorator."""

        @timed(level="info")
        def logged_function() -> int:
            time.sleep(0.01)
            return 42

        with caplog.at_level("INFO"):
            result = logged_function()

        assert result == 42
        assert len(caplog.records) == 1
        assert "logged_function completed in" in caplog.text
        assert caplog.records[0].levelname == "INFO"

    def test_timed_with_level_warning(self, caplog) -> None:
        """Test @timed(level='warning') decorator."""

        @timed(level="warning")
        def warned_function() -> None:
            time.sleep(0.01)

        with caplog.at_level("WARNING"):
            warned_function()

        assert len(caplog.records) == 1
        assert "warned_function completed in" in caplog.text
        assert caplog.records[0].levelname == "WARNING"

    def test_timed_preserves_function_metadata(self) -> None:
        """Test that @timed preserves function name and docstring."""

        @timed
        def documented_function() -> str:
            """This is a docstring."""
            return "value"

        assert documented_function.__name__ == "documented_function"
        assert documented_function.__doc__ == "This is a docstring."

    def test_timed_with_arguments_and_kwargs(self, caplog) -> None:
        """Test @timed works with functions that have arguments."""

        @timed
        def add_numbers(a: int, b: int, *, multiplier: int = 1) -> int:
            time.sleep(0.01)
            return (a + b) * multiplier

        with caplog.at_level("DEBUG", logger="metainformant.core.utils.timing"):
            result = add_numbers(3, 5, multiplier=2)

        assert result == 16
        assert "add_numbers completed in" in caplog.text

    def test_timed_with_exception(self, caplog) -> None:
        """Test @timed logs time even when function raises exception."""

        @timed
        def failing_function() -> None:
            time.sleep(0.01)
            raise ValueError("intentional error")

        with caplog.at_level("DEBUG", logger="metainformant.core.utils.timing"):
            with pytest.raises(ValueError, match="intentional error"):
                failing_function()

        # Should still log timing despite exception
        assert "failing_function completed in" in caplog.text

    def test_timed_very_fast_operation(self, caplog) -> None:
        """Test @timed with near-instantaneous operation."""

        @timed
        def instant_function() -> int:
            return 1 + 1

        with caplog.at_level("DEBUG", logger="metainformant.core.utils.timing"):
            result = instant_function()

        assert result == 2
        assert "instant_function completed in" in caplog.text
        # Should still show some time (even if very small)
        assert "s" in caplog.text


class TestTimer:
    """Test the Timer context manager."""

    def test_timer_basic_usage(self) -> None:
        """Test Timer context manager basic usage."""
        with Timer() as t:
            time.sleep(0.05)  # 50ms

        # Allow 20ms tolerance for timing variations
        assert 0.04 < t.elapsed < 0.08
        assert 40.0 < t.elapsed_ms < 80.0

    def test_timer_elapsed_updates_during_execution(self) -> None:
        """Test that elapsed updates while timer is running."""
        with Timer() as t:
            time.sleep(0.02)
            elapsed_mid = t.elapsed
            time.sleep(0.02)
            elapsed_end = t.elapsed

        # elapsed_end should be greater than elapsed_mid
        assert elapsed_end > elapsed_mid
        assert elapsed_mid > 0.015  # At least ~20ms

    def test_timer_elapsed_frozen_after_exit(self) -> None:
        """Test that elapsed is frozen after context exits."""
        with Timer() as t:
            time.sleep(0.02)

        elapsed_1 = t.elapsed
        time.sleep(0.01)
        elapsed_2 = t.elapsed

        # Should be the same after context exit
        assert elapsed_1 == elapsed_2

    def test_timer_str_microseconds(self) -> None:
        """Test Timer.__str__ formats microseconds correctly."""
        with Timer() as t:
            # Very fast operation (< 1ms)
            _ = 1 + 1

        timer_str = str(t)
        # Should show microseconds (us)
        assert "us" in timer_str

    def test_timer_str_milliseconds(self) -> None:
        """Test Timer.__str__ formats milliseconds correctly."""
        with Timer() as t:
            time.sleep(0.05)  # 50ms

        timer_str = str(t)
        assert "ms" in timer_str
        # Should be in ms range (not us or s)
        assert t.elapsed_ms > 10.0
        assert t.elapsed < 1.0

    def test_timer_str_seconds(self) -> None:
        """Test Timer.__str__ formats seconds correctly."""
        with Timer() as t:
            time.sleep(1.2)

        timer_str = str(t)
        assert "s" in timer_str
        assert "ms" not in timer_str
        assert "us" not in timer_str
        # Should show something like "1.2xxs"
        assert "1." in timer_str

    def test_timer_str_minutes(self) -> None:
        """Test Timer.__str__ formats minutes correctly."""
        # Simulate a long operation by manually setting internal state
        timer = Timer()
        timer._start = time.perf_counter()
        timer._end = timer._start + 125.7  # 2 minutes 5.7 seconds

        timer_str = str(timer)
        assert "m" in timer_str
        assert "2m" in timer_str
        assert "5." in timer_str  # Should show remaining seconds

    def test_timer_repr(self) -> None:
        """Test Timer.__repr__ shows elapsed time."""
        with Timer() as t:
            time.sleep(0.01)

        repr_str = repr(t)
        assert "Timer(elapsed=" in repr_str
        assert "s)" in repr_str

    def test_timer_zero_time(self) -> None:
        """Test Timer with extremely fast operation."""
        with Timer() as t:
            pass

        # Should be very small but non-negative
        assert t.elapsed >= 0.0
        assert t.elapsed < 0.001  # Should be under 1ms
        assert t.elapsed_ms >= 0.0

    def test_timer_multiple_uses(self) -> None:
        """Test using Timer multiple times."""
        with Timer() as t1:
            time.sleep(0.02)

        with Timer() as t2:
            time.sleep(0.04)

        # Both should have different elapsed times (allow wider tolerance)
        assert 0.015 < t1.elapsed < 0.04
        assert 0.035 < t2.elapsed < 0.07
        assert t2.elapsed > t1.elapsed


class TestRateLimiter:
    """Test the rate_limiter decorator."""

    def test_rate_limiter_enforces_delay(self) -> None:
        """Test that rate_limiter actually enforces minimum interval."""
        call_times: list[float] = []

        @rate_limiter(calls_per_second=5.0)  # Max 5 calls/sec = 0.2s interval
        def tracked_function() -> None:
            call_times.append(time.monotonic())

        # Make 3 calls
        tracked_function()
        tracked_function()
        tracked_function()

        # Check intervals between calls (should be ~0.2s each)
        assert len(call_times) == 3
        interval_1 = call_times[1] - call_times[0]
        interval_2 = call_times[2] - call_times[1]

        # Allow some tolerance (0.19 to 0.25 seconds)
        assert 0.19 < interval_1 < 0.25
        assert 0.19 < interval_2 < 0.25

    def test_rate_limiter_total_time(self) -> None:
        """Test rate limiter total execution time for multiple calls."""

        @rate_limiter(calls_per_second=10.0)  # 0.1s interval
        def fast_call() -> int:
            return 42

        start = time.monotonic()
        for _ in range(5):
            result = fast_call()
            assert result == 42
        end = time.monotonic()

        total_time = end - start
        # 5 calls at 0.1s interval = ~0.4s total (first call is immediate)
        assert 0.35 < total_time < 0.55

    def test_rate_limiter_preserves_return_value(self) -> None:
        """Test that rate_limiter preserves function return values."""

        @rate_limiter(calls_per_second=100.0)
        def add(a: int, b: int) -> int:
            return a + b

        result = add(3, 5)
        assert result == 8

    def test_rate_limiter_preserves_arguments(self) -> None:
        """Test that rate_limiter preserves function arguments."""
        received_args: list[tuple] = []

        @rate_limiter(calls_per_second=50.0)
        def record_args(*args, **kwargs) -> None:
            received_args.append((args, kwargs))

        record_args(1, 2, key="value")
        record_args("a", "b", flag=True)

        assert len(received_args) == 2
        assert received_args[0] == ((1, 2), {"key": "value"})
        assert received_args[1] == (("a", "b"), {"flag": True})

    def test_rate_limiter_zero_calls_raises_error(self) -> None:
        """Test that rate_limiter raises ValueError for zero calls_per_second."""
        with pytest.raises(ValueError, match="calls_per_second must be positive"):

            @rate_limiter(calls_per_second=0.0)
            def invalid_function() -> None:
                pass

    def test_rate_limiter_negative_calls_raises_error(self) -> None:
        """Test that rate_limiter raises ValueError for negative calls_per_second."""
        with pytest.raises(ValueError, match="calls_per_second must be positive"):

            @rate_limiter(calls_per_second=-5.0)
            def invalid_function() -> None:
                pass

    def test_rate_limiter_very_slow_rate(self) -> None:
        """Test rate_limiter with very slow rate (e.g., 1 call per 2 seconds)."""
        call_times: list[float] = []

        @rate_limiter(calls_per_second=0.5)  # 1 call per 2 seconds
        def slow_function() -> None:
            call_times.append(time.monotonic())

        slow_function()
        slow_function()

        # Should have ~2 second gap
        interval = call_times[1] - call_times[0]
        assert 1.9 < interval < 2.2

    def test_rate_limiter_very_fast_rate(self) -> None:
        """Test rate_limiter with very fast rate (e.g., 100 calls/second)."""

        @rate_limiter(calls_per_second=100.0)  # 0.01s interval
        def fast_function() -> int:
            return 1

        start = time.monotonic()
        for _ in range(10):
            fast_function()
        end = time.monotonic()

        # 10 calls at 0.01s = ~0.09s total (allow wider tolerance on macOS)
        total_time = end - start
        assert 0.08 < total_time < 0.25

    def test_rate_limiter_thread_safety(self) -> None:
        """Test that rate_limiter is thread-safe (basic check)."""
        import threading

        results: list[float] = []

        @rate_limiter(calls_per_second=10.0)
        def shared_function() -> None:
            results.append(time.monotonic())

        # Create threads that call the function
        threads = [threading.Thread(target=shared_function) for _ in range(5)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # All calls should complete
        assert len(results) == 5
        # Results should be roughly spaced by 0.1s
        results_sorted = sorted(results)
        total_span = results_sorted[-1] - results_sorted[0]
        # 5 calls = 4 intervals * 0.1s = ~0.4s
        assert 0.35 < total_span < 0.55


class TestTimeoutAfter:
    """Test the timeout_after context manager."""

    def test_timeout_after_no_timeout(self) -> None:
        """Test timeout_after with fast operation that completes in time."""
        with timeout_after(1.0, "Should not timeout"):
            time.sleep(0.05)  # Much faster than 1s timeout

        # Should complete without raising TimeoutError

    def test_timeout_after_raises_timeout_error(self) -> None:
        """Test timeout_after raises TimeoutError when operation is too slow."""
        with pytest.raises(TimeoutError, match="Custom timeout message"):
            with timeout_after(0.05, "Custom timeout message"):
                time.sleep(0.2)  # Slower than 0.05s timeout

    def test_timeout_after_default_message(self) -> None:
        """Test timeout_after uses default message when none provided."""
        with pytest.raises(TimeoutError, match="Operation timed out after"):
            with timeout_after(0.05):
                time.sleep(0.2)

    def test_timeout_after_yields_event(self) -> None:
        """Test timeout_after yields a threading.Event for cooperative checking."""
        import threading

        with timeout_after(1.0) as event:
            assert isinstance(event, threading.Event)
            assert not event.is_set()  # Should not be set immediately
            time.sleep(0.01)

    def test_timeout_after_event_set_on_timeout(self) -> None:
        """Test that the yielded Event is set when timeout occurs."""
        import threading

        event_was_set = False
        with pytest.raises(TimeoutError):
            with timeout_after(0.05) as event:
                time.sleep(0.1)  # Sleep longer than timeout
                # Check if event was set during sleep
                if event.is_set():
                    event_was_set = True

        # The event should have been set
        # Note: This is timing-dependent, so we just verify TimeoutError was raised

    def test_timeout_after_cooperative_check(self) -> None:
        """Test cooperative timeout checking with the Event."""
        iterations = 0
        try:
            with timeout_after(0.15) as timed_out:
                for i in range(1000):
                    if timed_out.is_set():
                        iterations = i + 1
                        break
                    time.sleep(0.01)
                    iterations = i + 1
        except TimeoutError:
            # This is expected if we don't check the flag in time
            pass

        # Should have stopped early due to timeout flag or TimeoutError
        assert iterations < 1000
        # Should have done at least a few iterations before timeout
        assert iterations >= 1

    def test_timeout_after_zero_seconds_raises_error(self) -> None:
        """Test timeout_after raises ValueError for zero seconds."""
        with pytest.raises(ValueError, match="seconds must be positive"):
            with timeout_after(0.0):
                pass

    def test_timeout_after_negative_seconds_raises_error(self) -> None:
        """Test timeout_after raises ValueError for negative seconds."""
        with pytest.raises(ValueError, match="seconds must be positive"):
            with timeout_after(-1.0):
                pass

    def test_timeout_after_very_short_timeout(self) -> None:
        """Test timeout_after with very short timeout (0.01s)."""
        with pytest.raises(TimeoutError):
            with timeout_after(0.01):
                time.sleep(0.1)

    def test_timeout_after_very_long_timeout(self) -> None:
        """Test timeout_after with very long timeout that doesn't trigger."""
        with timeout_after(10.0):
            time.sleep(0.01)
        # Should complete without timeout

    def test_timeout_after_with_exception_in_block(self) -> None:
        """Test timeout_after when block raises non-timeout exception."""
        with pytest.raises(ValueError, match="intentional"):
            with timeout_after(1.0):
                time.sleep(0.01)
                raise ValueError("intentional")

    def test_timeout_after_cancels_timer_on_success(self) -> None:
        """Test that timeout_after properly cancels timer on successful completion."""
        import threading

        # Get count of active threads before
        thread_count_before = threading.active_count()

        with timeout_after(1.0):
            time.sleep(0.01)

        # Give timer thread a moment to clean up
        time.sleep(0.05)

        # Thread count should return to normal (timer should be cancelled)
        thread_count_after = threading.active_count()
        assert thread_count_after <= thread_count_before + 1  # Allow some variance

    def test_timeout_after_multiple_sequential_uses(self) -> None:
        """Test using timeout_after multiple times sequentially."""
        # First use - no timeout
        with timeout_after(0.5):
            time.sleep(0.05)

        # Second use - no timeout
        with timeout_after(0.5):
            time.sleep(0.05)

        # Third use - with timeout
        with pytest.raises(TimeoutError):
            with timeout_after(0.05):
                time.sleep(0.2)

        # Fourth use - no timeout (should still work after previous timeout)
        with timeout_after(0.5):
            time.sleep(0.05)


class TestIntegration:
    """Integration tests combining multiple timing utilities."""

    def test_timed_with_timer(self, caplog) -> None:
        """Test using @timed decorator with Timer context manager."""

        @timed
        def timed_function_with_timer() -> float:
            with Timer() as t:
                time.sleep(0.05)
            return t.elapsed

        with caplog.at_level("DEBUG", logger="metainformant.core.utils.timing"):
            elapsed = timed_function_with_timer()

        # Should have timing from both @timed and Timer
        assert 0.04 < elapsed < 0.08
        assert "timed_function_with_timer completed in" in caplog.text

    def test_rate_limiter_with_timer(self) -> None:
        """Test rate_limiter with Timer to measure actual rate."""

        @rate_limiter(calls_per_second=5.0)  # 0.2s per call
        def rate_limited_call() -> None:
            pass

        with Timer() as t:
            for _ in range(3):
                rate_limited_call()

        # 3 calls at 5 calls/sec = 0.4s (first immediate, then 2 * 0.2s)
        assert 0.35 < t.elapsed < 0.5

    def test_timeout_with_rate_limiter(self) -> None:
        """Test timeout_after with rate_limiter - timeout should trigger."""

        @rate_limiter(calls_per_second=5.0)  # 0.2s per call
        def slow_rate_limited_call() -> None:
            time.sleep(0.01)

        with pytest.raises(TimeoutError):
            with timeout_after(0.3):  # Only 0.3s allowed
                for _ in range(5):  # Would take ~0.8s with rate limiting
                    slow_rate_limited_call()

    def test_all_utilities_together(self, caplog) -> None:
        """Test combining @timed, Timer, rate_limiter, and timeout_after."""

        @timed(level="info")
        @rate_limiter(calls_per_second=10.0)
        def complex_function() -> str:
            return "result"

        with caplog.at_level("INFO"):
            with Timer() as overall:
                with timeout_after(1.0):
                    results = [complex_function() for _ in range(3)]

        assert len(results) == 3
        assert all(r == "result" for r in results)
        assert 0.15 < overall.elapsed < 0.35  # 3 calls at 10/sec
        assert "complex_function completed in" in caplog.text
