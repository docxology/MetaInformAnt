"""Timing and performance utilities for METAINFORMANT."""

from __future__ import annotations

import functools
import threading
import time
from collections.abc import Callable
from contextlib import contextmanager
from typing import Any, ParamSpec, TypeVar

from metainformant.core.utils.logging import get_logger

P = ParamSpec("P")
T = TypeVar("T")

logger = get_logger(__name__)


def timed(
    func: Callable[P, T] | None = None, *, level: str = "debug"
) -> Callable[P, T] | Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator that logs execution time of a function.

    Can be used with or without arguments:
        @timed
        def my_func(): ...

        @timed(level="info")
        def my_func(): ...

    Args:
        func: The function to wrap (when used without parentheses).
        level: Log level for the timing message. Default "debug".

    Returns:
        Wrapped function that logs its execution time.
    """

    def decorator(fn: Callable[P, T]) -> Callable[P, T]:
        @functools.wraps(fn)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            start = time.perf_counter()
            try:
                result = fn(*args, **kwargs)
                return result
            finally:
                elapsed = time.perf_counter() - start
                log_fn = getattr(logger, level.lower(), logger.debug)
                log_fn("%s completed in %.4fs", fn.__qualname__, elapsed)

        return wrapper

    if func is not None:
        return decorator(func)
    return decorator


class Timer:
    """Context manager for timing code blocks.

    Usage:
        with Timer() as t:
            do_something()
        print(t)  # "1.234s" or "56.7ms"
        print(t.elapsed)  # seconds as float
        print(t.elapsed_ms)  # milliseconds as float
    """

    def __init__(self) -> None:
        self._start: float = 0.0
        self._end: float | None = None

    def __enter__(self) -> Timer:
        self._start = time.perf_counter()
        self._end = None
        return self

    def __exit__(self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: Any) -> None:
        self._end = time.perf_counter()

    @property
    def elapsed(self) -> float:
        """Elapsed time in seconds."""
        end = self._end if self._end is not None else time.perf_counter()
        return end - self._start

    @property
    def elapsed_ms(self) -> float:
        """Elapsed time in milliseconds."""
        return self.elapsed * 1000.0

    def __str__(self) -> str:
        """Human-readable duration string."""
        secs = self.elapsed
        if secs < 0.001:
            return f"{secs * 1_000_000:.1f}us"
        if secs < 1.0:
            return f"{secs * 1000:.1f}ms"
        if secs < 60.0:
            return f"{secs:.3f}s"
        minutes = int(secs // 60)
        remaining = secs % 60
        return f"{minutes}m {remaining:.1f}s"

    def __repr__(self) -> str:
        return f"Timer(elapsed={self.elapsed:.6f}s)"


def rate_limiter(calls_per_second: float) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator factory that enforces minimum interval between calls.

    Useful for API rate limiting. Blocks the calling thread until the minimum
    interval has elapsed since the last call.

    Args:
        calls_per_second: Maximum number of calls allowed per second.

    Returns:
        Decorator that rate-limits the wrapped function.

    Example:
        @rate_limiter(calls_per_second=2.0)
        def fetch_api(url: str) -> dict: ...
    """
    if calls_per_second <= 0:
        raise ValueError(f"calls_per_second must be positive, got {calls_per_second}")

    min_interval = 1.0 / calls_per_second
    lock = threading.Lock()
    last_call_time: list[float] = [0.0]  # mutable container for closure

    def decorator(func: Callable[P, T]) -> Callable[P, T]:
        @functools.wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            with lock:
                now = time.monotonic()
                elapsed_since_last = now - last_call_time[0]
                if elapsed_since_last < min_interval:
                    sleep_time = min_interval - elapsed_since_last
                    time.sleep(sleep_time)
                last_call_time[0] = time.monotonic()
            return func(*args, **kwargs)

        return wrapper

    return decorator


@contextmanager
def timeout_after(seconds: float, message: str = ""):
    """Context manager that raises TimeoutError after N seconds.

    Uses threading.Timer (not signals) so it works in any thread and on any OS.
    Note: This sets a flag and checks it; it cannot interrupt blocking I/O.
    The TimeoutError is raised when the context block completes if the timer
    has already fired, or the timer callback sets a flag that cooperative code
    can check.

    For truly preemptive timeout on blocking operations, consider using
    concurrent.futures with a timeout instead.

    Args:
        seconds: Maximum time in seconds before timeout.
        message: Optional message for the TimeoutError.

    Raises:
        TimeoutError: If the code block exceeds the specified duration.
        ValueError: If seconds is not positive.

    Example:
        with timeout_after(5.0, "API call took too long"):
            response = slow_api_call()
    """
    if seconds <= 0:
        raise ValueError(f"seconds must be positive, got {seconds}")

    timed_out = threading.Event()
    error_message = message or f"Operation timed out after {seconds:.1f}s"

    def _on_timeout() -> None:
        timed_out.set()

    timer = threading.Timer(seconds, _on_timeout)
    timer.daemon = True
    timer.start()
    try:
        yield timed_out
        if timed_out.is_set():
            raise TimeoutError(error_message)
    finally:
        timer.cancel()
