"""Graceful stop of a writer that runs worker threads and child process groups.

A producer that dies on SIGTERM (or leaves a ``ThreadPoolExecutor`` on a
second Ctrl-C) releases the campaign writer lock while its download, extraction
and ``kallisto`` process groups, which run in their own sessions, keep writing
the same sample directories. A relaunch then works on the same samples
concurrently. :class:`ProducerStop` replaces the default SIGINT, SIGTERM and
SIGHUP dispositions for the duration of a run:

* the handler only records the request (it never raises ``KeyboardInterrupt``
  and never takes a lock, so it cannot deadlock the main thread);
* the scheduling loop polls :attr:`ProducerStop.requested`, stops admitting
  work, cancels queued futures and calls :meth:`ProducerStop.stop_children`,
  which refuses new process groups and sends ``SIGTERM`` to the running ones,
  escalating to ``SIGKILL`` after ``grace_seconds`` or on a second signal;
* the writer releases its lock only after every worker thread has returned
  and every registered group has been reaped.

Handlers are installed only in the main thread; elsewhere (tests calling the
producer from a worker thread) the object is inert but still usable through
:meth:`request`.
"""

from __future__ import annotations

import os
import signal
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any

from metainformant.core.utils.logging import get_logger
from metainformant.rna.retrieval import process_groups

logger = get_logger(__name__)

STOP_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
DEFAULT_STOP_GRACE_SECONDS = 30.0


def stop_grace_seconds() -> float:
    """``AMALGKIT_STOP_GRACE_SECONDS`` (default 30): SIGTERM-to-SIGKILL grace for child process groups."""

    raw = os.environ.get("AMALGKIT_STOP_GRACE_SECONDS", "").strip()
    try:
        value = float(raw) if raw else DEFAULT_STOP_GRACE_SECONDS
    except ValueError:
        logger.warning("Invalid AMALGKIT_STOP_GRACE_SECONDS=%r; using %.0f", raw, DEFAULT_STOP_GRACE_SECONDS)
        return DEFAULT_STOP_GRACE_SECONDS
    return max(0.0, value)


class ProducerStop:
    """A stop request (from a signal or a caller) and the escalation of child process groups."""

    def __init__(self, grace_seconds: float | None = None) -> None:
        self.grace_seconds = stop_grace_seconds() if grace_seconds is None else max(0.0, grace_seconds)
        self._signals: list[int] = []
        self._requested_at: float | None = None
        self._terminated_at: float | None = None
        self._killed = False

    # -- request -----------------------------------------------------------------

    def _handle(self, signum: int, _frame: Any) -> None:
        # Async-signal context: record only (list append and float stores are atomic under the GIL).
        self._signals.append(signum)
        if self._requested_at is None:
            self._requested_at = time.monotonic()

    def request(self, signum: int = signal.SIGTERM) -> None:
        """Ask the run to stop as if ``signum`` had been received (``0``: an internal error, not a signal)."""

        self._handle(signum, None)

    @property
    def requested(self) -> bool:
        return self._requested_at is not None

    @property
    def signal_count(self) -> int:
        return len(self._signals)

    @property
    def signal_name(self) -> str | None:
        if not self._signals:
            return None
        if self._signals[0] == 0:
            return "internal"
        try:
            return signal.Signals(self._signals[0]).name
        except ValueError:
            return str(self._signals[0])

    @contextmanager
    def installed(self) -> Iterator["ProducerStop"]:
        """Replace the stop signals' handlers while the block runs (main thread only)."""

        if threading.current_thread() is not threading.main_thread():
            yield self
            return
        previous = {}
        try:
            for signum in STOP_SIGNALS:
                previous[signum] = signal.signal(signum, self._handle)
            yield self
        finally:
            for signum, handler in previous.items():
                signal.signal(signum, handler)

    # -- children ----------------------------------------------------------------

    def stop_children(self) -> None:
        """Refuse new process groups, TERM the running ones, and KILL them after the grace or a second signal.

        Call it repeatedly from the scheduling loop once :attr:`requested`;
        each call escalates as needed and is cheap otherwise.
        """

        now = time.monotonic()
        if process_groups.begin_shutdown(f"stop requested by {self.signal_name or 'caller'}"):
            self._terminated_at = now
            groups = process_groups.signal_all(signal.SIGTERM)
            logger.warning(
                "Stop requested (%s): no new samples are admitted; sent SIGTERM to %d child process group(s); "
                "the writer lock is held until every worker and child has exited",
                self.signal_name or "caller",
                len(groups),
            )
            return
        terminated_at = self._terminated_at if self._terminated_at is not None else now
        if self.signal_count >= 2 or now - terminated_at >= self.grace_seconds:
            groups = process_groups.signal_all(signal.SIGKILL)
            if groups and not self._killed:
                logger.warning("Sent SIGKILL to %d child process group(s) still running", len(groups))
            self._killed = self._killed or bool(groups)
