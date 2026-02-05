"""Core Terminal UI components for MetaInformAnt.

This module provides a robust, dependency-free (standard library only)
Terminal User Interface for tracking parallel processes.
It uses ANSI escape codes for "oldschool" elegance.
"""

import atexit
import shutil
import signal
import sys
import threading
import time
from dataclasses import dataclass
from types import FrameType
from typing import ClassVar, Dict, List, Optional

# ANSI Escape Codes
CSI = "\033["
RESET = f"{CSI}0m"
BOLD = f"{CSI}1m"

# Colors
RED = f"{CSI}31m"
GREEN = f"{CSI}32m"
YELLOW = f"{CSI}33m"
BLUE = f"{CSI}34m"
MAGENTA = f"{CSI}35m"
CYAN = f"{CSI}36m"
WHITE = f"{CSI}37m"

# Cursor Controls
HIDE_CURSOR = f"{CSI}?25l"
SHOW_CURSOR = f"{CSI}?25h"
CLEAR_LINE = f"{CSI}2K"
MOVE_UP = f"{CSI}A"


@dataclass
class ProgressState:
    """State of a single progress bar."""

    label: str
    current: float
    total: float
    unit: str = "B"
    status: str = "Pending"
    start_time: float = 0.0
    speed: str = "0 B/s"
    color: str = BLUE
    stage: str = ""  # For multi-stage workflows (e.g., "Download -> Quant")


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds into a compact human-readable string."""
    if seconds < 0:
        return "0s"
    total_secs = int(seconds)
    if total_secs < 60:
        return f"{total_secs}s"
    minutes = total_secs // 60
    secs = total_secs % 60
    if minutes < 60:
        return f"{minutes}m{secs:02d}s"
    hours = minutes // 60
    mins = minutes % 60
    return f"{hours}h{mins:02d}m"


class TerminalInterface:
    """Manages terminal output for parallel tasks.

    Supports context-manager usage::

        with TerminalInterface() as ui:
            ui.add_bar("t1", "Task 1", 100)
            ui.update("t1", current=50)

    Or traditional start/stop::

        ui = TerminalInterface()
        ui.start()
        ...
        ui.stop()
    """

    # Track all live instances for atexit cleanup
    _live_instances: ClassVar[List["TerminalInterface"]] = []

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._bars: Dict[str, ProgressState] = {}
        self._order: List[str] = []
        self._stop_event = threading.Event()
        self._last_lines = 0
        self._footer = ""
        self._monitor_thread: Optional[threading.Thread] = None
        self._started = False
        self._stopped = False

        # Cached terminal size with TTL
        self._cached_term_size: Optional[tuple[int, int]] = None
        self._term_size_ts: float = 0.0
        self._term_size_ttl: float = 1.0  # seconds

        # Stash original SIGINT handler so we can restore it
        self._original_sigint: Optional[signal.Handlers] = None

    # ------------------------------------------------------------------
    # Terminal size caching
    # ------------------------------------------------------------------

    def _get_terminal_width(self) -> int:
        """Return terminal column count, cached for up to ``_term_size_ttl`` seconds."""
        now = time.monotonic()
        if self._cached_term_size is None or (now - self._term_size_ts) >= self._term_size_ttl:
            size = shutil.get_terminal_size()
            self._cached_term_size = (size.columns, size.lines)
            self._term_size_ts = now
        return self._cached_term_size[0]

    # ------------------------------------------------------------------
    # Terminal state save / restore
    # ------------------------------------------------------------------

    @staticmethod
    def _is_tty() -> bool:
        """Return True if stdout is connected to a terminal."""
        try:
            return sys.stdout.isatty()
        except Exception:
            return False

    def _save_terminal_state(self) -> None:
        """Hide cursor and register cleanup handlers."""
        if self._is_tty():
            sys.stdout.write(HIDE_CURSOR)
            sys.stdout.flush()

        # Register atexit only once per instance
        TerminalInterface._live_instances.append(self)

        # Install a SIGINT handler that restores the terminal
        try:
            self._original_sigint = signal.getsignal(signal.SIGINT)
            signal.signal(signal.SIGINT, self._sigint_handler)
        except (OSError, ValueError):
            # Cannot set signal handler from non-main thread -- that is fine
            pass

    def _restore_terminal_state(self) -> None:
        """Show cursor and remove cleanup registrations."""
        if self._is_tty():
            try:
                sys.stdout.write(SHOW_CURSOR)
                sys.stdout.write("\n")
                sys.stdout.flush()
            except (OSError, ValueError):
                pass

        # Remove from live instances
        try:
            TerminalInterface._live_instances.remove(self)
        except ValueError:
            pass

        # Restore original SIGINT handler
        if self._original_sigint is not None:
            try:
                signal.signal(signal.SIGINT, self._original_sigint)
            except (OSError, ValueError):
                pass
            self._original_sigint = None

    def _sigint_handler(self, signum: int, frame: Optional[FrameType]) -> None:
        """Handle Ctrl+C: restore terminal, then re-raise."""
        self.stop()
        # Re-raise as KeyboardInterrupt so callers can catch it normally
        raise KeyboardInterrupt

    # ------------------------------------------------------------------
    # Bar management
    # ------------------------------------------------------------------

    def add_bar(self, task_id: str, label: str, total: float, unit: str = "B") -> None:
        """Register a new progress bar."""
        with self._lock:
            if task_id not in self._bars:
                self._order.append(task_id)
            self._bars[task_id] = ProgressState(
                label=label, current=0, total=total, unit=unit, status="Starting", start_time=time.time()
            )

    def remove_bar(self, task_id: str) -> None:
        """Remove a progress bar by task id.

        If the task_id does not exist this is a no-op.
        """
        with self._lock:
            if task_id in self._bars:
                del self._bars[task_id]
                self._order.remove(task_id)

    def update(
        self,
        task_id: str,
        current: float = None,
        total: float = None,
        status: str = None,
        speed: str = None,
        color: str = None,
        stage: str = None,
    ) -> None:
        """Update a specific bar's state."""
        with self._lock:
            if task_id in self._bars:
                state = self._bars[task_id]
                if current is not None:
                    state.current = current
                if total is not None:
                    state.total = total
                if status is not None:
                    state.status = status
                if speed is not None:
                    state.speed = speed
                if color is not None:
                    state.color = color
                if stage is not None:
                    state.stage = stage

    def clear(self) -> None:
        """Remove all bars and reset render state."""
        with self._lock:
            self._bars.clear()
            self._order.clear()
            self._last_lines = 0
            self._footer = ""

    def set_footer(self, text: str) -> None:
        """Update the footer text."""
        with self._lock:
            self._footer = text

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self) -> None:
        """Start the rendering thread."""
        if self._started and not self._stopped:
            return  # Already running
        self._stop_event.clear()
        self._started = True
        self._stopped = False
        self._save_terminal_state()
        self._monitor_thread = threading.Thread(target=self._render_loop, daemon=True, name="tui-render")
        self._monitor_thread.start()

    def stop(self) -> None:
        """Stop rendering and cleanup.

        Safe to call multiple times. Joins the render thread with a timeout
        to avoid blocking indefinitely.
        """
        if self._stopped:
            return
        self._stopped = True
        self._stop_event.set()
        if self._monitor_thread is not None:
            self._monitor_thread.join(timeout=2.0)
            self._monitor_thread = None
        self._restore_terminal_state()

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def __enter__(self) -> "TerminalInterface":
        self.start()
        return self

    def __exit__(self, exc_type: type | None, exc_val: BaseException | None, exc_tb: object) -> None:
        self.stop()

    def __del__(self) -> None:
        """Safety net: restore terminal if the instance is garbage-collected without stop()."""
        if self._started and not self._stopped:
            try:
                self.stop()
            except Exception:
                # Best-effort in __del__; swallow any errors
                pass

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def _render_loop(self) -> None:
        """Loop to render the UI at 10Hz. Uses Event.wait for clean shutdown."""
        while not self._stop_event.is_set():
            self._render()
            # Event.wait returns True when the event is set (stop requested).
            # Timeout controls the refresh rate.
            self._stop_event.wait(timeout=0.1)

    def _render(self) -> None:
        """Render all bars to stdout."""
        if not self._is_tty():
            return

        with self._lock:
            # Move cursor up to overwrite previous items
            if self._last_lines > 0:
                sys.stdout.write(f"{CSI}{self._last_lines}A")

            lines: list[str] = []
            width = self._get_terminal_width()

            # Header
            header = f"{BOLD}{CYAN}MetaInformAnt SRA Downloader {RESET}"
            lines.append(header.center(width + len(BOLD) + len(CYAN) + len(RESET)))  # approximate centering fix
            lines.append("-" * width)

            now = time.time()

            for task_id in self._order:
                state = self._bars[task_id]

                # Calculate percentage
                pct = 0.0
                if state.total > 0:
                    pct = (state.current / state.total) * 100

                # Elapsed time
                elapsed_secs = now - state.start_time if state.start_time > 0 else 0.0
                elapsed_str = _format_elapsed(elapsed_secs)

                # Determine bar width
                # Layout: [Label] [Bar] [Pct] [Elapsed] [Speed] [Status]
                label_w = 20
                pct_w = 8
                elapsed_w = 7
                speed_w = 12
                status_w = 25  # Increased for TID info
                padding = 8  # Spaces between columns

                remaining_w = width - (label_w + pct_w + elapsed_w + speed_w + status_w + padding)
                if remaining_w < 10:
                    remaining_w = 10

                # Create Bar
                eff_pct = min(pct, 100.0)  # Cap for visual bar
                filled_len = int(remaining_w * eff_pct / 100)
                bar_str = "\u2588" * filled_len + "\u2591" * (remaining_w - filled_len)

                # Format strings
                label_str = state.label[:label_w].ljust(label_w)
                display_pct = min(pct, 100.0)
                pct_str = f"{display_pct:6.1f}%"
                elapsed_fmt = elapsed_str.rjust(elapsed_w)
                speed_str = state.speed.rjust(speed_w)
                status_str = state.status[:status_w].ljust(status_w)

                line = (
                    f"{state.color}{label_str}{RESET} "
                    f"{state.color}{bar_str}{RESET} "
                    f"{BOLD}{pct_str}{RESET} "
                    f"{MAGENTA}{elapsed_fmt}{RESET} "
                    f"{YELLOW}{speed_str}{RESET} "
                    f"{WHITE}{status_str}{RESET}"
                )
                lines.append(line)

            # Footer
            if self._footer:
                lines.append("-" * width)
                lines.append(f"{self._footer}".center(width))

            # Output everything at once
            output = "\n".join(lines) + "\n"
            try:
                sys.stdout.write(output)
                sys.stdout.flush()
            except (OSError, ValueError):
                # stdout closed or not writable -- silently skip
                pass

            self._last_lines = len(lines)


# ---------------------------------------------------------------------------
# Module-level atexit handler: ensure all live instances are cleaned up
# ---------------------------------------------------------------------------


def _atexit_restore_all() -> None:
    """Restore terminal state for any TerminalInterface instances still alive at exit."""
    for instance in list(TerminalInterface._live_instances):
        try:
            instance.stop()
        except Exception:
            # Last-resort: try to show cursor directly
            try:
                sys.stdout.write(SHOW_CURSOR)
                sys.stdout.flush()
            except Exception:
                pass


atexit.register(_atexit_restore_all)
