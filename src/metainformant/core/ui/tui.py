"""Core Terminal UI components for MetaInformAnt.

This module provides a robust, dependency-free (standard library only)
Terminal User Interface for tracking parallel processes.
It uses ANSI escape codes for "oldschool" elegance.
"""

import sys
import threading
import time
import shutil
from dataclasses import dataclass
from typing import Dict, List, Optional

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
    stage: str = ""  # For multi-stage workflows (e.g., "Download → Quant")


class TerminalInterface:
    """Manages terminal output for parallel tasks."""

    def __init__(self):
        self._lock = threading.Lock()
        self._bars: Dict[str, ProgressState] = {}
        self._order: List[str] = []
        self._running = False
        self._last_lines = 0
        self._footer = ""
        self._monitor_thread: Optional[threading.Thread] = None

    def add_bar(self, task_id: str, label: str, total: float, unit: str = "B"):
        """Register a new progress bar."""
        with self._lock:
            if task_id not in self._bars:
                self._order.append(task_id)
            self._bars[task_id] = ProgressState(
                label=label, current=0, total=total, unit=unit, status="Starting", start_time=time.time()
            )

    def update(
        self,
        task_id: str,
        current: float = None,
        total: float = None,
        status: str = None,
        speed: str = None,
        color: str = None,
        stage: str = None,
    ):
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

    def start(self):
        """Start the rendering thread."""
        self._running = True
        sys.stdout.write(HIDE_CURSOR)
        self._monitor_thread = threading.Thread(target=self._render_loop, daemon=True)
        self._monitor_thread.start()

    def stop(self):
        """Stop rendering and cleanup."""
        self._running = False
        if self._monitor_thread:
            self._monitor_thread.join()
        sys.stdout.write(SHOW_CURSOR)
        sys.stdout.write("\n")

    def _render_loop(self):
        """Loop to render the UI at 10Hz."""
        while self._running:
            self._render()
            time.sleep(0.1)

    def _render(self):
        """Render all bars to stdout."""
        with self._lock:
            # Move cursor up to overwrite previous items
            if self._last_lines > 0:
                sys.stdout.write(f"{CSI}{self._last_lines}A")

            lines = []
            width = shutil.get_terminal_size().columns

            # Header
            header = f"{BOLD}{CYAN}MetaInformAnt SRA Downloader {RESET}"
            lines.append(header.center(width + len(BOLD) + len(CYAN) + len(RESET)))  # approximate centering fix
            lines.append("-" * width)

            for task_id in self._order:
                state = self._bars[task_id]

                # Calculate percentage
                pct = 0.0
                if state.total > 0:
                    pct = (state.current / state.total) * 100

                # Determine bar width
                # Layout: [Label] [Bar] [Pct] [Speed] [Status]
                # Fixed widths: Label(20), Pct(8), Speed(12), Status(15)
                # Remainder for bar

                label_w = 20
                pct_w = 8
                speed_w = 12
                status_w = 25  # Increased for TID info
                padding = 6  # Spaces between columns

                remaining_w = width - (label_w + pct_w + speed_w + status_w + padding)
                if remaining_w < 10:
                    remaining_w = 10

                # Create Bar
                eff_pct = min(pct, 100.0)  # Cap for visual bar
                filled_len = int(remaining_w * eff_pct / 100)
                bar_str = "█" * filled_len + "░" * (remaining_w - filled_len)

                # Format strings
                label_str = state.label[:label_w].ljust(label_w)
                display_pct = min(pct, 100.0)
                pct_str = f"{display_pct:6.1f}%"
                speed_str = state.speed.rjust(speed_w)
                status_str = state.status[:status_w].ljust(status_w)

                line = (
                    f"{state.color}{label_str}{RESET} "
                    f"{state.color}{bar_str}{RESET} "
                    f"{BOLD}{pct_str}{RESET} "
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
            sys.stdout.write(output)
            sys.stdout.flush()

            self._last_lines = len(lines)

    def set_footer(self, text: str):
        """Update the footer text."""
        with self._lock:
            self._footer = text
