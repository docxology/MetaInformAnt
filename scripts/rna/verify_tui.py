"""Verification script for Terminal UI.

Simulates parallel downloads to demonstrate the TUI.
"""

import time
from pathlib import Path

from metainformant.core.ui.tui import BLUE, GREEN, YELLOW, TerminalInterface


def test_tui():
    ui = TerminalInterface()

    # Add fake tasks
    ui.add_bar("task1", "SRR12345678", 100)
    ui.add_bar("task2", "SRR87654321", 100)
    ui.add_bar("task3", "SRR00000000", 100)

    ui.start()

    try:
        for i in range(101):
            ui.update("task1", current=float(i), status="Downloading", speed=f"{i*0.5:.1f} MB/s", color=BLUE)
            if i > 20:
                ui.update(
                    "task2", current=float(i - 20), status="Downloading", speed=f"{(i-20)*0.8:.1f} MB/s", color=YELLOW
                )
            if i > 50:
                ui.update(
                    "task3", current=float(i - 50), status="Downloading", speed=f"{(i-50)*1.2:.1f} MB/s", color=GREEN
                )

            time.sleep(0.05)

        ui.update("task1", status="Done", color=GREEN)
        ui.update("task2", status="Done", color=GREEN)
        ui.update("task3", status="Done", color=GREEN)
        time.sleep(1)

    finally:
        ui.stop()


if __name__ == "__main__":
    test_tui()
