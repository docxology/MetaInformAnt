"""Progress tracking for RNA-seq workflows.

This module provides utilities for tracking progress and status of RNA-seq analysis workflows.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)


class ProgressTracker:
    """Tracks progress of RNA-seq workflow steps.

    Provides methods for monitoring workflow execution, estimating completion times,
    and reporting status updates.
    """

    def __init__(self, total_steps: Optional[int] = None, description: str = "RNA-seq workflow"):
        """Initialize progress tracker.

        Args:
            total_steps: Total number of steps (if known)
            description: Description of the workflow being tracked
        """
        self.total_steps = total_steps
        self.description = description
        self.completed_steps = 0
        self.start_time = time.time()
        self.step_times = []
        self.current_step_start = None
        self.step_descriptions = []

    def start_step(self, step_name: str) -> None:
        """Mark the start of a workflow step.

        Args:
            step_name: Name of the step being started
        """
        if self.current_step_start is not None:
            # End previous step
            step_time = time.time() - self.current_step_start
            self.step_times.append(step_time)

        self.current_step_start = time.time()
        self.step_descriptions.append(step_name)

        logger.info(f"Started step: {step_name}")

    def complete_step(self, step_name: Optional[str] = None) -> None:
        """Mark the completion of a workflow step.

        Args:
            step_name: Name of the completed step (optional, for logging)
        """
        if self.current_step_start is not None:
            step_time = time.time() - self.current_step_start
            self.step_times.append(step_time)
            self.completed_steps += 1

            step_desc = step_name or f"Step {self.completed_steps}"
            logger.info(".2f")

            self.current_step_start = None

    def get_progress(self) -> Dict[str, Any]:
        """Get current progress information.

        Returns:
            Dictionary with progress statistics
        """
        elapsed_time = time.time() - self.start_time
        progress_info = {
            'description': self.description,
            'completed_steps': self.completed_steps,
            'total_steps': self.total_steps,
            'elapsed_time': elapsed_time,
            'average_step_time': sum(self.step_times) / len(self.step_times) if self.step_times else 0,
        }

        if self.total_steps:
            progress_info['percent_complete'] = (self.completed_steps / self.total_steps) * 100

            # Estimate remaining time
            if self.step_times:
                avg_step_time = sum(self.step_times) / len(self.step_times)
                remaining_steps = self.total_steps - self.completed_steps
                progress_info['estimated_remaining_time'] = remaining_steps * avg_step_time

                if progress_info.get('percent_complete', 0) > 0:
                    total_estimated_time = elapsed_time / (progress_info['percent_complete'] / 100)
                    progress_info['estimated_total_time'] = total_estimated_time

        return progress_info

    def log_progress(self) -> None:
        """Log current progress information."""
        progress = self.get_progress()

        if self.total_steps:
            percent = progress.get('percent_complete', 0)
            logger.info(f"Progress: {percent:.1f}% complete "
                       f"({progress['completed_steps']}/{self.total_steps} steps, "
                       f"elapsed: {progress['elapsed_time']:.1f}s)")
        else:
            logger.info(f"Progress: {progress['completed_steps']} steps completed "
                       f"(elapsed: {progress['elapsed_time']:.1f}s)")

        # Log step timing information
        if self.step_times:
            total_step_time = sum(self.step_times)
            logger.debug(f"Total step time: {total_step_time:.2f}s, "
                        f"Average step time: {progress['average_step_time']:.2f}s")

    def reset(self) -> None:
        """Reset the progress tracker."""
        self.completed_steps = 0
        self.start_time = time.time()
        self.step_times = []
        self.current_step_start = None
        self.step_descriptions = []

        logger.info("Progress tracker reset")

    def __str__(self) -> str:
        progress = self.get_progress()
        if self.total_steps:
            return ".1f"
        else:
            return f"ProgressTracker({progress['completed_steps']} steps completed, " \
                   f"{progress['elapsed_time']:.1f}s elapsed)"


def get_tracker(total_steps: Optional[int] = None, description: str = "RNA-seq workflow") -> ProgressTracker:
    """Factory function to create a ProgressTracker instance.

    Args:
        total_steps: Total number of steps (if known)
        description: Description of the workflow being tracked

    Returns:
        ProgressTracker instance
    """
    return ProgressTracker(total_steps=total_steps, description=description)
