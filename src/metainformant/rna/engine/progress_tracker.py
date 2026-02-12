"""Progress tracking for RNA-seq workflows.

This module provides utilities for tracking progress and status of RNA-seq analysis workflows.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional, Set

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


class ProgressTracker:
    """Tracks progress of RNA-seq workflow steps.

    Provides methods for monitoring workflow execution, estimating completion times,
    and reporting status updates. Supports both high-level step tracking and
    granular per-sample tracking for SRA extraction and quantification.
    """

    def __init__(
        self,
        total_steps: Optional[int] = None,
        description: str = "RNA-seq workflow",
        state_file: Optional[Path] = None,
        dashboard_file: Optional[Path] = None,
    ):
        """Initialize progress tracker.

        Args:
            total_steps: Total number of steps (if known)
            description: Description of the workflow being tracked
            state_file: Path to JSON file for persisting state
            dashboard_file: Path to text file for status dashboard
        """
        self.total_steps = total_steps
        self.description = description
        self.state_file = Path(state_file) if state_file else None
        self.dashboard_file = Path(dashboard_file) if dashboard_file else None

        self.completed_steps = 0
        self.start_time = time.time()
        self.step_times = []
        self.current_step_start = None
        self.step_descriptions = []

        # State for granular tracking (per-species, per-sample)
        self.state = {}
        self._load_state()

    def _load_state(self):
        """Load state from file if it exists."""
        if self.state_file and self.state_file.exists():
            try:
                content = self.state_file.read_text()
                loaded = json.loads(content)
                # Convert list back to sets for need_download etc
                for species, data in loaded.items():
                    for key in ["need_download", "ongoing_download", "needs_quant", "needs_delete", "completed"]:
                        if key in data:
                            data[key] = set(data[key])
                self.state = loaded
            except Exception as e:
                logger.warning(f"Failed to load progress state: {e}")

    def _save_state(self):
        """Save state to file for persistence."""
        if not self.state_file:
            return

        try:
            self.state_file.parent.mkdir(parents=True, exist_ok=True)
            # Convert sets to lists for JSON serialization
            serialized = {}
            for species, data in self.state.items():
                species_data = data.copy()
                for key in ["need_download", "ongoing_download", "needs_quant", "needs_delete", "completed"]:
                    if key in species_data:
                        species_data[key] = list(species_data[key])
                serialized[species] = species_data

            self.state_file.write_text(json.dumps(serialized, indent=2))
        except Exception as e:
            logger.warning(f"Failed to save progress state: {e}")

    def initialize_species(self, species: str, total_samples: int, sample_ids: list[str]):
        """Initialize tracking for a specific species."""
        self.state[species] = {
            "total": total_samples,
            "need_download": set(sample_ids),
            "ongoing_download": set(),
            "needs_quant": set(),
            "needs_delete": set(),
            "completed": set(),
        }
        self._save_state()

    def on_download_start(self, species: str, sample_id: str):
        if species in self.state:
            self.state[species]["need_download"].discard(sample_id)
            self.state[species]["ongoing_download"].add(sample_id)
            self._save_state()

    def on_download_complete(self, species: str, sample_id: str):
        if species in self.state:
            self.state[species]["ongoing_download"].discard(sample_id)
            self.state[species]["needs_quant"].add(sample_id)
            self._save_state()

    def on_quant_complete(self, species: str, sample_id: str):
        if species in self.state:
            self.state[species]["needs_quant"].discard(sample_id)
            self.state[species]["needs_delete"].add(sample_id)
            self._save_state()

    def on_delete_complete(self, species: str, sample_id: str):
        if species in self.state:
            self.state[species]["needs_delete"].discard(sample_id)
            self.state[species]["completed"].add(sample_id)
            self._save_state()

    def get_species_state(self, species: str) -> Dict[str, Any]:
        if species not in self.state:
            return {}
        data = self.state[species]
        return {
            "total": data["total"],
            "completed": len(data["completed"]),
            "need_download": len(data["need_download"]),
            "ongoing_download": len(data["ongoing_download"]),
            "needs_quant": len(data["needs_quant"]),
            "needs_delete": len(data["needs_delete"]),
        }

    def start_step(self, step_name: str) -> None:
        """Mark the start of a workflow step."""
        if self.current_step_start is not None:
            # End previous step
            step_time = time.time() - self.current_step_start
            self.step_times.append(step_time)

        self.current_step_start = time.time()
        self.step_descriptions.append(step_name)

        logger.info(f"Started step: {step_name}")

    def complete_step(self, step_name: Optional[str] = None) -> None:
        """Mark the completion of a workflow step."""
        if self.current_step_start is not None:
            step_time = time.time() - self.current_step_start
            self.step_times.append(step_time)
            self.completed_steps += 1

            step_desc = step_name or f"Step {self.completed_steps}"
            logger.info(f"Completed {step_desc} in {step_time:.2f}s")

            self.current_step_start = None

    def get_progress(self) -> Dict[str, Any]:
        """Get current progress information."""
        elapsed_time = time.time() - self.start_time
        progress_info = {
            "description": self.description,
            "completed_steps": self.completed_steps,
            "total_steps": self.total_steps,
            "elapsed_time": elapsed_time,
            "average_step_time": sum(self.step_times) / len(self.step_times) if self.step_times else 0,
        }

        if self.total_steps:
            progress_info["percent_complete"] = (self.completed_steps / self.total_steps) * 100

            # Estimate remaining time
            if self.step_times:
                avg_step_time = sum(self.step_times) / len(self.step_times)
                remaining_steps = self.total_steps - self.completed_steps
                progress_info["estimated_remaining_time"] = remaining_steps * avg_step_time

                if progress_info.get("percent_complete", 0) > 0:
                    total_estimated_time = elapsed_time / (progress_info["percent_complete"] / 100)
                    progress_info["estimated_total_time"] = total_estimated_time

        return progress_info

    def log_progress(self) -> None:
        """Log current progress information."""
        progress = self.get_progress()

        if self.total_steps:
            percent = progress.get("percent_complete", 0)
            logger.info(
                f"Progress: {percent:.1f}% complete "
                f"({progress['completed_steps']}/{self.total_steps} steps, "
                f"elapsed: {progress['elapsed_time']:.1f}s)"
            )
        else:
            logger.info(
                f"Progress: {progress['completed_steps']} steps completed "
                f"(elapsed: {progress['elapsed_time']:.1f}s)"
            )

        # Log step timing information
        if self.step_times:
            total_step_time = sum(self.step_times)
            logger.debug(
                f"Total step time: {total_step_time:.2f}s, " f"Average step time: {progress['average_step_time']:.2f}s"
            )

    def reset(self) -> None:
        """Reset the progress tracker."""
        self.completed_steps = 0
        self.start_time = time.time()
        self.step_times = []
        self.current_step_start = None
        self.step_descriptions = []
        self.state = {}
        if self.state_file and self.state_file.exists():
            self.state_file.unlink()

        logger.info("Progress tracker reset")

    def __str__(self) -> str:
        progress = self.get_progress()
        if self.total_steps:
            return f"Progress: {progress.get('percent_complete', 0):.1f}% complete"
        else:
            return (
                f"ProgressTracker({progress['completed_steps']} steps completed, "
                f"{progress['elapsed_time']:.1f}s elapsed)"
            )


_tracker_instance = None


def get_tracker(
    total_steps: Optional[int] = None,
    description: str = "RNA-seq workflow",
    state_file: Optional[Path] = None,
    dashboard_file: Optional[Path] = None,
) -> ProgressTracker:
    """Factory function to create or get a ProgressTracker instance.

    Uses a singleton pattern for the global tracker if no paths are provided.
    """
    global _tracker_instance
    if state_file or dashboard_file:
        return ProgressTracker(
            total_steps=total_steps, description=description, state_file=state_file, dashboard_file=dashboard_file
        )

    if _tracker_instance is None:
        _tracker_instance = ProgressTracker(total_steps=total_steps, description=description)
    return _tracker_instance
