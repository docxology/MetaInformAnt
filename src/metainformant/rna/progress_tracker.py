"""Real-time progress tracking and visualization for RNA-seq workflows.

This module provides event-driven progress tracking that reports sample completion
events (download, quantification, deletion) and generates a continuously updated
visualization showing per-species sample counts in each category.
"""

from __future__ import annotations

import json
import threading
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    from ...core.io import read_delimited
    from ...core.logging import get_logger
except ImportError:
    from metainformant.core.io import read_delimited
    from metainformant.core.logging import get_logger

logger = get_logger("progress_tracker")


class ProgressTracker:
    """Tracks RNA-seq workflow progress per species with categorized sample states.
    
    Maintains state for each species:
    - need_download: Samples not yet downloaded
    - ongoing_download: Samples currently downloading
    - failed_download: Samples that failed to download
    - needs_quant: Samples downloaded but not quantified
    - needs_delete: Samples quantified but FASTQs still present
    - completed: Samples quantified and deleted
    """
    
    def __init__(self, state_file: Path | None = None, dashboard_file: Path | None = None):
        """Initialize progress tracker.
        
        Args:
            state_file: Path to JSON file for state persistence (default: output/amalgkit/progress_state.json)
            dashboard_file: Path to dashboard visualization file (default: output/amalgkit/progress_dashboard.txt)
        """
        repo_root = Path(__file__).parent.parent.parent.parent.resolve()
        output_dir = repo_root / "output"
        amalgkit_dir = output_dir / "amalgkit"
        
        self.state_file = state_file or (amalgkit_dir / "progress_state.json")
        self.dashboard_file = dashboard_file or (amalgkit_dir / "progress_dashboard.txt")
        self.log_file = amalgkit_dir / "progress_tracker.log"
        
        # Per-species state: species_name -> dict with category sets
        self.state: dict[str, dict[str, Any]] = defaultdict(lambda: {
            "need_download": set(),
            "ongoing_download": set(),
            "failed_download": set(),
            "needs_quant": set(),
            "needs_delete": set(),
            "completed": set(),
            "total": 0,
            "last_update": None,
        })
        
        # Thread lock for concurrent access
        self.lock = threading.Lock()
        
        # Load existing state if available
        self._load_state()
    
    def _load_state(self) -> None:
        """Load state from disk if available."""
        if self.state_file.exists():
            try:
                with open(self.state_file) as f:
                    data = json.load(f)
                
                # Convert list sets back to sets and handle datetime strings
                for species, species_data in data.items():
                    for key in ["need_download", "ongoing_download", "failed_download", 
                               "needs_quant", "needs_delete", "completed"]:
                        if key in species_data:
                            species_data[key] = set(species_data[key])
                    # Handle datetime string
                    if "last_update" in species_data and species_data["last_update"]:
                        if isinstance(species_data["last_update"], str):
                            try:
                                from datetime import datetime
                                species_data["last_update"] = datetime.fromisoformat(species_data["last_update"])
                            except Exception:
                                species_data["last_update"] = None
                    self.state[species] = species_data
                
                logger.info(f"Loaded progress state from {self.state_file}")
            except Exception as e:
                logger.warning(f"Failed to load progress state: {e}")
    
    def _save_state(self) -> None:
        """Save state to disk."""
        try:
            # Convert sets to lists for JSON serialization
            data = {}
            for species, species_data in self.state.items():
                data[species] = {
                    "need_download": list(species_data["need_download"]),
                    "ongoing_download": list(species_data["ongoing_download"]),
                    "failed_download": list(species_data["failed_download"]),
                    "needs_quant": list(species_data["needs_quant"]),
                    "needs_delete": list(species_data["needs_delete"]),
                    "completed": list(species_data["completed"]),
                    "total": species_data["total"],
                    "last_update": species_data["last_update"].isoformat() if species_data["last_update"] else None,
                }
            
            # Ensure directory exists
            self.state_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(self.state_file, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save progress state: {e}")
    
    def initialize_species(self, species: str, total_samples: int, sample_ids: list[str]) -> None:
        """Initialize tracking for a species.
        
        Args:
            species: Species identifier (e.g., "cfloridanus")
            total_samples: Total number of samples for this species
            sample_ids: List of all sample IDs (SRA run IDs)
        """
        with self.lock:
            self.state[species]["total"] = total_samples
            self.state[species]["need_download"] = set(sample_ids)
            self.state[species]["last_update"] = datetime.now()
            self._save_state()
            logger.debug(f"Initialized {species}: {total_samples} samples")
    
    def on_download_start(self, species: str, sample_id: str) -> None:
        """Report that a download has started.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID (e.g., "SRR1234567")
        """
        with self.lock:
            if sample_id in self.state[species]["need_download"]:
                self.state[species]["need_download"].discard(sample_id)
            self.state[species]["ongoing_download"].add(sample_id)
            self.state[species]["last_update"] = datetime.now()
            self._log_event(species, sample_id, "download_start")
    
    def on_download_complete(self, species: str, sample_id: str) -> None:
        """Report that a download has completed.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID
        """
        with self.lock:
            self.state[species]["ongoing_download"].discard(sample_id)
            self.state[species]["failed_download"].discard(sample_id)
            self.state[species]["needs_quant"].add(sample_id)
            self.state[species]["last_update"] = datetime.now()
            self._log_event(species, sample_id, "download_complete")
            self._save_state()
            self.update_dashboard()
    
    def on_download_failed(self, species: str, sample_id: str) -> None:
        """Report that a download has failed.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID
        """
        with self.lock:
            self.state[species]["ongoing_download"].discard(sample_id)
            self.state[species]["failed_download"].add(sample_id)
            self.state[species]["last_update"] = datetime.now()
            self._log_event(species, sample_id, "download_failed")
            self._save_state()
            self.update_dashboard()
    
    def on_quant_complete(self, species: str, sample_id: str) -> None:
        """Report that quantification has completed.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID
        """
        with self.lock:
            self.state[species]["needs_quant"].discard(sample_id)
            self.state[species]["needs_delete"].add(sample_id)
            self.state[species]["last_update"] = datetime.now()
            self._log_event(species, sample_id, "quant_complete")
            self._save_state()
            self.update_dashboard()
    
    def on_delete_complete(self, species: str, sample_id: str) -> None:
        """Report that FASTQ deletion has completed.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID
        """
        with self.lock:
            self.state[species]["needs_delete"].discard(sample_id)
            self.state[species]["completed"].add(sample_id)
            self.state[species]["last_update"] = datetime.now()
            self._log_event(species, sample_id, "delete_complete")
            self._save_state()
            self.update_dashboard()
    
    def get_species_state(self, species: str) -> dict[str, Any]:
        """Get current state for a species.
        
        Args:
            species: Species identifier
            
        Returns:
            Dictionary with category counts and sets
        """
        with self.lock:
            state = self.state[species]
            return {
                "need_download": len(state["need_download"]),
                "ongoing_download": len(state["ongoing_download"]),
                "failed_download": len(state["failed_download"]),
                "needs_quant": len(state["needs_quant"]),
                "needs_delete": len(state["needs_delete"]),
                "completed": len(state["completed"]),
                "total": state["total"],
                "last_update": state["last_update"],
            }
    
    def get_all_species_state(self) -> dict[str, dict[str, Any]]:
        """Get current state for all species.
        
        Returns:
            Dictionary mapping species to their state
        """
        with self.lock:
            return {species: self.get_species_state(species) for species in self.state.keys()}
    
    def _log_event(self, species: str, sample_id: str, event_type: str) -> None:
        """Log an event to the progress tracker log file with enhanced context.
        
        Args:
            species: Species identifier
            sample_id: SRA run ID
            event_type: Type of event (e.g., "download_complete")
        """
        try:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Get current state for context
            state = self.get_species_state(species)
            completed = state.get("completed", 0)
            total = state.get("total", 0)
            pct = (completed / total * 100) if total > 0 else 0.0
            
            log_entry = (f"{timestamp} | INFO | progress_tracker | {species:20s} | {sample_id:12s} | "
                        f"{event_type:20s} | Progress: {completed:4d}/{total:4d} ({pct:5.1f}%)\n")
            
            # Ensure directory exists
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(self.log_file, 'a') as f:
                f.write(log_entry)
            
            # Also log to main logger at debug level
            logger.debug(f"{species}/{sample_id}: {event_type} ({completed}/{total} = {pct:.1f}%)")
        except Exception as e:
            logger.debug(f"Failed to log event: {e}")
    
    def generate_dashboard(self) -> str:
        """Generate ASCII dashboard visualization with enhanced statistics.
        
        Returns:
            Dashboard text as string
        """
        lines = []
        lines.append("=" * 100)
        lines.append("RNA-SEQ WORKFLOW PROGRESS DASHBOARD")
        lines.append(f"Updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("=" * 100)
        lines.append("")
        
        # Calculate overall statistics
        total_samples = 0
        total_completed = 0
        total_ongoing = 0
        total_failed = 0
        total_needs_quant = 0
        total_needs_delete = 0
        
        with self.lock:
            for species_data in self.state.values():
                total_samples += species_data["total"]
                total_completed += len(species_data["completed"])
                total_ongoing += len(species_data["ongoing_download"])
                total_failed += len(species_data["failed_download"])
                total_needs_quant += len(species_data["needs_quant"])
                total_needs_delete += len(species_data["needs_delete"])
        
        overall_pct = (total_completed / total_samples * 100) if total_samples > 0 else 0.0
        in_progress = total_ongoing + total_needs_quant + total_needs_delete
        
        lines.append("OVERALL STATISTICS")
        lines.append("-" * 100)
        lines.append(f"Total Samples: {total_samples}")
        lines.append(f"Completed: {total_completed} ({overall_pct:.1f}%)")
        lines.append(f"In Progress: {in_progress} (downloading: {total_ongoing}, needs quant: {total_needs_quant}, needs delete: {total_needs_delete})")
        lines.append(f"Failed Downloads: {total_failed}")
        lines.append(f"Remaining: {total_samples - total_completed - in_progress}")
        lines.append("")
        lines.append("=" * 100)
        lines.append("")
        
        # Table header
        lines.append(f"{'Species':<35} {'Total':>6} {'Need DL':>9} {'Downloading':>12} {'Failed':>7} "
                    f"{'Needs Quant':>12} {'Needs Delete':>13} {'Completed':>10} {'%':>6}")
        lines.append("-" * 100)
        
        # Per-species rows
        species_states = []
        with self.lock:
            for species in sorted(self.state.keys()):
                state = self.get_species_state(species)
                species_states.append((species, state))
        
        for species, state in species_states:
            pct = (state["completed"] / state["total"] * 100) if state["total"] > 0 else 0.0
            species_display = species.replace("_", " ").title()
            
            # Color coding symbols
            status_symbol = "✅" if pct == 100 else "🔄" if state["ongoing_download"] > 0 or state["needs_quant"] > 0 else "⏳"
            
            lines.append(f"{status_symbol} {species_display:<33} {state['total']:>6} {state['need_download']:>9} "
                        f"{state['ongoing_download']:>12} {state['failed_download']:>7} "
                        f"{state['needs_quant']:>12} {state['needs_delete']:>13} "
                        f"{state['completed']:>10} {pct:>5.1f}%")
        
        lines.append("")
        lines.append("=" * 100)
        
        # Progress bars for all species (sorted by completion)
        if species_states:
            lines.append("")
            lines.append("PROGRESS BARS (sorted by completion %):")
            lines.append("")
            
            # Sort by completion percentage
            sorted_species = sorted(species_states, key=lambda x: x[1]["completed"] / max(1, x[1]["total"]), reverse=True)
            
            for species, state in sorted_species:
                pct = (state["completed"] / state["total"] * 100) if state["total"] > 0 else 0.0
                bar_length = 50
                filled = int(bar_length * state["completed"] / max(1, state["total"]))
                bar = "█" * filled + "░" * (bar_length - filled)
                species_display = species.replace("_", " ").title()
                lines.append(f"{species_display:<35} [{bar}] {pct:>5.1f}% ({state['completed']}/{state['total']})")
        
        lines.append("")
        lines.append("=" * 100)
        lines.append("")
        lines.append("LEGEND:")
        lines.append("  ✅ = Complete (100%)")
        lines.append("  🔄 = In Progress (downloading, quantifying, or deleting)")
        lines.append("  ⏳ = Waiting (not started or paused)")
        lines.append("")
        lines.append("CATEGORIES:")
        lines.append("  Need DL = Samples not yet downloaded")
        lines.append("  Downloading = Currently downloading from ENA/SRA")
        lines.append("  Failed = Downloads that failed (will retry)")
        lines.append("  Needs Quant = Downloaded but not yet quantified")
        lines.append("  Needs Delete = Quantified but FASTQs still present")
        lines.append("  Completed = Quantified and FASTQs deleted")
        lines.append("")
        lines.append(f"Dashboard Location: {self.dashboard_file}")
        lines.append(f"State File: {self.state_file}")
        lines.append(f"Log File: {self.log_file}")
        lines.append("")
        
        return "\n".join(lines)
    
    def generate_png_visualization(self) -> Path | None:
        """Generate PNG visualization of workflow progress.
        
        Creates a multi-panel visualization showing:
        - Overall progress bar
        - Per-species progress bars
        - Category breakdown (completed, in progress, failed)
        
        Returns:
            Path to generated PNG file, or None if generation failed
        """
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.debug("matplotlib not available, skipping PNG visualization")
            return None
        
        try:
            # Get all species states
            species_states = self.get_all_species_states()
            if not species_states:
                return None
            
            # Prepare data
            species_names = sorted(species_states.keys())
            totals = [species_states[s]["total"] for s in species_names]
            completed = [species_states[s]["completed"] for s in species_names]
            percentages = [(c / t * 100) if t > 0 else 0.0 for c, t in zip(completed, totals)]
            
            # Calculate category totals
            total_samples = sum(totals)
            total_completed = sum(completed)
            total_ongoing = sum(species_states[s]["ongoing_download"] + 
                               species_states[s]["needs_quant"] + 
                               species_states[s]["needs_delete"] for s in species_names)
            total_failed = sum(species_states[s]["failed_download"] for s in species_names)
            total_remaining = total_samples - total_completed - total_ongoing
            
            # Create figure with subplots
            fig = plt.figure(figsize=(14, 10))
            gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
            
            # Panel 1: Overall progress bar
            ax1 = fig.add_subplot(gs[0, :])
            overall_pct = (total_completed / total_samples * 100) if total_samples > 0 else 0.0
            ax1.barh([0], [overall_pct], color='#2ecc71', height=0.5, label='Completed')
            ax1.barh([0], [(total_ongoing / total_samples * 100) if total_samples > 0 else 0.0], 
                    left=[overall_pct], color='#f39c12', height=0.5, label='In Progress')
            ax1.barh([0], [(total_failed / total_samples * 100) if total_samples > 0 else 0.0],
                    left=[overall_pct + (total_ongoing / total_samples * 100) if total_samples > 0 else 0.0],
                    color='#e74c3c', height=0.5, label='Failed')
            ax1.set_xlim(0, 100)
            ax1.set_ylim(-0.5, 0.5)
            ax1.set_xlabel('Progress (%)', fontsize=12, fontweight='bold')
            ax1.set_title(f'Overall Progress: {total_completed}/{total_samples} samples ({overall_pct:.1f}%)', 
                         fontsize=14, fontweight='bold')
            ax1.legend(loc='upper right')
            ax1.set_yticks([])
            ax1.grid(axis='x', alpha=0.3)
            
            # Panel 2: Per-species progress bars
            ax2 = fig.add_subplot(gs[1, :])
            y_pos = np.arange(len(species_names))
            colors = ['#2ecc71' if p == 100 else '#f39c12' if p > 0 else '#95a5a6' 
                     for p in percentages]
            bars = ax2.barh(y_pos, percentages, color=colors, alpha=0.8)
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels([s.replace('_', ' ').title() for s in species_names], fontsize=9)
            ax2.set_xlabel('Progress (%)', fontsize=11)
            ax2.set_title('Per-Species Progress', fontsize=12, fontweight='bold')
            ax2.set_xlim(0, 100)
            ax2.grid(axis='x', alpha=0.3)
            
            # Add percentage labels on bars
            for i, (bar, pct) in enumerate(zip(bars, percentages)):
                if pct > 5:  # Only label if bar is visible
                    ax2.text(pct / 2, i, f'{pct:.1f}%', 
                            ha='center', va='center', fontweight='bold', fontsize=8)
            
            # Panel 3: Category breakdown (pie chart)
            ax3 = fig.add_subplot(gs[2, 0])
            categories = ['Completed', 'In Progress', 'Failed', 'Remaining']
            sizes = [total_completed, total_ongoing, total_failed, total_remaining]
            colors_pie = ['#2ecc71', '#f39c12', '#e74c3c', '#95a5a6']
            # Only show non-zero categories
            non_zero = [(c, s, col) for c, s, col in zip(categories, sizes, colors_pie) if s > 0]
            if non_zero:
                cats, szs, cols = zip(*non_zero)
                ax3.pie(szs, labels=cats, colors=cols, autopct='%1.1f%%', startangle=90)
                ax3.set_title('Sample Status Distribution', fontsize=11, fontweight='bold')
            
            # Panel 4: Statistics table
            ax4 = fig.add_subplot(gs[2, 1])
            ax4.axis('off')
            stats_text = [
                f"Total Samples: {total_samples:,}",
                f"Completed: {total_completed:,} ({total_completed/total_samples*100:.1f}%)" if total_samples > 0 else "Completed: 0",
                f"In Progress: {total_ongoing:,}",
                f"Failed Downloads: {total_failed:,}",
                f"Remaining: {total_remaining:,}",
                "",
                f"Last Update: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
            ]
            ax4.text(0.1, 0.9, '\n'.join(stats_text), transform=ax4.transAxes,
                    fontsize=10, verticalalignment='top', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Overall title
            fig.suptitle('RNA-Seq Workflow Progress Summary', fontsize=16, fontweight='bold', y=0.98)
            
            # Save PNG
            png_file = self.dashboard_file.with_suffix('.png')
            png_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(png_file, dpi=150, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
            logger.debug(f"PNG visualization generated: {png_file}")
            return png_file
            
        except Exception as e:
            logger.warning(f"Failed to generate PNG visualization: {e}")
            return None
    
    def update_dashboard(self) -> None:
        """Update the dashboard file and PNG visualization with current state."""
        import time
        try:
            update_start = time.time()
            dashboard_text = self.generate_dashboard()
            
            # Ensure directory exists
            self.dashboard_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Update text dashboard
            with open(self.dashboard_file, 'w') as f:
                f.write(dashboard_text)
            
            # Generate PNG visualization
            png_file = self.generate_png_visualization()
            if png_file:
                logger.debug(f"PNG visualization updated: {png_file}")
            
            update_time = time.time() - update_start
            if update_time > 1.0:
                logger.debug(f"Dashboard update took {update_time:.1f}s (slow)")
        except Exception as e:
            logger.warning(f"Failed to update dashboard: {e}")


# Global tracker instance (singleton pattern)
_tracker_instance: ProgressTracker | None = None
_tracker_lock = threading.Lock()


def get_tracker(state_file: Path | None = None, dashboard_file: Path | None = None) -> ProgressTracker:
    """Get or create the global progress tracker instance.
    
    Args:
        state_file: Optional path to state file
        dashboard_file: Optional path to dashboard file
        
    Returns:
        ProgressTracker instance
    """
    global _tracker_instance
    
    with _tracker_lock:
        if _tracker_instance is None:
            _tracker_instance = ProgressTracker(state_file=state_file, dashboard_file=dashboard_file)
        return _tracker_instance

