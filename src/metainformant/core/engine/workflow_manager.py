"""Workflow Manager for orchestrating multi-stage RNA-seq workflows with TUI.

This module provides a WorkflowManager that coordinates:
- Download phase (using DownloadManager)
- getfastq phase (FASTQ extraction via amalgkit)
- quant phase (quantification via amalgkit)

Each sample tracks its progress through all stages in the terminal UI.
"""

import threading
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

from metainformant.core.ui.tui import TerminalInterface, GREEN, YELLOW, RED, BLUE, CYAN, MAGENTA
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


class SampleStage(Enum):
    """Stages a sample goes through in the workflow."""

    PENDING = "Pending"
    DOWNLOADING = "Downloading"
    DOWNLOADED = "Downloaded"
    EXTRACTING = "Extracting"
    EXTRACTED = "Extracted"
    QUANTIFYING = "Quantifying"
    DONE = "Done"
    FAILED = "Failed"
    TRUNCATED = "Truncated"


@dataclass
class SampleState:
    """Tracks a sample's state through the workflow."""

    sample_id: str
    sra_url: str
    dest_path: Path
    stage: SampleStage = SampleStage.PENDING
    current_bytes: int = 0
    total_bytes: int = 0
    error: str = ""


class WorkflowManager:
    """Manages the complete RNA-seq workflow with TUI visualization."""

    STAGE_COLORS = {
        SampleStage.PENDING: BLUE,
        SampleStage.DOWNLOADING: CYAN,
        SampleStage.DOWNLOADED: YELLOW,
        SampleStage.EXTRACTING: MAGENTA,
        SampleStage.EXTRACTED: YELLOW,
        SampleStage.QUANTIFYING: CYAN,
        SampleStage.DONE: GREEN,
        SampleStage.FAILED: RED,
        SampleStage.TRUNCATED: RED,
    }

    def __init__(self, config_path: Path, max_threads: int = 5):
        self.config_path = Path(config_path)
        self.max_threads = max_threads
        self.ui = TerminalInterface()
        self.samples: Dict[str, SampleState] = {}
        self.executor = ThreadPoolExecutor(max_workers=max_threads)
        self._running = False

        # Load config to get paths
        self._load_config()

    def _load_config(self):
        """Load workflow configuration."""
        import yaml

        with open(self.config_path) as f:
            self.config = yaml.safe_load(f)

        self.work_dir = Path(self.config.get("work_dir", "output/amalgkit"))
        self.species = self.config.get("species", "unknown")

    def add_sample(self, sample_id: str, sra_url: str, dest_path: Path):
        """Add a sample to track through the workflow."""
        self.samples[sample_id] = SampleState(
            sample_id=sample_id,
            sra_url=sra_url,
            dest_path=dest_path,
        )
        self.ui.add_bar(sample_id, sample_id[:20], total=0.0, unit="MB")  # Truncate long IDs
        self.ui.update(sample_id, stage="Pending", status="Queued", color=BLUE)

    def run(self) -> Dict[str, bool]:
        """Execute the full workflow for all samples.

        Returns:
            Dict mapping sample_id to success status.
        """
        import logging

        # Silence verbose loggers during TUI
        robust_logger = logging.getLogger("metainformant.core.io.download_robust")
        original_level = robust_logger.getEffectiveLevel()
        robust_logger.setLevel(logging.ERROR)  # Suppress all INFO/WARNING during TUI

        results = {}
        self._running = True
        self.ui.start()

        try:
            # Phase 1: Download all SRA files
            self._run_download_phase()

            # Phase 2: Run getfastq for extracted samples
            self._run_getfastq_phase()

            # Phase 3: Run quant for all extracted samples
            self._run_quant_phase()

            # Collect results
            for sample_id, state in self.samples.items():
                results[sample_id] = state.stage == SampleStage.DONE

        finally:
            self.ui.stop()
            self.executor.shutdown(wait=True)
            robust_logger.setLevel(original_level)

        return results

    def _run_download_phase(self):
        """Download all SRA files in parallel."""
        from metainformant.core.io.download_robust import robust_download_url, get_remote_file_size

        futures = {}

        for sample_id, state in self.samples.items():
            state.stage = SampleStage.DOWNLOADING
            self.ui.update(sample_id, stage="Download", status="Starting...", color=CYAN)

            future = self.executor.submit(self._download_sample, sample_id, state.sra_url, state.dest_path)
            futures[future] = sample_id

        # Monitor progress
        while any(not f.done() for f in futures):
            self._update_download_progress()
            time.sleep(0.5)

        # Collect download results
        for f in as_completed(futures):
            sample_id = futures[f]
            try:
                success = f.result()
                if success:
                    self.samples[sample_id].stage = SampleStage.DOWNLOADED
                    self.ui.update(sample_id, stage="Download", status="Downloaded ✓", color=YELLOW)
                else:
                    self.samples[sample_id].stage = SampleStage.FAILED
            except Exception as e:
                self.samples[sample_id].stage = SampleStage.FAILED
                self.samples[sample_id].error = str(e)
                logger.error(f"Download failed for {sample_id}: {e}")

    def _download_sample(self, sample_id: str, url: str, dest: Path) -> bool:
        """Download a single sample."""
        from metainformant.core.io.download_robust import robust_download_url, get_remote_file_size

        tid = threading.get_native_id()
        self.ui.update(sample_id, status=f"Init (TID:{tid})", color=CYAN)

        # Get expected file size
        total_bytes = get_remote_file_size(url)
        total_mb = total_bytes / 1024 / 1024 if total_bytes else 0
        self.samples[sample_id].total_bytes = total_bytes

        # Check if file already exists with correct size (skip re-download)
        if dest.exists():
            existing_size = dest.stat().st_size
            existing_mb = existing_size / 1024 / 1024

            # If file exists and is >= 99% of expected size, skip download
            if total_bytes > 0 and existing_size >= (total_bytes * 0.99):
                self.ui.update(sample_id, current=existing_mb, total=existing_mb, status="Cached ✓", color=GREEN)
                return True
            elif total_bytes == 0:
                # Unknown size, trust existing file
                self.ui.update(sample_id, current=existing_mb, total=existing_mb, status="Cached ✓", color=GREEN)
                return True

        if total_mb > 0:
            self.ui.update(sample_id, total=total_mb, status=f"Downloading (TID:{tid})")
        else:
            self.ui.update(sample_id, status=f"Downloading (TID:{tid})")

        success = robust_download_url(url, dest)

        if success and dest.exists():
            final_size = dest.stat().st_size
            final_mb = final_size / 1024 / 1024

            # Check for truncation
            if total_bytes > 0 and final_size < (total_bytes * 0.99):
                self.ui.update(sample_id, status="Truncated", color=RED)
                self.samples[sample_id].stage = SampleStage.TRUNCATED
                return False

            self.ui.update(sample_id, current=final_mb, total=max(total_mb, final_mb))
            return True
        else:
            self.ui.update(sample_id, status="Failed", color=RED)
            return False

    def _update_download_progress(self):
        """Poll file sizes during download phase."""
        active = 0
        done = 0

        for sample_id, state in self.samples.items():
            if state.stage == SampleStage.DOWNLOADING:
                dest = state.dest_path
                temp = dest.with_suffix(dest.suffix + ".part")

                current_size = 0
                if dest.exists():
                    current_size = dest.stat().st_size
                    done += 1
                elif temp.exists():
                    current_size = temp.stat().st_size
                    active += 1

                self.ui.update(
                    sample_id, current=current_size / 1024 / 1024, speed=f"{current_size / 1024 / 1024:.1f} MB"
                )
            elif state.stage == SampleStage.DOWNLOADED:
                done += 1

        self.ui.set_footer(f"Phase: Download | Active: {active} | Done: {done} | Total: {len(self.samples)}")

    def _run_getfastq_phase(self):
        """Run getfastq (FASTQ extraction) for downloaded samples."""
        from metainformant.rna.amalgkit import run_amalgkit

        # Filter to downloaded samples
        downloaded = [s for s in self.samples.values() if s.stage == SampleStage.DOWNLOADED]

        for state in downloaded:
            state.stage = SampleStage.EXTRACTING
            self.ui.update(
                state.sample_id, stage="Extract", status="Extracting FASTQ...", color=MAGENTA, current=0.0, total=100.0
            )

        # Run amalgkit getfastq (runs on all samples at once)
        # getfastq handles extraction internally
        self.ui.set_footer(f"Phase: getfastq | Processing {len(downloaded)} samples...")

        try:
            # Build params from config - amalgkit getfastq accepts out_dir and metadata, NOT work_dir
            metadata_path = self.work_dir / "metadata" / "metadata_selected.tsv"

            params = {
                "out_dir": str(self.work_dir),
                "metadata": str(metadata_path),
                "threads": self.config.get("threads", 8),
            }

            result = run_amalgkit("getfastq", params)

            if result.returncode == 0:
                for state in downloaded:
                    state.stage = SampleStage.EXTRACTED
                    self.ui.update(
                        state.sample_id, current=100.0, total=100.0, stage="Extract", status="Extracted ✓", color=YELLOW
                    )
            else:
                for state in downloaded:
                    state.stage = SampleStage.FAILED
                    state.error = "getfastq failed"
                    self.ui.update(state.sample_id, status="Extract Failed", color=RED)

        except Exception as e:
            logger.error(f"getfastq phase failed: {e}")
            for state in downloaded:
                state.stage = SampleStage.FAILED
                state.error = str(e)
                self.ui.update(state.sample_id, status="Error", color=RED)

    def _run_quant_phase(self):
        """Run quantification for extracted samples."""
        from metainformant.rna.amalgkit import run_amalgkit

        # Filter to extracted samples
        extracted = [s for s in self.samples.values() if s.stage == SampleStage.EXTRACTED]

        for state in extracted:
            state.stage = SampleStage.QUANTIFYING
            self.ui.update(
                state.sample_id, stage="Quant", status="Quantifying...", color=CYAN, current=0.0, total=100.0
            )

        self.ui.set_footer(f"Phase: quant | Processing {len(extracted)} samples...")

        try:
            # Build params from config - amalgkit quant accepts out_dir and metadata, NOT work_dir
            metadata_path = self.work_dir / "metadata" / "metadata_selected.tsv"

            params = {
                "out_dir": str(self.work_dir),
                "metadata": str(metadata_path),
                "threads": self.config.get("threads", 8),
            }

            result = run_amalgkit("quant", params)

            if result.returncode == 0:
                for state in extracted:
                    state.stage = SampleStage.DONE
                    self.ui.update(
                        state.sample_id, current=100.0, total=100.0, stage="Done", status="Quantified ✓", color=GREEN
                    )
            else:
                for state in extracted:
                    state.stage = SampleStage.FAILED
                    state.error = "quant failed"
                    self.ui.update(state.sample_id, status="Quant Failed", color=RED)

        except Exception as e:
            logger.error(f"quant phase failed: {e}")
            for state in extracted:
                state.stage = SampleStage.FAILED
                state.error = str(e)
                self.ui.update(state.sample_id, status="Error", color=RED)

        # Final summary
        done = sum(1 for s in self.samples.values() if s.stage == SampleStage.DONE)
        failed = sum(1 for s in self.samples.values() if s.stage in (SampleStage.FAILED, SampleStage.TRUNCATED))
        self.ui.set_footer(f"✓ Complete | Done: {done} | Failed: {failed} | Total: {len(self.samples)}")
