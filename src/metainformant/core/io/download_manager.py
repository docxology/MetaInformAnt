"""Download Manager module.

Orchestrates parallel downloads and updates the Terminal UI.
"""

import threading
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional, Tuple

from metainformant.core.ui.tui import TerminalInterface, GREEN, YELLOW, RED, BLUE
from metainformant.core.io.download_robust import robust_download_url
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


class DownloadManager:
    """Manages parallel file downloads with TUI visualization."""

    def __init__(self, max_threads: int = 4):
        self.max_threads = max_threads
        self.ui = TerminalInterface()
        self.executor = ThreadPoolExecutor(max_workers=max_threads)
        self.tasks: List[Dict] = []

    def add_download(self, url: str, dest_path: Path, label: str):
        """Queue a download task."""
        self.tasks.append({"url": url, "dest": dest_path, "label": label, "id": f"task_{len(self.tasks)}"})
        self.ui.add_bar(f"task_{len(self.tasks)-1}", label, total=0.0, unit="%")

    def start(self) -> Dict[str, bool]:
        """Execute all queued downloads and block until complete."""
        import logging

        # Silence the robust logger to prevent console interference with TUI
        robust_logger = logging.getLogger("metainformant.core.io.download_robust")
        original_level = robust_logger.getEffectiveLevel()
        robust_logger.setLevel(logging.WARNING)

        results = {}
        futures = {}

        self.ui.start()

        try:
            for task in self.tasks:
                task_id = task["id"]
                dest = task["dest"]
                url = task["url"]

                future = self.executor.submit(self._download_worker, task_id, url, dest)
                futures[future] = task["label"]

            # Monitor loop for progress
            while any(not f.done() for f in futures):
                self._update_progress()
                time.sleep(0.5)

            # Collect results
            for f in as_completed(futures):
                label = futures[f]
                try:
                    success = f.result()
                    results[label] = success
                except Exception as e:
                    logger.error(f"Download {label} failed: {e}")
                    results[label] = False

        finally:
            self.ui.stop()
            self.executor.shutdown()
            # Restore logging
            robust_logger.setLevel(original_level)

        return results

    def _download_worker(self, task_id: str, url: str, dest: Path) -> bool:
        """Worker thread for a single download."""

        # 1. Fetch size first for progress bar
        from metainformant.core.io.download_robust import get_remote_file_size

        tid = threading.get_native_id()
        self.ui.update(task_id, status=f"Init (TID:{tid})", color=BLUE)

        total_bytes = get_remote_file_size(url)
        total_mb = total_bytes / 1024 / 1024 if total_bytes else 0

        if total_mb > 0:
            self.ui.update(task_id, total=total_mb, status=f"Downloading (TID:{tid})")
        else:
            self.ui.update(task_id, status=f"Downloading (TID:{tid}) [Size Unknown]")

        try:
            success = robust_download_url(url, dest)

            if success:
                # If we didn't know total, set it now to match current
                final_size_mb = 0
                if dest.exists():
                    final_size_mb = dest.stat().st_size / 1024 / 1024

                # Verify size if we knew it
                if total_mb > 0 and final_size_mb < (total_mb * 0.99):  # Allow 1% variance
                    self.ui.update(task_id, status="Truncated", color=RED)
                    logger.error(f"Download {task_id} truncated: {final_size_mb:.1f}MB < {total_mb:.1f}MB")
                    return False

                self.ui.update(
                    task_id, current=final_size_mb, total=max(total_mb, final_size_mb), status="Done", color=GREEN
                )
            else:
                self.ui.update(task_id, status="Failed", color=RED)
            return success

        except Exception as e:
            self.ui.update(task_id, status="Error", color=RED)
            logger.error(f"Worker error: {e}")
            return False

    def _update_progress(self):
        """Poll file sizes to update UI speed/activity."""

        active_count = 0
        done_count = 0

        for task in self.tasks:
            task_id = task["id"]
            dest = task["dest"]
            temp_part = dest.with_suffix(dest.suffix + ".part")

            current_size = 0
            if dest.exists():
                done_count += 1
                current_size = dest.stat().st_size
            elif temp_part.exists():
                active_count += 1
                current_size = temp_part.stat().st_size

            self.ui.update(
                task_id, current=current_size / 1024 / 1024, speed=f"{current_size / 1024 / 1024:.1f} MB"  # MB
            )

        # Update Footer
        self.ui.set_footer(f"Active Threads: {self.max_threads} | Running: {active_count} | Done: {done_count}")

    def start(self) -> Dict[str, bool]:
        """Execute all queued downloads and block until complete."""
        import logging

        # Silence the robust logger to prevent console interference with TUI
        robust_logger = logging.getLogger("metainformant.core.io.download_robust")
        original_level = robust_logger.getEffectiveLevel()
        robust_logger.setLevel(logging.WARNING)

        results = {}
        futures = {}

        self.ui.start()

        try:
            for task in self.tasks:
                task_id = task["id"]
                dest = task["dest"]
                url = task["url"]

                # Wrap the robust download to provide updates
                # Note: valid robust_download logic blocks, so we do it in threads.
                # However, robust_download_url uses subprocess.run which blocks.
                # To get progress we have to approximate or file-poll.
                # For this version, we will poll the file size in a monitor loop.

                future = self.executor.submit(self._download_worker, task_id, url, dest)
                futures[future] = task["label"]

            # Monitor loop for progress
            while any(not f.done() for f in futures):
                self._update_progress()
                time.sleep(0.5)

            # Collect results
            for f in as_completed(futures):
                label = futures[f]
                try:
                    success = f.result()
                    results[label] = success
                except Exception as e:
                    logger.error(f"Download {label} failed: {e}")
                    results[label] = False

        finally:
            self.ui.stop()
            self.executor.shutdown()
            # Restore logging
            robust_logger.setLevel(original_level)

        return results

    def _update_progress(self):
        """Poll file sizes to update UI speed/activity."""

        for task in self.tasks:
            task_id = task["id"]
            dest = task["dest"]
            temp_part = dest.with_suffix(dest.suffix + ".part")

            current_size = 0
            if dest.exists():
                current_size = dest.stat().st_size
            elif temp_part.exists():
                current_size = temp_part.stat().st_size

            # Update Total if not set (lazy init if not fetched yet)
            # The worker attempts to set it, but we can double check or just rely on state

            # Calculate speed (rough)
            # For now, just showing size

            self.ui.update(
                task_id,
                current=current_size / 1024 / 1024,  # MB
                # Total is managed by worker setting it
                speed=f"{current_size / 1024 / 1024:.1f} MB",
            )

        try:
            success = robust_download_url(url, dest)

            if success:
                final_size_mb = 0
                if dest.exists():
                    final_size_mb = dest.stat().st_size / 1024 / 1024

                # Verify size if we knew it
                if total_mb > 0 and final_size_mb < (total_mb * 0.99):  # Allow 1% variance
                    self.ui.update(task_id, status="Truncated", color=RED)
                    logger.error(f"Download {task_id} truncated: {final_size_mb:.1f}MB < {total_mb:.1f}MB")
                    return False

                # Update final state
                self.ui.update(
                    task_id, current=final_size_mb, total=max(total_mb, final_size_mb), status="Done", color=GREEN
                )
            else:
                self.ui.update(task_id, status="Failed", color=RED)
            return success

        except Exception as e:
            self.ui.update(task_id, status="Error", color=RED)
            logger.error(f"Worker error: {e}")
            return False

    def _update_progress(self):
        """Poll file sizes to update UI speed/activity."""

        active_count = 0
        done_count = 0
        failed_count = 0
        total_speed_mb = 0.0  # Placeholder for future

        for task in self.tasks:
            task_id = task["id"]
            dest = task["dest"]
            temp_part = dest.with_suffix(dest.suffix + ".part")

            # Check status from UI state to count
            # Use lock to access robustly? Or just infer from files
            # Better to use UI state if accessible, but files work

            if dest.exists():
                done_count += 1
                current_size = dest.stat().st_size
                # If we are done, no need to update continuous progress unless needed
            elif temp_part.exists():
                active_count += 1
                current_size = temp_part.stat().st_size
            else:
                # Pending or just starting
                current_size = 0

            # Update Total if current > total (in case header was wrong)
            # This requires access to current total in UI, or we blind update
            # Since we can't easily read back from UI here without lock, we just push update

            self.ui.update(
                task_id,
                current=current_size / 1024 / 1024,  # MB
                # If current > total, UI should ideally respect that.
                # We can't update total here easily without knowing it.
                # But TUI display is clamped now.
                speed=f"{current_size / 1024 / 1024:.1f} MB",
            )

        # Update Footer
        self.ui.set_footer(f"Active Threads: {self.max_threads} | Running: {active_count} | Done: {done_count}")
