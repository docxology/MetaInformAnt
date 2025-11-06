"""Real-time download progress tracking for getfastq step.

This module provides comprehensive progress monitoring for FASTQ downloads:
- Per-thread download progress tracking
- File size monitoring in real-time
- Download rate calculations (MB/s)
- Elapsed time tracking
- Live terminal updates with progress bars
"""

from __future__ import annotations

import logging
import threading
import time
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

from ...core.progress import progress_bar

logger = logging.getLogger(__name__)


class FileSizeMonitor:
    """Monitors file size changes in a directory for download progress."""
    
    def __init__(self, sample_dir: Path, update_interval: float = 2.0):
        """Initialize file size monitor.
        
        Args:
            sample_dir: Directory to monitor (e.g., out_dir/getfastq/SRR1234567/)
            update_interval: How often to check file sizes in seconds (default: 2.0)
        """
        self.sample_dir = Path(sample_dir)
        self.update_interval = update_interval
        self.last_size = 0
        self.current_size = 0
        self.last_update = time.time()
        self._lock = threading.Lock()
    
    def get_size(self) -> int:
        """Get current total size of all files in the directory.
        
        Returns:
            Total size in bytes, or 0 if directory doesn't exist
        """
        if not self.sample_dir.exists():
            return 0
        
        total_size = 0
        try:
            # Count both .sra files (during prefetch) and .fastq* files (after conversion)
            for file_path in self.sample_dir.rglob("*"):
                if file_path.is_file():
                    try:
                        total_size += file_path.stat().st_size
                    except (OSError, FileNotFoundError):
                        # File may have been deleted or moved
                        pass
        except Exception as e:
            logger.debug(f"Error calculating size for {self.sample_dir}: {e}")
        
        return total_size
    
    def update(self) -> tuple[int, float, float]:
        """Update size tracking and calculate rate.
        
        Returns:
            Tuple of (current_size_bytes, size_delta_bytes, rate_mbps)
            rate_mbps is 0.0 if not enough time has passed or size decreased
        """
        with self._lock:
            now = time.time()
            self.current_size = self.get_size()
            
            elapsed = now - self.last_update
            size_delta = self.current_size - self.last_size
            
            # Calculate rate (MB/s) if size increased and enough time passed
            rate_mbps = 0.0
            if elapsed > 0 and size_delta > 0:
                rate_mbps = (size_delta / (1024 * 1024)) / elapsed
            
            self.last_size = self.current_size
            self.last_update = now
            
            return self.current_size, size_delta, rate_mbps


class ThreadProgressTracker:
    """Tracks progress for a single download thread."""
    
    def __init__(
        self,
        thread_id: int,
        run_id: str,
        sample_dir: Path,
        update_interval: float = 2.0,
        use_progress_bar: bool = True,
    ):
        """Initialize thread progress tracker.
        
        Args:
            thread_id: Unique thread identifier
            run_id: SRA run ID being downloaded (e.g., "SRR1234567")
            sample_dir: Directory where files are being downloaded
            update_interval: How often to update progress (seconds)
            use_progress_bar: If True, use tqdm progress bar; otherwise use text updates
        """
        self.thread_id = thread_id
        self.run_id = run_id
        self.sample_dir = Path(sample_dir)
        self.update_interval = update_interval
        self.use_progress_bar = use_progress_bar
        
        self.file_monitor = FileSizeMonitor(sample_dir, update_interval)
        self.start_time = time.time()
        self.last_update_time = time.time()
        self.is_complete = False
        self.failed = False
        self._lock = threading.Lock()
        
        # Progress bar (optional)
        self.progress_bar = None
        if use_progress_bar:
            try:
                # Use tqdm if available
                self.progress_bar = progress_bar(
                    total=None,  # Indeterminate size
                    desc=f"Thread {thread_id}: {run_id}",
                    unit="B",
                    unit_scale=True,
                    leave=False,
                )
            except Exception:
                self.use_progress_bar = False
                logger.debug("Progress bar not available, using text updates")
    
    def update(self) -> dict[str, Any]:
        """Update progress tracking.
        
        Returns:
            Dictionary with current progress state:
            - size_bytes: Current total size in bytes
            - size_mb: Current total size in MB
            - rate_mbps: Current download rate in MB/s
            - elapsed_seconds: Time elapsed since start
            - is_complete: Whether download is complete
        """
        with self._lock:
            if self.is_complete:
                return self._get_state()
            
            size_bytes, size_delta, rate_mbps = self.file_monitor.update()
            elapsed = time.time() - self.start_time
            
            # Update progress bar if available
            if self.progress_bar:
                try:
                    # Update description with current stats
                    elapsed_str = str(timedelta(seconds=int(elapsed)))
                    desc = f"Thread {self.thread_id}: {self.run_id} [{size_bytes / (1024*1024):.1f}MB @ {rate_mbps:.2f}MB/s] ({elapsed_str})"
                    self.progress_bar.set_description(desc)
                    self.progress_bar.update(size_delta)
                except Exception as e:
                    logger.debug(f"Error updating progress bar: {e}")
            
            self.last_update_time = time.time()
            
            return {
                "size_bytes": size_bytes,
                "size_mb": size_bytes / (1024 * 1024),
                "rate_mbps": rate_mbps,
                "elapsed_seconds": elapsed,
                "is_complete": self.is_complete,
                "failed": self.failed,
            }
    
    def mark_complete(self, success: bool = True) -> None:
        """Mark download as complete.
        
        Args:
            success: If True, download succeeded; if False, it failed
        """
        with self._lock:
            self.is_complete = True
            self.failed = not success
            
            if self.progress_bar:
                try:
                    if success:
                        self.progress_bar.set_description(f"Thread {self.thread_id}: {self.run_id} [Complete]")
                    else:
                        self.progress_bar.set_description(f"Thread {self.thread_id}: {self.run_id} [Failed]")
                    self.progress_bar.close()
                except Exception:
                    pass
    
    def _get_state(self) -> dict[str, Any]:
        """Get current state without updating."""
        elapsed = time.time() - self.start_time
        size_bytes = self.file_monitor.current_size
        
        return {
            "size_bytes": size_bytes,
            "size_mb": size_bytes / (1024 * 1024),
            "rate_mbps": 0.0,
            "elapsed_seconds": elapsed,
            "is_complete": self.is_complete,
            "failed": self.failed,
        }
    
    def close(self) -> None:
        """Close progress bar and cleanup."""
        if self.progress_bar:
            try:
                self.progress_bar.close()
            except Exception:
                pass


class DownloadProgressMonitor:
    """Main progress monitor that tracks multiple concurrent downloads."""
    
    def __init__(
        self,
        out_dir: Path,
        update_interval: float = 2.0,
        use_progress_bars: bool = True,
        show_summary: bool = True,
    ):
        """Initialize download progress monitor.
        
        Args:
            out_dir: Base output directory (e.g., output/amalgkit/work/fastq)
            update_interval: How often to update progress (seconds, default: 2.0)
            use_progress_bars: If True, show tqdm progress bars; otherwise use text
            show_summary: If True, periodically print summary statistics
        """
        self.out_dir = Path(out_dir)
        self.update_interval = update_interval
        self.use_progress_bars = use_progress_bars
        self.show_summary = show_summary
        
        # Track active downloads: thread_id -> ThreadProgressTracker
        self.active_trackers: dict[int, ThreadProgressTracker] = {}
        self._lock = threading.Lock()
        
        # Summary statistics
        self.total_samples = 0
        self.completed_samples = 0
        self.failed_samples = 0
        self.start_time = time.time()
        
        # Background monitoring thread
        self._monitoring = False
        self._monitor_thread: threading.Thread | None = None
    
    def register_thread(
        self,
        thread_id: int,
        run_id: str,
        sample_dir: Path | None = None,
    ) -> None:
        """Register a new download thread for tracking.
        
        Args:
            thread_id: Unique thread identifier
            run_id: SRA run ID being downloaded
            sample_dir: Optional specific directory; if None, uses out_dir/getfastq/run_id
        """
        if sample_dir is None:
            sample_dir = self.out_dir / "getfastq" / run_id
        
        with self._lock:
            tracker = ThreadProgressTracker(
                thread_id=thread_id,
                run_id=run_id,
                sample_dir=sample_dir,
                update_interval=self.update_interval,
                use_progress_bar=self.use_progress_bars,
            )
            self.active_trackers[thread_id] = tracker
            self.total_samples += 1
    
    def unregister_thread(self, thread_id: int, success: bool = True) -> None:
        """Unregister a download thread (download complete).
        
        Args:
            thread_id: Thread identifier to unregister
            success: Whether download succeeded
        """
        with self._lock:
            if thread_id in self.active_trackers:
                tracker = self.active_trackers[thread_id]
                tracker.mark_complete(success)
                
                if success:
                    self.completed_samples += 1
                else:
                    self.failed_samples += 1
                
                # Don't remove immediately - keep for final summary
                # Will be cleaned up by stop_monitoring()
    
    def update_thread(self, thread_id: int) -> dict[str, Any] | None:
        """Update progress for a specific thread.
        
        Args:
            thread_id: Thread identifier
            
        Returns:
            Current progress state dict, or None if thread not found
        """
        with self._lock:
            if thread_id not in self.active_trackers:
                return None
            return self.active_trackers[thread_id].update()
    
    def get_summary(self) -> dict[str, Any]:
        """Get overall progress summary.
        
        Returns:
            Dictionary with summary statistics:
            - total_samples: Total number of samples
            - completed: Number completed
            - failed: Number failed
            - active: Number currently downloading
            - total_size_mb: Total size of all downloads
            - avg_rate_mbps: Average download rate
            - elapsed_seconds: Total elapsed time
        """
        with self._lock:
            active_count = len([t for t in self.active_trackers.values() if not t.is_complete])
            
            total_size = 0
            total_rate = 0.0
            active_rates = []
            
            for tracker in self.active_trackers.values():
                state = tracker._get_state()
                total_size += state["size_bytes"]
                if not tracker.is_complete and state["rate_mbps"] > 0:
                    active_rates.append(state["rate_mbps"])
            
            avg_rate = sum(active_rates) / len(active_rates) if active_rates else 0.0
            elapsed = time.time() - self.start_time
            
            return {
                "total_samples": self.total_samples,
                "completed": self.completed_samples,
                "failed": self.failed_samples,
                "active": active_count,
                "total_size_mb": total_size / (1024 * 1024),
                "avg_rate_mbps": avg_rate,
                "elapsed_seconds": elapsed,
            }
    
    def start_monitoring(self) -> None:
        """Start background monitoring thread that updates progress and shows summaries."""
        if self._monitoring:
            return
        
        self._monitoring = True
        
        def _monitor_loop():
            while self._monitoring:
                try:
                    # Update all active trackers
                    with self._lock:
                        for tracker in self.active_trackers.values():
                            if not tracker.is_complete:
                                tracker.update()
                    
                    # Show summary if enabled
                    if self.show_summary and not self.use_progress_bars:
                        summary = self.get_summary()
                        if summary["active"] > 0 or summary["completed"] < summary["total_samples"]:
                            elapsed_str = str(timedelta(seconds=int(summary["elapsed_seconds"])))
                            logger.info(
                                f"📊 Progress: {summary['completed']}/{summary['total_samples']} samples | "
                                f"{summary['total_size_mb']:.1f}MB downloaded | "
                                f"Avg: {summary['avg_rate_mbps']:.2f}MB/s | "
                                f"Elapsed: {elapsed_str} | "
                                f"Active: {summary['active']}"
                            )
                    
                    time.sleep(self.update_interval)
                except Exception as e:
                    logger.debug(f"Error in monitoring loop: {e}")
                    time.sleep(self.update_interval)
        
        self._monitor_thread = threading.Thread(target=_monitor_loop, daemon=True)
        self._monitor_thread.start()
    
    def stop_monitoring(self) -> None:
        """Stop background monitoring and cleanup."""
        self._monitoring = False
        
        if self._monitor_thread:
            self._monitor_thread.join(timeout=5.0)
        
        # Close all progress bars
        with self._lock:
            for tracker in self.active_trackers.values():
                tracker.close()
        
        # Print final summary
        summary = self.get_summary()
        elapsed_str = str(timedelta(seconds=int(summary["elapsed_seconds"])))
        logger.info("=" * 80)
        logger.info("DOWNLOAD PROGRESS SUMMARY")
        logger.info(f"Total samples: {summary['total_samples']}")
        logger.info(f"Completed: {summary['completed']}")
        logger.info(f"Failed: {summary['failed']}")
        logger.info(f"Total size: {summary['total_size_mb']:.1f} MB")
        logger.info(f"Average rate: {summary['avg_rate_mbps']:.2f} MB/s")
        logger.info(f"Total time: {elapsed_str}")
        logger.info("=" * 80)

