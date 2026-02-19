"""
Process Watchdog Utility.

This module provides the `ProcessWatchdog` class, designed to monitor specific
processes (or process trees) and trigger callbacks if they stall (e.g., low CPU
usage for an extended period).

Usage:
    watchdog = ProcessWatchdog(pid=1234, cpu_threshold=1.0, timeout_seconds=3600)
    watchdog.start()
"""

import logging
import threading
import time
from typing import Callable, Optional

import psutil

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


class ProcessWatchdog:
    """
    Monitors a process for stalls based on CPU usage.
    Running in a separate thread.
    """

    def __init__(
        self,
        pid: int,
        cpu_threshold: float = 1.0,
        timeout_seconds: int = 3600,
        check_interval: int = 60,
        on_stall: Optional[Callable[[int], None]] = None,
    ):
        """
        Initialize the Watchdog.

        Args:
            pid: Process ID to monitor.
            cpu_threshold: Minimum average CPU percent required to be considered "active".
            timeout_seconds: Duration (in seconds) of low activity before triggering on_stall.
            check_interval: How often (in seconds) to query process status.
            on_stall: Callback function to execute when a stall is detected.
                      Receives the PID as an argument.
        """
        self.pid = pid
        self.cpu_threshold = cpu_threshold
        self.timeout_seconds = timeout_seconds
        self.check_interval = check_interval
        self.on_stall = on_stall
        
        self.low_activity_start_time: Optional[float] = None
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def start(self):
        """Start the monitoring thread."""
        if self._thread and self._thread.is_alive():
            logger.warning(f"Watchdog for PID {self.pid} is already running.")
            return

        self._stop_event.clear()
        self._thread = threading.Thread(target=self._monitor_loop, daemon=True, name=f"Watchdog-{self.pid}")
        self._thread.start()
        logger.info(f"Watchdog started for PID {self.pid} (Threshold: <{self.cpu_threshold}% CPU for {self.timeout_seconds}s)")

    def stop(self):
        """Stop the monitoring thread."""
        if self._thread and self._thread.is_alive():
            self._stop_event.set()
            self._thread.join(timeout=5)
            logger.info(f"Watchdog stopped for PID {self.pid}")

    @staticmethod
    def kill_process_tree(pid: int):
        """Kill a process and its children."""
        try:
            parent = psutil.Process(pid)
            children = parent.children(recursive=True)
            for child in children:
                try:
                    child.kill()
                except psutil.NoSuchProcess:
                    pass
            parent.kill()
            logger.info(f"Killed process tree for PID {pid}")
        except psutil.NoSuchProcess:
            pass

    def _monitor_loop(self):
        """Internal monitoring loop."""
        try:
            process = psutil.Process(self.pid)
        except psutil.NoSuchProcess:
            logger.debug(f"Watchdog failed to start: PID {self.pid} not found.")
            return

        while not self._stop_event.is_set():
            try:
                # Check if process is still running
                if not process.is_running():
                    break
                
                # Get CPU usage (non-blocking)
                # psutil.cpu_percent with interval=None returns usage since last call
                # We need an interval to get an instant reading, but blocking the thread is okay here
                # (interval=1.0)
                try:
                    cpu_pct = process.cpu_percent(interval=1.0)
                except psutil.NoSuchProcess:
                    break

                # Check activity
                if cpu_pct < self.cpu_threshold:
                    if self.low_activity_start_time is None:
                        self.low_activity_start_time = time.time()
                    
                    elapsed = time.time() - self.low_activity_start_time
                    if elapsed > self.timeout_seconds:
                        logger.warning(f"STALL DETECTED: PID {self.pid} has {cpu_pct}% CPU for {elapsed:.1f}s.")
                        if self.on_stall:
                            try:
                                self.on_stall(self.pid)
                            except Exception as e:
                                logger.error(f"Error in watchdog on_stall callback: {e}")
                        
                        # Monitor stopped after firing once? 
                        # Or reset? If we killed it, the loop will break next iteration.
                        # If we just warned, we reset.
                        self.low_activity_start_time = time.time() 
                else:
                    # Healthy activity
                    if self.low_activity_start_time is not None:
                        # logger.debug(f"PID {self.pid} recovered activity ({cpu_pct}% CPU).")
                        self.low_activity_start_time = None

            except psutil.NoSuchProcess:
                break
            except Exception as e:
                logger.error(f"Watchdog monitoring error: {e}")
                time.sleep(5)

            # Wait for next check
            if self._stop_event.wait(timeout=self.check_interval):
                break
