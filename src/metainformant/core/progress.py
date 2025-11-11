"""Progress tracking utilities for METAINFORMANT.

Provides progress bars and task tracking for long-running operations.
"""

from __future__ import annotations

import contextlib
from typing import Any, Iterator

from .logging import get_logger

logger = get_logger(__name__)

# Try to import tqdm, but make it optional
try:
    from tqdm import tqdm

    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False
    # Create a minimal tqdm-like interface
    class tqdm:  # type: ignore
        def __init__(self, iterable=None, total=None, desc=None, **kwargs):
            self.iterable = iterable
            self.total = total
            self.desc = desc or ""
            self.n = 0

        def __enter__(self):
            return self

        def __exit__(self, *args):
            if self.desc:
                logger.info(f"{self.desc}: completed {self.n} items")

        def __iter__(self):
            if self.iterable:
                for item in self.iterable:
                    self.n += 1
                    yield item
            else:
                return self

        def update(self, n=1):
            self.n += n

        def set_description(self, desc):
            self.desc = desc


def progress_bar(iterable: Iterator[Any] | None = None, total: int | None = None, desc: str | None = None, **kwargs: Any) -> tqdm:
    """Create a progress bar for iterating over items.

    Args:
        iterable: Iterable to wrap with progress bar
        total: Total number of items (if not inferrable from iterable)
        desc: Description text for progress bar
        **kwargs: Additional arguments for tqdm

    Returns:
        Progress bar object (tqdm-like interface)

    Example:
        for item in progress_bar(items, desc="Processing"):
            process(item)
    """
    if not _HAS_TQDM:
        logger.debug("tqdm not available, using basic progress tracking")
    return tqdm(iterable=iterable, total=total, desc=desc, **kwargs)


@contextlib.contextmanager
def task_context(task_name: str, total_steps: int | None = None) -> Iterator[Any]:
    """Context manager for tracking multi-step tasks.

    Args:
        task_name: Name of the task
        total_steps: Total number of steps (optional)

    Yields:
        Task tracker object

    Example:
        with task_context("Processing data", total_steps=10) as task:
            for i in range(10):
                process_step(i)
                task.update(1)
    """
    logger.info(f"Starting task: {task_name}" + (f" ({total_steps} steps)" if total_steps else ""))

    class TaskTracker:
        def __init__(self, name: str, total: int | None):
            self.name = name
            self.total = total
            self.current = 0

        def update(self, n: int = 1):
            self.current += n
            if self.total:
                pct = (self.current / self.total) * 100
                logger.debug(f"{self.name}: {self.current}/{self.total} ({pct:.1f}%)")
            else:
                logger.debug(f"{self.name}: {self.current} steps completed")

        def set_description(self, desc: str):
            logger.info(f"{self.name}: {desc}")

    tracker = TaskTracker(task_name, total_steps)
    try:
        yield tracker
        logger.info(f"Completed task: {task_name}")
    except Exception as e:
        logger.error(f"Task failed: {task_name}: {e}")
        raise


def log_progress(current: int, total: int | None, message: str = "") -> None:
    """Log progress for an operation.

    Args:
        current: Current progress count
        total: Total count (None for indeterminate)
        message: Optional message to include

    Example:
        log_progress(5, 10, "Processing items")
    """
    if total:
        pct = (current / total) * 100
        logger.info(f"{message}: {current}/{total} ({pct:.1f}%)")
    else:
        logger.info(f"{message}: {current} items")


