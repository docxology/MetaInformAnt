"""Batch processing utilities for long-read sequencing data.

Provides concurrent processing of multiple samples/files with progress
tracking, error handling, and result aggregation. Uses the core
parallel execution infrastructure.

Optional dependencies:
    - concurrent.futures (stdlib): For thread-based parallelism
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_ENV_PREFIX = "LR_"


@dataclass
class BatchResult:
    """Result from a batch processing operation.

    Attributes:
        total_items: Total number of items processed.
        successful: Number of successfully processed items.
        failed: Number of failed items.
        results: Dict mapping item IDs to their results.
        errors: Dict mapping item IDs to error messages.
        duration_seconds: Total elapsed time.
        items_per_second: Processing throughput.
    """

    total_items: int = 0
    successful: int = 0
    failed: int = 0
    results: dict[str, Any] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)
    duration_seconds: float = 0.0
    items_per_second: float = 0.0

    @property
    def success_rate(self) -> float:
        """Fraction of items processed successfully."""
        return self.successful / self.total_items if self.total_items > 0 else 0.0


def process_batch(
    items: dict[str, Any],
    processor: Callable[[str, Any], Any],
    max_workers: int | None = None,
    continue_on_error: bool = True,
    timeout_per_item: float | None = None,
) -> BatchResult:
    """Process multiple items using a provided function.

    Runs the processor function on each item, optionally in parallel
    using a thread pool. Collects results and errors, tracks timing.

    Args:
        items: Dict mapping item IDs to item data.
        processor: Function taking (item_id, item_data) -> result.
        max_workers: Maximum parallel workers. None = sequential.
            Set from LR_BATCH_WORKERS env var if not specified.
        continue_on_error: If True, continue processing after failures.
        timeout_per_item: Maximum seconds per item (None = no timeout).

    Returns:
        BatchResult with aggregated results and statistics.
    """
    if max_workers is None:
        env_workers = os.environ.get(f"{_ENV_PREFIX}BATCH_WORKERS")
        if env_workers:
            max_workers = int(env_workers)

    start_time = time.monotonic()
    result = BatchResult(total_items=len(items))

    if max_workers is not None and max_workers > 1:
        result = _process_parallel(items, processor, max_workers, continue_on_error)
    else:
        result = _process_sequential(items, processor, continue_on_error)

    result.duration_seconds = time.monotonic() - start_time
    result.items_per_second = (
        result.total_items / result.duration_seconds if result.duration_seconds > 0 else 0.0
    )

    logger.info(
        f"Batch processing complete: {result.successful}/{result.total_items} succeeded "
        f"({result.success_rate:.1%}) in {result.duration_seconds:.1f}s"
    )
    return result


def _process_sequential(
    items: dict[str, Any],
    processor: Callable[[str, Any], Any],
    continue_on_error: bool,
) -> BatchResult:
    """Process items sequentially."""
    result = BatchResult(total_items=len(items))

    for item_id, item_data in items.items():
        try:
            item_result = processor(item_id, item_data)
            result.results[item_id] = item_result
            result.successful += 1
            logger.debug(f"Processed item '{item_id}' successfully")
        except Exception as e:
            result.errors[item_id] = str(e)
            result.failed += 1
            logger.warning(f"Failed to process item '{item_id}': {e}")
            if not continue_on_error:
                break

    return result


def _process_parallel(
    items: dict[str, Any],
    processor: Callable[[str, Any], Any],
    max_workers: int,
    continue_on_error: bool,
) -> BatchResult:
    """Process items in parallel using thread pool."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    result = BatchResult(total_items=len(items))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_id = {
            executor.submit(processor, item_id, item_data): item_id
            for item_id, item_data in items.items()
        }

        for future in as_completed(future_to_id):
            item_id = future_to_id[future]
            try:
                item_result = future.result()
                result.results[item_id] = item_result
                result.successful += 1
            except Exception as e:
                result.errors[item_id] = str(e)
                result.failed += 1
                logger.warning(f"Failed to process item '{item_id}': {e}")
                if not continue_on_error:
                    executor.shutdown(wait=False, cancel_futures=True)
                    break

    return result


def batch_filter_reads(
    read_sets: dict[str, list[dict[str, Any]]],
    min_length: int = 1000,
    min_quality: float = 7.0,
    max_workers: int | None = None,
) -> BatchResult:
    """Filter reads across multiple samples in batch.

    Convenience function that applies standard QC filtering to multiple
    read sets simultaneously.

    Args:
        read_sets: Dict mapping sample IDs to lists of read dicts.
            Each read dict must have 'sequence' and optionally 'quality_string'.
        min_length: Minimum read length.
        min_quality: Minimum mean quality score.
        max_workers: Parallel workers.

    Returns:
        BatchResult with filtered read lists per sample.
    """
    from metainformant.longread.quality.filtering import ReadRecord, filter_by_length, filter_by_quality

    def _filter_sample(sample_id: str, reads: Any) -> dict[str, Any]:
        records = [
            ReadRecord(
                read_id=r.get("read_id", f"read_{i}"),
                sequence=r["sequence"],
                quality_string=r.get("quality_string", ""),
            )
            for i, r in enumerate(reads)
        ]

        length_filtered = filter_by_length(records, min_length=min_length)
        if min_quality > 0 and any(r.quality_string for r in length_filtered):
            quality_filtered = filter_by_quality(length_filtered, min_q=min_quality)
        else:
            quality_filtered = length_filtered

        return {
            "input_count": len(records),
            "output_count": len(quality_filtered),
            "pass_rate": len(quality_filtered) / len(records) if records else 0.0,
            "reads": quality_filtered,
        }

    return process_batch(read_sets, _filter_sample, max_workers=max_workers)


def batch_compute_metrics(
    read_sets: dict[str, list[int]],
    max_workers: int | None = None,
) -> BatchResult:
    """Compute QC metrics across multiple samples in batch.

    Args:
        read_sets: Dict mapping sample IDs to lists of read lengths.
        max_workers: Parallel workers.

    Returns:
        BatchResult with ReadLengthStatistics per sample.
    """
    from metainformant.longread.quality.metrics import read_length_stats

    def _compute_metrics(sample_id: str, lengths: Any) -> Any:
        return read_length_stats(lengths)

    return process_batch(read_sets, _compute_metrics, max_workers=max_workers)


__all__ = [
    "BatchResult",
    "process_batch",
    "batch_filter_reads",
    "batch_compute_metrics",
]
