#!/usr/bin/env python3
"""Unified SRA to FASTQ conversion script.

This script finds SRA files that need conversion to FASTQ format and converts them.
Supports both sequential and parallel processing modes.

Usage:
    # Sequential conversion (default)
    python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

    # Parallel conversion (faster for multiple samples)
    python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --parallel

    # Custom number of parallel workers
    python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --parallel --workers 5
"""

from __future__ import annotations

import argparse
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.steps.getfastq import convert_sra_to_fastq
from metainformant.rna.workflow import load_workflow_config
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def find_sra_files_needing_conversion(fastq_dir: Path) -> list[tuple[str, Path, Path]]:
    """Find SRA files that need conversion to FASTQ.

    Args:
        fastq_dir: Base directory containing getfastq subdirectories

    Returns:
        List of tuples: (sample_id, sra_file_path, output_dir)
    """
    sra_files_to_convert: list[tuple[str, Path, Path]] = []
    getfastq_dir = fastq_dir / "getfastq"

    if not getfastq_dir.exists():
        logger.warning(f"getfastq directory not found: {getfastq_dir}")
        return sra_files_to_convert

    # Find all SRA files
    for sra_file in getfastq_dir.rglob("*.sra"):
        sample_id = sra_file.stem  # e.g., "SRR14740513" from "SRR14740513.sra"
        output_dir = sra_file.parent

        # Check if FASTQ files already exist
        has_fastq = False
        for pattern in [
            f"{sample_id}_*.fastq.gz",
            f"{sample_id}_*.fastq",
            f"{sample_id}.fastq.gz",
            f"{sample_id}.fastq",
        ]:
            if list(output_dir.glob(pattern)):
                has_fastq = True
                break

        if not has_fastq:
            sra_files_to_convert.append((sample_id, sra_file, output_dir))
            logger.info(f"Found SRA file needing conversion: {sample_id} at {sra_file}")
        else:
            logger.debug(f"Skipping {sample_id}: FASTQ files already exist")

    return sra_files_to_convert


def convert_sequential(sra_files: list[tuple[str, Path, Path]], threads: int, log_dir: Path | None) -> tuple[int, int]:
    """Convert SRA files sequentially (one at a time)."""
    successful = 0
    failed = 0

    logger.info("Starting sequential conversion...")

    for idx, (sample_id, sra_file, output_dir) in enumerate(sra_files, 1):
        logger.info(f"[{idx}/{len(sra_files)}] Converting {sample_id}...")

        success, message, fastq_files = convert_sra_to_fastq(
            sample_id,
            sra_file,
            output_dir,
            threads=threads,
            log_dir=log_dir,
        )

        if success:
            successful += 1
            logger.info(f"  ✅ Success: {message}")
            if fastq_files:
                total_size_mb = sum(f.stat().st_size for f in fastq_files) / 1e6
                logger.info(f"  Created {len(fastq_files)} FASTQ files ({total_size_mb:.1f} MB total)")
        else:
            failed += 1
            logger.error(f"  ❌ Failed: {message}")
        logger.info("")

    return successful, failed


def convert_parallel(sra_files: list[tuple[str, Path, Path]], threads: int, workers: int, log_dir: Path | None) -> tuple[int, int]:
    """Convert SRA files in parallel using multiple workers."""
    logger.info(f"Starting parallel conversion with {workers} workers...")

    # Progress tracking
    successful = 0
    failed = 0
    start_time = time.time()
    completed = 0
    lock = threading.Lock()

    def convert_sample(sample_id: str, sra_file: Path, idx: int, total: int) -> tuple[str, bool, str, int]:
        """Convert a single sample and return results."""
        nonlocal completed, successful, failed

        # Verify SRA file still exists before starting
        if not sra_file.exists():
            logger.warning(f"   [{idx}/{total}] ⚠️  {sample_id}: SRA file not found, skipping")
            with lock:
                completed += 1
                failed += 1
            return sample_id, False, f"SRA file not found: {sra_file}", 0

        sample_start = time.time()
        logger.info(f"   [{idx}/{total}] 🔄 Starting conversion for {sample_id}...")

        success, message, fastq_files = convert_sra_to_fastq(
            sample_id,
            sra_file,
            sra_file.parent,
            threads=threads,
            log_dir=log_dir,
        )

        elapsed = time.time() - sample_start
        with lock:
            completed += 1

            if success and fastq_files:
                successful += 1
                total_elapsed = time.time() - start_time
                avg_time = total_elapsed / completed
                remaining = (total - completed) * avg_time
                logger.info(
                    f"   [{idx}/{total}] ✅ {sample_id}: Successfully converted "
                    f"({len(fastq_files)} files, {elapsed:.1f}s elapsed, "
                    f"~{remaining/60:.1f}min remaining)"
                )
                return sample_id, True, message, len(fastq_files)
            else:
                failed += 1
                total_elapsed = time.time() - start_time
                avg_time = total_elapsed / completed if completed > 0 else 0
                remaining = (total - completed) * avg_time if completed > 0 else 0
                logger.error(
                    f"   [{idx}/{total}] ❌ {sample_id}: Conversion failed: {message} "
                    f"({elapsed:.1f}s elapsed, ~{remaining/60:.1f}min remaining)"
                )
                return sample_id, False, message, 0

    # Process in parallel
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(convert_sample, sample_id, sra_file, idx + 1, len(sra_files))
            for idx, (sample_id, sra_file, _) in enumerate(sra_files)
        }

        # Wait for all to complete
        for future in as_completed(futures):
            try:
                sample_id, success, message, file_count = future.result()
            except Exception as e:
                logger.error(f"   💥 Exception during conversion: {e}", exc_info=True)
                with lock:
                    failed += 1
                    completed += 1

    return successful, failed


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Unified SRA to FASTQ conversion script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Sequential conversion (default)
  python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

  # Parallel conversion (faster for multiple samples)
  python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --parallel

  # Custom number of parallel workers
  python3 scripts/rna/convert_sra_to_fastq.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --parallel --workers 5
        """
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to amalgkit config YAML file",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing (default: sequential)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: min(10, num_samples))",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads per conversion (default: from config)",
    )

    args = parser.parse_args()

    # Load config
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1

    cfg = load_workflow_config(config_path)

    # Determine FASTQ directory
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))

    # Determine threads
    threads = args.threads or cfg.per_step.get("getfastq", {}).get("threads", 24)

    # Determine log directory
    log_dir = Path(cfg.work_dir) / "logs" if cfg.work_dir else None

    logger.info("=" * 80)
    logger.info("Unified SRA to FASTQ Conversion")
    logger.info("=" * 80)
    logger.info(f"Config: {config_path}")
    logger.info(f"FASTQ directory: {fastq_dir}")
    logger.info(f"Threads per conversion: {threads}")
    logger.info(f"Mode: {'Parallel' if args.parallel else 'Sequential'}")
    if args.parallel:
        num_workers = args.workers or min(10, 20)  # Will be updated based on actual samples
        logger.info(f"Parallel workers: {num_workers}")
    logger.info("")

    # Find SRA files needing conversion
    sra_files = find_sra_files_needing_conversion(fastq_dir)

    if not sra_files:
        logger.info("✅ No SRA files found that need conversion.")
        return 0

    logger.info(f"Found {len(sra_files)} SRA files needing conversion:")
    for sample_id, sra_file, _ in sra_files:
        sra_size_gb = sra_file.stat().st_size / 1e9
        logger.info(f"  - {sample_id}: {sra_file.name} ({sra_size_gb:.2f} GB)")
    logger.info("")

    # Convert files
    start_time = time.time()

    if args.parallel:
        # Parallel mode
        num_workers = args.workers or min(10, len(sra_files))
        successful, failed = convert_parallel(sra_files, threads, num_workers, log_dir)
    else:
        # Sequential mode (default)
        successful, failed = convert_sequential(sra_files, threads, log_dir)

    # Summary
    total_time = time.time() - start_time
    logger.info("=" * 80)
    logger.info("Conversion Summary")
    logger.info("=" * 80)
    logger.info(f"Total: {len(sra_files)}")
    logger.info(f"✅ Successful: {successful}")
    if failed > 0:
        logger.info(f"❌ Failed: {failed}")
    logger.info(f"⏱️  Total time: {total_time/60:.1f} minutes ({total_time:.1f} seconds)")
    if successful > 0:
        logger.info(f"📈 Average time per sample: {total_time/successful:.1f} seconds")
    logger.info("=" * 80)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

