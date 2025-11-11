#!/usr/bin/env python3
"""Convert existing SRA files to FASTQ format.

This script finds all SRA files that don't have corresponding FASTQ files
and automatically converts them using convert_sra_to_fastq().

Supports parallel processing for faster conversion.
"""

from __future__ import annotations

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

# Suppress optional dependency warnings for this script
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="metainformant.rna.discovery")

logger = get_logger(__name__)


def main() -> int:
    """Convert existing SRA files to FASTQ."""
    config_path = Path("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
    cfg = load_workflow_config(config_path)
    
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    getfastq_dir = fastq_dir / "getfastq"
    if not getfastq_dir.exists():
        getfastq_dir = fastq_dir
    
    log_dir = Path(cfg.work_dir) / "logs" if cfg.work_dir else None
    
    # Find samples with SRA but no FASTQ
    conversion_needed = []
    for sample_dir in getfastq_dir.iterdir():
        if sample_dir.is_dir() and sample_dir.name.startswith("SRR"):
            sra_files = list(sample_dir.glob("*.sra"))
            fastq_files = list(sample_dir.glob("*.fastq*"))
            
            # Only add if SRA file exists and has no FASTQ files
            if sra_files and not fastq_files:
                sra_file = sra_files[0]
                # Double-check file still exists (might have been deleted)
                if sra_file.exists() and sra_file.is_file():
                    conversion_needed.append((sample_dir.name, sra_file))
                else:
                    logger.warning(f"   ⚠️  {sample_dir.name}: SRA file listed but not found, skipping")
    
    if not conversion_needed:
        logger.info("✅ No SRA files need conversion (all have FASTQ or are empty)")
        return 0
    
    logger.info(f"🔄 Found {len(conversion_needed)} samples with SRA but no FASTQ")
    
    # Get configuration
    conversion_threads = cfg.per_step.get("getfastq", {}).get("threads", 24)
    num_workers = min(10, len(conversion_needed))  # Process up to 10 samples in parallel
    
    logger.info(f"   Converting with {num_workers} parallel workers")
    logger.info(f"   Each conversion uses {conversion_threads} threads")
    
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
            threads=conversion_threads,
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
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(convert_sample, sample_id, sra_file, idx + 1, len(conversion_needed))
            for idx, (sample_id, sra_file) in enumerate(conversion_needed)
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
    
    total_time = time.time() - start_time
    logger.info(f"\n📊 Conversion Summary:")
    logger.info(f"   ✅ Successful: {successful}")
    logger.info(f"   ❌ Failed: {failed}")
    logger.info(f"   📦 Total: {len(conversion_needed)}")
    logger.info(f"   ⏱️  Total time: {total_time/60:.1f} minutes ({total_time:.1f} seconds)")
    if successful > 0:
        logger.info(f"   📈 Average time per sample: {total_time/successful:.1f} seconds")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

