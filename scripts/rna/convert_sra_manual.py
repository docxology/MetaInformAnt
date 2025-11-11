#!/usr/bin/env python3
"""Manually convert SRA files to FASTQ format.

This script finds SRA files that have been downloaded but not yet converted
to FASTQ format and converts them using the convert_sra_to_fastq function.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.config import load_mapping_from_file
from metainformant.core.logging import get_logger
from metainformant.rna.steps.getfastq import convert_sra_to_fastq

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


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Manually convert SRA files to FASTQ format"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to amalgkit config YAML file",
    )
    parser.add_argument(
        "--fastq-dir",
        type=str,
        help="Override fastq directory from config (default: from config)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads for conversion (default: from config)",
    )
    args = parser.parse_args()
    
    # Load config
    config = load_mapping_from_file(args.config)
    getfastq_params = config.get("steps", {}).get("getfastq", {})
    
    # Determine fastq directory
    if args.fastq_dir:
        fastq_dir = Path(args.fastq_dir).expanduser().resolve()
    else:
        fastq_dir = Path(getfastq_params.get("out_dir", "output/amalgkit/pogonomyrmex_barbatus/fastq")).expanduser().resolve()
    
    # Determine threads
    threads = args.threads or getfastq_params.get("threads", 24)
    
    # Determine log directory
    log_dir = Path(config.get("log_dir", "output/amalgkit/pogonomyrmex_barbatus/logs")).expanduser().resolve()
    
    logger.info("=" * 80)
    logger.info("SRA to FASTQ Manual Conversion")
    logger.info("=" * 80)
    logger.info(f"FastQ directory: {fastq_dir}")
    logger.info(f"Threads: {threads}")
    logger.info(f"Log directory: {log_dir}")
    logger.info("")
    
    # Find SRA files needing conversion
    sra_files = find_sra_files_needing_conversion(fastq_dir)
    
    if not sra_files:
        logger.info("No SRA files found that need conversion.")
        return 0
    
    logger.info(f"Found {len(sra_files)} SRA files needing conversion:")
    for sample_id, sra_file, _ in sra_files:
        sra_size_gb = sra_file.stat().st_size / 1e9
        logger.info(f"  - {sample_id}: {sra_file.name} ({sra_size_gb:.2f} GB)")
    logger.info("")
    
    # Convert each SRA file
    successful = 0
    failed = 0
    
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
    
    # Summary
    logger.info("=" * 80)
    logger.info("Conversion Summary")
    logger.info("=" * 80)
    logger.info(f"Total: {len(sra_files)}")
    logger.info(f"✅ Successful: {successful}")
    if failed > 0:
        logger.info(f"❌ Failed: {failed}")
    logger.info("=" * 80)
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

