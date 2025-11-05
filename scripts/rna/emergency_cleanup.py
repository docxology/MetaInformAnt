#!/usr/bin/env python3
"""Emergency cleanup of quantified FASTQ files to free disk space.

This script finds all samples that have been quantified (have abundance.tsv)
but still have FASTQ files, and deletes the FASTQ files to free space.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.logging import get_logger
from metainformant.rna.progress_tracker import get_tracker

logger = get_logger("emergency_cleanup")


def cleanup_quantified_fastqs(dry_run: bool = False) -> dict[str, int]:
    """Delete FASTQ files for samples that have been quantified.
    
    Args:
        dry_run: If True, only report what would be deleted
        
    Returns:
        Dictionary with stats: {'deleted': count, 'freed_mb': size_in_mb}
    """
    output_dir = Path("output/amalgkit")
    stats = {"deleted": 0, "freed_mb": 0.0, "errors": 0}
    tracker = get_tracker()
    
    for species_dir in sorted(output_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        
        species_id = species_dir.name
        quant_dir = species_dir / "quant"
        fastq_dir = species_dir / "fastq"
        
        if not quant_dir.exists() or not fastq_dir.exists():
            continue
        
        logger.info(f"Checking {species_id}...")
        
        # Find quantified samples
        for sample_dir in fastq_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            
            sample_id = sample_dir.name
            abundance_file = quant_dir / sample_id / "abundance.tsv"
            
            if abundance_file.exists():
                # This sample is quantified - delete FASTQs
                fastq_files = list(sample_dir.glob("*.fastq.gz")) + list(sample_dir.glob("*.fastq"))
                
                if fastq_files:
                    # Calculate size before deletion
                    total_size = sum(f.stat().st_size for f in fastq_files)
                    size_mb = total_size / (1024 * 1024)
                    
                    if dry_run:
                        logger.info(f"  Would delete {sample_id}: {size_mb:.1f} MB")
                        stats["freed_mb"] += size_mb
                    else:
                        try:
                            # Delete the entire sample directory
                            import shutil
                            shutil.rmtree(sample_dir)
                            stats["deleted"] += 1
                            stats["freed_mb"] += size_mb
                            
                            # Update tracker
                            tracker.on_delete_complete(species_id, sample_id)
                            
                            logger.info(f"  ✅ Deleted {sample_id}: {size_mb:.1f} MB")
                        except Exception as e:
                            stats["errors"] += 1
                            logger.error(f"  ❌ Failed to delete {sample_id}: {e}")
    
    if not dry_run:
        tracker._save_state()
        tracker.update_dashboard()
    
    return stats


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Emergency cleanup of quantified FASTQ files")
    parser.add_argument("--dry-run", action="store_true", help="Only report what would be deleted")
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("EMERGENCY CLEANUP: Quantified FASTQ Files")
    logger.info("=" * 80)
    
    if args.dry_run:
        logger.info("DRY RUN MODE - no files will be deleted")
    else:
        logger.warning("DELETING FASTQ files for quantified samples!")
    
    logger.info("")
    
    stats = cleanup_quantified_fastqs(dry_run=args.dry_run)
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("CLEANUP SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Samples processed: {stats['deleted']}")
    logger.info(f"Space freed: {stats['freed_mb']:.1f} MB ({stats['freed_mb']/1024:.2f} GB)")
    if stats["errors"] > 0:
        logger.warning(f"Errors: {stats['errors']}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()

