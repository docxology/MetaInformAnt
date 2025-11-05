#!/usr/bin/env python3
"""
Clean up partial and failed downloads.

Removes samples that have partial FASTQ/SRA files but are not quantified.
This frees up disk space for retrying downloads.
"""

import os
import sys
import shutil
from pathlib import Path
from datetime import datetime

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.core.logging import get_logger

logger = get_logger("cleanup_partial")


def find_partial_downloads(fastq_dir: Path, quant_dir: Path) -> list[tuple[str, Path, int]]:
    """Find samples with partial downloads that aren't quantified.
    
    Returns:
        List of (sample_id, sample_dir, size_mb) tuples
    """
    partial_samples = []
    
    if not fastq_dir.exists():
        return []
    
    # Check getfastq subdirectory
    for sample_dir in fastq_dir.glob("getfastq/SRR*"):
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        
        # Check if quantified
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if abundance_file.exists():
            continue
        
        # Check for files
        has_files = False
        for pattern in ["*.fastq*", "*.sra"]:
            if list(sample_dir.glob(pattern)):
                has_files = True
                break
        
        if has_files:
            # Calculate size
            size_mb = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file()) / (1024 * 1024)
            partial_samples.append((sample_id, sample_dir, int(size_mb)))
    
    # Check direct structure
    for sample_dir in fastq_dir.glob("SRR*"):
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        
        # Check if quantified
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if abundance_file.exists():
            continue
        
        # Check for files
        has_files = False
        for pattern in ["*.fastq*", "*.sra"]:
            if list(sample_dir.glob(pattern)):
                has_files = True
                break
        
        if has_files:
            # Calculate size
            size_mb = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file()) / (1024 * 1024)
            partial_samples.append((sample_id, sample_dir, int(size_mb)))
    
    return partial_samples


def cleanup_species(species_name: str, config_path: Path, dry_run: bool = False) -> dict:
    """Clean up partial downloads for a species."""
    try:
        cfg = load_workflow_config(config_path)
    except Exception as e:
        logger.warning(f"Failed to load config for {species_name}: {e}")
        return {"deleted": 0, "freed_mb": 0, "errors": 1}
    
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    partial_samples = find_partial_downloads(fastq_dir, quant_dir)
    
    if not partial_samples:
        return {"deleted": 0, "freed_mb": 0, "errors": 0}
    
    logger.info(f"Found {len(partial_samples)} partial downloads for {species_name}")
    
    deleted_count = 0
    freed_mb = 0
    errors = 0
    
    for sample_id, sample_dir, size_mb in partial_samples:
        try:
            if dry_run:
                logger.info(f"  [DRY RUN] Would delete {sample_id} ({size_mb}MB)")
            else:
                shutil.rmtree(sample_dir)
                logger.info(f"  🗑️  Deleted {sample_id} ({size_mb}MB)")
                deleted_count += 1
                freed_mb += size_mb
        except Exception as e:
            logger.warning(f"  ⚠️  Failed to delete {sample_dir}: {e}")
            errors += 1
    
    return {
        "deleted": deleted_count,
        "freed_mb": freed_mb,
        "errors": errors,
    }


def main():
    """Clean up partial downloads across all species."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Clean up partial and failed downloads",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without actually deleting"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually perform the cleanup (required for safety)"
    )
    
    args = parser.parse_args()
    
    if not args.dry_run and not args.execute:
        print("=" * 80)
        print("CLEANUP PARTIAL DOWNLOADS")
        print("=" * 80)
        print("\nThis script will delete partial/failed downloads that are not quantified.")
        print("\nUse --dry-run to see what would be deleted:")
        print("  python3 scripts/rna/cleanup_partial_downloads.py --dry-run")
        print("\nUse --execute to actually perform cleanup:")
        print("  python3 scripts/rna/cleanup_partial_downloads.py --execute")
        print("=" * 80)
        return 0
    
    repo_root = Path(__file__).parent.parent.parent.resolve()
    config_dir = repo_root / "config" / "amalgkit"
    
    if not config_dir.exists() or not list(config_dir.glob("amalgkit_*.yaml")):
        config_dir = repo_root / "config"
    
    # Discover all config files
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        display_name = species_code.replace("_", " ").title()
        species_configs.append((display_name, path))
    
    if not species_configs:
        logger.error("⚠️  No species configs found")
        return 1
    
    print("\n" + "=" * 80)
    if args.dry_run:
        print("DRY RUN: PARTIAL DOWNLOAD CLEANUP")
    else:
        print("CLEANING UP PARTIAL DOWNLOADS")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    print(f"Analyzing {len(species_configs)} species")
    print("=" * 80 + "\n")
    
    total_deleted = 0
    total_freed_mb = 0
    total_errors = 0
    
    for species_name, config_path in species_configs:
        if not config_path.exists():
            continue
        
        logger.info(f"Processing {species_name}...")
        result = cleanup_species(species_name, config_path, dry_run=args.dry_run)
        
        total_deleted += result["deleted"]
        total_freed_mb += result["freed_mb"]
        total_errors += result["errors"]
        
        if result["deleted"] > 0 or args.dry_run:
            logger.info(f"  {species_name}: {result['deleted']} deleted, {result['freed_mb']}MB freed")
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total samples deleted: {total_deleted}")
    print(f"Total space freed: {total_freed_mb}MB ({total_freed_mb/1024:.2f}GB)")
    print(f"Errors: {total_errors}")
    print("=" * 80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

