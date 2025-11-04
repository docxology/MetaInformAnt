#!/usr/bin/env python3
"""Fix abundance file naming for amalgkit merge compatibility.

Our workflow creates abundance.tsv but amalgkit merge expects {SRR}_abundance.tsv.
This script creates symlinks with the expected naming convention.
"""

import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
logger = logging.getLogger(__name__)


def fix_abundance_naming(species: str, base_dir: Path = Path("output/amalgkit")) -> tuple[int, int]:
    """Create symlinks from abundance.tsv to {SRR}_abundance.tsv for all samples.
    
    Args:
        species: Species identifier (e.g., 'sinvicta', 'mpharaonis')
        base_dir: Base amalgkit output directory
        
    Returns:
        Tuple of (created_count, already_exists_count)
    """
    quant_dir = base_dir / species / "quant"
    
    if not quant_dir.exists():
        logger.error(f"Quant directory not found: {quant_dir}")
        return 0, 0
    
    created = 0
    already_exists = 0
    
    # Find all SRR* directories
    srr_dirs = sorted([d for d in quant_dir.iterdir() if d.is_dir() and d.name.startswith("SRR")])
    
    logger.info(f"Processing {len(srr_dirs)} samples for {species}")
    
    for srr_dir in srr_dirs:
        srr_id = srr_dir.name
        source = srr_dir / "abundance.tsv"
        target = srr_dir / f"{srr_id}_abundance.tsv"
        
        # Check if source exists
        if not source.exists():
            continue
            
        # Check if target already exists
        if target.exists():
            if target.is_symlink():
                already_exists += 1
            continue
        
        # Create symlink
        try:
            target.symlink_to("abundance.tsv")  # Relative symlink within same directory
            created += 1
            if created <= 5:  # Log first few
                logger.debug(f"  Created: {target.name} -> abundance.tsv")
        except Exception as e:
            logger.warning(f"  Failed to create symlink for {srr_id}: {e}")
    
    logger.info(f"  ✅ Created {created} symlinks, {already_exists} already existed")
    return created, already_exists


def main():
    """Fix abundance naming for all species."""
    species_list = ["sinvicta", "mpharaonis", "cfloridanus", "pbarbatus"]
    
    logger.info("=" * 80)
    logger.info("FIXING ABUNDANCE FILE NAMING FOR AMALGKIT MERGE")
    logger.info("=" * 80)
    
    total_created = 0
    total_exists = 0
    
    for species in species_list:
        logger.info(f"\nProcessing {species}...")
        created, exists = fix_abundance_naming(species)
        total_created += created
        total_exists += exists
    
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"COMPLETE: {total_created} created, {total_exists} already existed")
    logger.info("=" * 80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())



