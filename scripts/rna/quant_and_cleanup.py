#!/usr/bin/env python3
"""
Quantify all downloaded samples and immediately cleanup FASTQs.

This script:
1. Finds all downloaded samples (have FASTQs but no quantification)
2. Runs kallisto quant on each
3. Deletes FASTQs immediately after successful quantification
4. Processes all species in parallel for speed
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_config(config_path: Path) -> dict:
    """Load amalgkit config to get species info and kallisto index."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def find_samples_needing_quantification(species_dir: Path, species_name: str) -> list[tuple[str, Path, Path]]:
    """Find samples that have FASTQs but no quantification."""
    fastq_dir = species_dir / "fastq"
    quant_dir = species_dir / "quant"
    
    if not fastq_dir.exists():
        return []
    
    needs_quant = []
    
    # Check each sample directory
    for sample_dir in fastq_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        quant_sample_dir = quant_dir / sample_id
        abundance_file = quant_sample_dir / "abundance.tsv"
        
        # Skip if already quantified
        if abundance_file.exists():
            continue
        
        # Check if FASTQs exist
        fastq_files = list(sample_dir.glob("*.fastq.gz"))
        if not fastq_files:
            continue
        
        needs_quant.append((species_name, sample_id, sample_dir))
    
    return needs_quant


def quantify_sample(species_name: str, sample_id: str, sample_dir: Path, 
                   kallisto_index: Path, output_dir: Path, threads: int = 12) -> tuple[bool, str, Path]:
    """
    Run kallisto quantification on a sample.
    
    Returns: (success, sample_id, sample_dir)
    """
    try:
        quant_output = output_dir / sample_id
        quant_output.mkdir(parents=True, exist_ok=True)
        
        # Find FASTQ files
        fastq_files = sorted(sample_dir.glob("*.fastq.gz"))
        if not fastq_files:
            logger.warning(f"  {sample_id}: No FASTQ files found")
            return (False, sample_id, sample_dir)
        
        # Build kallisto command
        cmd = [
            "kallisto", "quant",
            "-i", str(kallisto_index),
            "-o", str(quant_output),
            "-t", str(threads)
        ]
        
        # Add --single flag for single-end data
        if len(fastq_files) == 1:
            cmd.extend(["--single", "-l", "200", "-s", "20"])
        
        # Add FASTQ files
        cmd.extend([str(f) for f in fastq_files])
        
        # Run kallisto
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout per sample
        )
        
        if result.returncode == 0:
            # Verify output
            abundance_file = quant_output / "abundance.tsv"
            if abundance_file.exists():
                logger.info(f"  ✅ {species_name}/{sample_id}: Quantified successfully")
                return (True, sample_id, sample_dir)
            else:
                logger.error(f"  ❌ {species_name}/{sample_id}: No abundance.tsv created")
                return (False, sample_id, sample_dir)
        else:
            logger.error(f"  ❌ {species_name}/{sample_id}: kallisto failed: {result.stderr[:200]}")
            return (False, sample_id, sample_dir)
    
    except subprocess.TimeoutExpired:
        logger.error(f"  ❌ {species_name}/{sample_id}: Timeout after 10 minutes")
        return (False, sample_id, sample_dir)
    except Exception as e:
        logger.error(f"  ❌ {species_name}/{sample_id}: Error: {e}")
        return (False, sample_id, sample_dir)


def cleanup_fastqs(sample_dir: Path, sample_id: str, species_name: str) -> bool:
    """Delete FASTQ files after successful quantification."""
    try:
        # Get size before deletion
        size = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file())
        size_mb = size / (1024 * 1024)
        
        # Delete the entire sample directory
        shutil.rmtree(sample_dir)
        
        logger.info(f"  🗑️  {species_name}/{sample_id}: Deleted {size_mb:.1f} MB of FASTQs")
        return True
    except Exception as e:
        logger.error(f"  ❌ {species_name}/{sample_id}: Cleanup failed: {e}")
        return False


def process_species(species_name: str, config_path: Path, max_workers: int = 6) -> tuple[int, int, int]:
    """
    Process all samples for a species.
    
    Returns: (quantified, cleaned, failed)
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing {species_name}")
    logger.info(f"{'='*80}")
    
    # Load config
    config = load_config(config_path)
    
    # Get paths
    repo_root = Path(__file__).parent.parent.parent
    species_dir = repo_root / "output" / "amalgkit" / species_name
    quant_dir = species_dir / "quant"
    
    # Get kallisto index
    index_dir = species_dir / "work" / "index"
    kallisto_index = None
    if index_dir.exists():
        indices = list(index_dir.glob("*.idx"))
        if indices:
            kallisto_index = indices[0]
    
    if not kallisto_index or not kallisto_index.exists():
        logger.error(f"❌ No kallisto index found for {species_name}")
        return (0, 0, 0)
    
    # Find samples needing quantification
    samples = find_samples_needing_quantification(species_dir, species_name)
    
    if not samples:
        logger.info(f"✅ No samples need quantification for {species_name}")
        return (0, 0, 0)
    
    logger.info(f"Found {len(samples)} samples needing quantification")
    
    quantified = 0
    cleaned = 0
    failed = 0
    
    # Process samples with thread pool
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                quantify_sample,
                sp_name, sample_id, sample_dir,
                kallisto_index, quant_dir, threads=2
            ): (sp_name, sample_id, sample_dir)
            for sp_name, sample_id, sample_dir in samples
        }
        
        for future in as_completed(futures):
            sp_name, sample_id, sample_dir = futures[future]
            success, returned_id, returned_dir = future.result()
            
            if success:
                quantified += 1
                # Immediately cleanup FASTQs
                if cleanup_fastqs(returned_dir, returned_id, sp_name):
                    cleaned += 1
            else:
                failed += 1
    
    logger.info(f"\n{species_name} Summary:")
    logger.info(f"  Quantified: {quantified}")
    logger.info(f"  Cleaned up: {cleaned}")
    logger.info(f"  Failed: {failed}")
    
    return (quantified, cleaned, failed)


def main():
    """Main function to process all species."""
    logger.info("="*80)
    logger.info("QUANTIFY AND CLEANUP - All Species")
    logger.info("="*80)
    
    repo_root = Path(__file__).parent.parent.parent
    config_dir = repo_root / "config" / "amalgkit"
    
    species_configs = {
        "cfloridanus": config_dir / "amalgkit_cfloridanus.yaml",
        "pbarbatus": config_dir / "amalgkit_pbarbatus.yaml",
        "mpharaonis": config_dir / "amalgkit_mpharaonis.yaml",
        "sinvicta": config_dir / "amalgkit_sinvicta.yaml"
    }
    
    total_quantified = 0
    total_cleaned = 0
    total_failed = 0
    
    # Process each species sequentially (to avoid overwhelming kallisto)
    for species_name, config_path in species_configs.items():
        if not config_path.exists():
            logger.warning(f"⚠️  Config not found: {config_path}")
            continue
        
        quantified, cleaned, failed = process_species(
            species_name, 
            config_path,
            max_workers=6  # 6 samples in parallel per species
        )
        
        total_quantified += quantified
        total_cleaned += cleaned
        total_failed += failed
    
    # Final summary
    logger.info("\n" + "="*80)
    logger.info("FINAL SUMMARY")
    logger.info("="*80)
    logger.info(f"Total quantified: {total_quantified}")
    logger.info(f"Total cleaned up: {total_cleaned} (freed disk space)")
    logger.info(f"Total failed: {total_failed}")
    logger.info("="*80)
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

