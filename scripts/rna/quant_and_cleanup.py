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
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.steps import quantify_sample, delete_sample_fastqs
from metainformant.rna import find_unquantified_samples
from metainformant.core.logging import get_logger
from metainformant.core.io import read_delimited

logger = get_logger("quant_cleanup")


# Functions now use metainformant - removed duplicate implementations


def process_species(species_name: str, config_path: Path, max_workers: int = 6) -> tuple[int, int, int]:
    """
    Process all samples for a species using metainformant functions.
    
    Returns: (quantified, cleaned, failed)
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing {species_name}")
    logger.info(f"{'='*80}")
    
    # Load config
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Find unquantified samples using metainformant function
    unquantified = find_unquantified_samples(config_path)
    
    if not unquantified:
        logger.info(f"✅ No samples need quantification for {species_name}")
        return (0, 0, 0)
    
    logger.info(f"Found {len(unquantified)} samples needing quantification")
    
    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        logger.error(f"❌ No metadata file found for {species_name}")
        return (0, 0, len(unquantified))
    
    rows = list(read_delimited(metadata_file, delimiter="\t"))
    
    # Get quant params from config
    quant_params = dict(cfg.per_step.get("quant", {}))
    quant_params["out_dir"] = str(quant_dir.absolute())
    quant_params["threads"] = cfg.threads or 12
    
    # Inject index_dir if needed
    if "index_dir" not in quant_params and "index-dir" not in quant_params:
        index_dir = quant_dir.parent / "work" / "index"
        if index_dir.exists():
            quant_params["index_dir"] = str(index_dir.absolute())
    
    quantified = 0
    cleaned = 0
    failed = 0
    
    # Process samples with thread pool
    def process_sample(sample_id: str) -> tuple[bool, str]:
        """Process a single sample."""
        try:
            # Get metadata rows for this sample
            sample_rows = [row for row in rows if row.get("run") == sample_id]
            
            if not sample_rows:
                logger.warning(f"  ⚠️  {sample_id} not in metadata, skipping")
                return (False, "Not in metadata")
            
            # Use metainformant quantify_sample function
            success, message, abundance_file = quantify_sample(
                sample_id=sample_id,
                metadata_rows=sample_rows,
                quant_params=quant_params,
                log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
                step_name=f"quant_{species_name}_{sample_id}",
            )
            
            if success:
                # Delete FASTQs immediately
                delete_sample_fastqs(sample_id, fastq_dir)
                return (True, message)
            else:
                # Still delete FASTQs to free space
                delete_sample_fastqs(sample_id, fastq_dir)
                return (False, message)
        
        except Exception as e:
            logger.error(f"  ❌ Error processing {sample_id}: {e}")
            # Attempt cleanup
            try:
                delete_sample_fastqs(sample_id, fastq_dir)
            except Exception:
                pass
            return (False, str(e))
    
    # Process with thread pool
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_sample, sample_id): sample_id
            for sample_id in unquantified
        }
        
        for future in as_completed(futures):
            sample_id = futures[future]
            success, message = future.result()
            
            if success:
                quantified += 1
                cleaned += 1
                logger.info(f"  ✅ {species_name}/{sample_id}: {message}")
            else:
                failed += 1
                logger.warning(f"  ⚠️  {species_name}/{sample_id}: {message}")
    
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
    
    # Discover all config files
    if not config_dir.exists():
        config_dir = repo_root / "config"
    
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = {}
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        species_configs[species_code] = path
    
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

