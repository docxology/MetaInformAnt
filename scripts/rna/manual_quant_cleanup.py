#!/usr/bin/env python3
"""
Manually quantify completed downloads and cleanup FASTQs.
Runs sequentially with proper threading to avoid overwhelming the system.
"""

import sys
from pathlib import Path
from glob import glob

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.steps import quantify_sample, delete_sample_fastqs
from metainformant.rna import find_unquantified_samples
from metainformant.core.logging import get_logger
from metainformant.core.io import read_delimited

logger = get_logger("manual_quant_cleanup")


# Functions now use metainformant - removed duplicate implementations


def process_species(species_name: str, config_path: Path):
    """Process one species using metainformant functions."""
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
        logger.info(f"✅ No complete samples need quantification")
        return
    
    logger.info(f"Found {len(unquantified)} samples ready for quantification")
    
    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        logger.error(f"❌ No metadata file found for {species_name}")
        return
    
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
    
    # Process sequentially to avoid overwhelming system
    for sample_id in unquantified:
        try:
            # Get metadata rows for this sample
            sample_rows = [row for row in rows if row.get("run") == sample_id]
            
            if not sample_rows:
                logger.warning(f"  ⚠️  {sample_id} not in metadata, skipping")
                continue
            
            # Use metainformant quantify_sample function
            logger.info(f"  Quantifying {species_name}/{sample_id} ...")
            success, message, abundance_file = quantify_sample(
                sample_id=sample_id,
                metadata_rows=sample_rows,
                quant_params=quant_params,
                log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
                step_name=f"quant_{species_name}_{sample_id}",
            )
            
            if success:
                logger.info(f"  ✅ {species_name}/{sample_id}: {message}")
                quantified += 1
            else:
                logger.error(f"  ❌ {species_name}/{sample_id}: {message}")
            
            # Always cleanup FASTQs to free space
            delete_sample_fastqs(sample_id, fastq_dir)
            cleaned += 1
        
        except Exception as e:
            logger.error(f"  ❌ {species_name}/{sample_id}: {e}")
            # Attempt cleanup
            try:
                delete_sample_fastqs(sample_id, fastq_dir)
                cleaned += 1
            except Exception:
                pass
    
    logger.info(f"\n{species_name} Summary: {quantified} quantified, {cleaned} cleaned")


def main():
    logger.info("="*80)
    logger.info("MANUAL QUANTIFICATION AND CLEANUP")
    logger.info("="*80)
    
    repo_root = Path(__file__).parent.parent.parent
    config_dir = repo_root / "config" / "amalgkit"
    
    # Discover all config files
    if not config_dir.exists():
        config_dir = repo_root / "config"
    
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_code = path.stem.replace("amalgkit_", "")
        species_configs.append((species_code, path))
    
    for species_name, config_path in species_configs:
        if config_path.exists():
            process_species(species_name, config_path)
    
    logger.info("\n" + "="*80)
    logger.info("COMPLETE")
    logger.info("="*80)


if __name__ == "__main__":
    main()

