#!/usr/bin/env python3
"""
Run amalgkit workflows for multiple species with sequential sample processing.

This script implements disk-space-friendly processing:
1. Download one SRA at a time
2. Quantify it immediately  
3. Delete FASTQ files to free space
4. Move to next sample

This prevents disk exhaustion when processing hundreds of samples.
"""

import sys
from pathlib import Path
from datetime import datetime
from glob import glob

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow
from metainformant.rna.steps import process_samples_sequentially
from metainformant.rna.amalgkit import run_amalgkit
from metainformant.core.logging import get_logger


def discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]:
    """Discover all species configuration files.
    
    Looks in config/amalgkit/ directory for amalgkit_*.yaml files.
    """
    # Find all amalgkit_*.yaml files in config/amalgkit/, excluding template
    pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(pattern))
    
    # Also check old location for backwards compatibility
    if not config_files:
        old_pattern = str(Path("config") / "amalgkit_*.yaml")
        config_files = sorted(glob(old_pattern))
        if config_files:
            logger = get_logger("amalgkit_discovery")
            logger.warning(f"Found configs in old location (config/), consider moving to config/amalgkit/")
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        if "template" in path.stem.lower():
            continue
        
        species_name = path.stem.replace("amalgkit_", "").replace("_", " ").title()
        species_configs.append((species_name, path))
    
    return species_configs


def run_species_with_sequential_processing(config_path: Path, species_name: str) -> tuple[bool, Path]:
    """Run workflow for a species using sequential sample processing."""
    logger = get_logger(f"amalgkit_{species_name}")
    
    logger.info("="*80)
    logger.info(f"Starting sequential workflow for {species_name}")
    logger.info("="*80)
    logger.info(f"Config: {config_path}")
    logger.info(f"Start time: {datetime.now()}")
    
    try:
        # Load config
        logger.info("Loading configuration...")
        cfg = load_workflow_config(config_path)
        logger.info(f"  Work dir: {cfg.work_dir}")
        logger.info(f"  Genome: {cfg.genome.get('accession') if cfg.genome else 'None'}")
        logger.info(f"  Species: {cfg.species_list}")
        
        # Phase 1: Run initial steps (metadata, config, select)
        logger.info("\n=== Phase 1: Initial Setup ===")
        initial_steps = ["metadata", "config", "select"]
        
        for step in initial_steps:
            logger.info(f"Running step: {step}")
            step_params = cfg.per_step.get(step, {})
            
            # Run amalgkit step
            result = run_amalgkit(
                step,
                step_params,
                work_dir=cfg.work_dir,
                log_dir=cfg.log_dir,
                step_name=step,
                check=False,
            )
            
            if result.returncode != 0 and result.returncode != 204:
                logger.warning(f"Step {step} returned code {result.returncode}")
        
        # Phase 2: Sequential sample processing (getfastq + quant)
        logger.info("\n=== Phase 2: Sequential Sample Processing ===")
        logger.info("Processing samples one-at-a-time: download → quantify → cleanup")
        
        metadata_path = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        if not metadata_path.exists():
            metadata_path = cfg.work_dir / "metadata" / "metadata.tsv"
        
        if not metadata_path.exists():
            logger.error(f"Metadata file not found: {metadata_path}")
            return False, cfg.work_dir
        
        # Get parameters for getfastq and quant
        getfastq_params = cfg.per_step.get("getfastq", {}).copy()
        quant_params = cfg.per_step.get("quant", {}).copy()
        
        # Ensure we have necessary paths
        fastq_dir = Path(getfastq_params.get("out_dir", cfg.work_dir.parent / "fastq"))
        quant_dir = Path(quant_params.get("out_dir", cfg.work_dir / "quant"))
        
        # Run sequential processing
        stats = process_samples_sequentially(
            metadata_path=metadata_path,
            fastq_dir=fastq_dir,
            quant_dir=quant_dir,
            work_dir=cfg.work_dir,
            log_dir=cfg.log_dir,
            getfastq_params=getfastq_params,
            quant_params=quant_params,
        )
        
        logger.info("\n=== Sequential Processing Results ===")
        logger.info(f"  Total samples: {stats['total_samples']}")
        logger.info(f"  Already quantified: {stats['already_quantified']}")
        logger.info(f"  Newly processed: {stats['newly_processed']}")
        logger.info(f"  Failed: {stats['failed']}")
        
        # Phase 3: Run final steps (integrate, merge, curate, sanity)
        logger.info("\n=== Phase 3: Final Processing ===")
        final_steps = ["integrate", "merge", "curate", "csca", "sanity"]
        
        for step in final_steps:
            if step not in cfg.per_step:
                continue
            
            logger.info(f"Running step: {step}")
            step_params = cfg.per_step.get(step, {})
            
            result = run_amalgkit(
                step,
                step_params,
                work_dir=cfg.work_dir,
                log_dir=cfg.log_dir,
                step_name=step,
                check=False,
            )
            
            if result.returncode != 0 and result.returncode != 204:
                logger.warning(f"Step {step} returned code {result.returncode}")
        
        logger.info("\n" + "="*80)
        logger.info(f"✅ {species_name} workflow completed!")
        logger.info(f"End time: {datetime.now()}")
        logger.info("="*80)
        
        return True, cfg.work_dir
        
    except Exception as e:
        logger.error(f"\n❌ Error running workflow for {species_name}: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def main():
    """Run workflows for all discovered species with sequential processing."""
    
    # Discover species configurations
    species_configs = discover_species_configs()
    
    if not species_configs:
        print("❌ No species configurations found in config/ directory!")
        sys.exit(1)
    
    print("\n" + "="*80)
    print("MULTI-SPECIES AMALGKIT WORKFLOW (SEQUENTIAL PROCESSING)")
    print("="*80)
    print(f"Date: {datetime.now()}")
    print(f"Species discovered: {len(species_configs)}")
    for name, config in species_configs:
        print(f"  - {name} ({config})")
    print("="*80 + "\n")
    
    print("Processing mode: Download → Quantify → Delete (one sample at a time)")
    print("This prevents disk space exhaustion with large datasets.\n")
    
    # Process each species
    results = {}
    work_dirs = []
    
    for species_name, config_path in species_configs:
        success, work_dir = run_species_with_sequential_processing(config_path, species_name)
        results[species_name] = success
        if success and work_dir:
            work_dirs.append(work_dir)
    
    # Final summary
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    
    print("\nIndividual Species:")
    for species, success in results.items():
        status = "✅ SUCCESS" if success else "❌ FAILED"
        print(f"  {species:30s}: {status}")
    
    print("="*80)
    
    # Exit code
    all_success = all(results.values())
    if all_success:
        print("\n✅ All species processed successfully!")
        sys.exit(0)
    else:
        print("\n⚠️  Some species had failures - check logs")
        sys.exit(1)


if __name__ == "__main__":
    main()

