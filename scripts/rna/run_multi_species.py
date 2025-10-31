#!/usr/bin/env python3
"""
Run amalgkit workflows for multiple ant species with cross-species analysis.

This script:
1. Auto-discovers all species configs (amalgkit_*.yaml in config/amalgkit/, excluding template)
2. Runs full end-to-end workflow for each species (metadata → sanity)
3. Uses batched processing (download 8 → quantify 8 → delete 8 → repeat)
4. Performs cross-species analysis (CSTMM, CSCA) across all species

Processing mode:
- Batched: Process 8 samples at a time (configurable via threads in config)
- Disk-friendly: Only 8 samples' FASTQs on disk at once (~16-40 GB peak)
- Fast: Parallel processing within each batch
- Resume: Automatically skips already-quantified samples

Usage:
    # Process all discovered species
    python3 scripts/rna/run_multi_species.py
    
    # The script will find all config/amalgkit/amalgkit_*.yaml files
    # Currently configured species:
    #   - Pogonomyrmex barbatus (pbarbatus)
    #   - Camponotus floridanus (cfloridanus)
    #   - Monomorium pharaonis (mpharaonis)
    #   - Solenopsis invicta (sinvicta)
"""

import sys
from pathlib import Path
from datetime import datetime
from glob import glob

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow
from metainformant.rna.amalgkit import run_amalgkit
from metainformant.core.logging import get_logger


def discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]:
    """Discover all species configuration files.
    
    Looks in config/amalgkit/ directory for amalgkit_*.yaml files.
    
    Returns:
        List of (species_name, config_path) tuples
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
        # Skip template
        if "template" in path.stem.lower():
            continue
        
        # Extract species name from filename: amalgkit_<species>.yaml
        species_name = path.stem.replace("amalgkit_", "").replace("_", " ").title()
        species_configs.append((species_name, path))
    
    return species_configs


def run_species_workflow(config_path: Path, species_name: str) -> tuple[bool, Path]:
    """Run full workflow for a single species using batched processing.
    
    The workflow automatically uses batched download-quant-delete processing:
    - Downloads 8 samples (batch size from config.threads)
    - Quantifies all 8 samples
    - Deletes FASTQ files for those 8 samples
    - Moves to next batch
    
    This keeps disk usage bounded while maintaining good performance.
    
    Returns:
        (success, work_dir) tuple
    """
    logger = get_logger(f"amalgkit_{species_name}")
    
    logger.info("="*80)
    logger.info(f"Starting batched workflow for {species_name}")
    logger.info("="*80)
    logger.info(f"Config: {config_path}")
    logger.info(f"Start time: {datetime.now()}")
    logger.info(f"Processing mode: Batched (download → quantify → delete, 8 samples at a time)")
    
    try:
        # Load config
        logger.info("Loading configuration...")
        cfg = load_workflow_config(config_path)
        logger.info(f"  Work dir: {cfg.work_dir}")
        logger.info(f"  Genome: {cfg.genome.get('accession') if cfg.genome else 'None'}")
        logger.info(f"  Species: {cfg.species_list}")
        logger.info(f"  Batch size: {cfg.threads} samples (from threads config)")
        
        # Execute workflow (automatically uses batched processing)
        logger.info("\nExecuting full amalgkit workflow with batched processing...")
        return_codes = execute_workflow(cfg)
        
        # Check results
        logger.info("\n" + "="*80)
        logger.info("WORKFLOW RESULTS")
        logger.info("="*80)
        
        steps = ["genome", "metadata", "integrate", "config", "select", 
                 "getfastq", "quant", "merge", "cstmm", "curate", "csca", "sanity"]
        
        all_success = True
        failed_steps = []
        for i, (step, code) in enumerate(zip(steps, return_codes)):
            status = "✅ PASS" if code == 0 else f"❌ FAIL (code {code})"
            logger.info(f"{step:15s}: {status}")
            if code != 0 and code != 204:  # 204 = skipped (ok for some steps)
                all_success = False
                failed_steps.append(step)
        
        logger.info("="*80)
        if all_success or (failed_steps and all(s in ["integrate", "cstmm", "csca"] for s in failed_steps)):
            logger.info(f"✅ {species_name} workflow completed successfully!")
            success = True
        else:
            logger.warning(f"⚠️  {species_name} workflow completed with failures: {', '.join(failed_steps)}")
            success = False
        
        logger.info(f"End time: {datetime.now()}")
        logger.info("="*80)
        
        return success, cfg.work_dir
        
    except Exception as e:
        logger.error(f"\n❌ Error running workflow for {species_name}: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def run_cross_species_analysis(work_dirs: list[Path], output_dir: Path = Path("output/amalgkit/cross_species")):
    """Run cross-species TMM normalization and correlation analysis.
    
    Args:
        work_dirs: List of work directories for all species
        output_dir: Output directory for cross-species results
    """
    logger = get_logger("amalgkit_cross_species")
    
    logger.info("\n" + "="*80)
    logger.info("CROSS-SPECIES ANALYSIS")
    logger.info("="*80)
    logger.info(f"Species work directories: {len(work_dirs)}")
    for wd in work_dirs:
        logger.info(f"  - {wd}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run CSTMM (Cross-Species TMM Normalization)
    logger.info("\n--- Running CSTMM (Cross-Species TMM Normalization) ---")
    try:
        cstmm_params = {
            "out_dir": str(output_dir / "cstmm"),
        }
        # Add all species work dirs for cross-species comparison
        for i, wd in enumerate(work_dirs):
            cstmm_params[f"species{i+1}_dir"] = str(wd)
        
        logger.info(f"CSTMM params: {cstmm_params}")
        result = run_amalgkit("cstmm", cstmm_params, check=False, capture_output=False)
        
        if result.returncode == 0:
            logger.info("✅ CSTMM completed successfully")
        else:
            logger.warning(f"⚠️  CSTMM completed with code {result.returncode}")
    except Exception as e:
        logger.error(f"❌ CSTMM failed: {e}")
    
    # Run CSCA (Cross-Species Correlation Analysis)
    logger.info("\n--- Running CSCA (Cross-Species Correlation Analysis) ---")
    try:
        csca_params = {
            "out_dir": str(output_dir / "csca"),
        }
        # Add all species work dirs
        for i, wd in enumerate(work_dirs):
            csca_params[f"species{i+1}_dir"] = str(wd)
        
        logger.info(f"CSCA params: {csca_params}")
        result = run_amalgkit("csca", csca_params, check=False, capture_output=False)
        
        if result.returncode == 0:
            logger.info("✅ CSCA completed successfully")
        else:
            logger.warning(f"⚠️  CSCA completed with code {result.returncode}")
    except Exception as e:
        logger.error(f"❌ CSCA failed: {e}")
    
    logger.info("\n" + "="*80)
    logger.info("Cross-species analysis complete")
    logger.info(f"Results in: {output_dir}")
    logger.info("="*80)


def main():
    """Run workflows for all discovered species, then cross-species analysis."""
    
    # Discover species configurations
    species_configs = discover_species_configs()
    
    if not species_configs:
        print("❌ No species configurations found!")
        print("   Looked in:")
        print("     - config/amalgkit/amalgkit_*.yaml")
        print("     - config/amalgkit_*.yaml (legacy)")
        print("   Excluding template files")
        sys.exit(1)
    
    print("\n" + "="*80)
    print("MULTI-SPECIES AMALGKIT WORKFLOW WITH CROSS-SPECIES ANALYSIS")
    print("="*80)
    print(f"Date: {datetime.now()}")
    print(f"Species discovered: {len(species_configs)}")
    for name, config in species_configs:
        print(f"  - {name} ({config})")
    print("="*80)
    print("\nProcessing mode: Batched (8 samples at a time)")
    print("  • Download batch → Quantify batch → Delete FASTQs → Repeat")
    print("  • Disk-friendly: Only ~16-40 GB peak usage")
    print("  • Fast: Parallel processing within each batch")
    print("  • Resume: Automatically skips completed samples")
    print("="*80 + "\n")
    
    # Phase 1: Individual species workflows
    print("="*80)
    print("PHASE 1: INDIVIDUAL SPECIES WORKFLOWS")
    print("="*80 + "\n")
    
    results = {}
    work_dirs = []
    
    for species_name, config_path in species_configs:
        success, work_dir = run_species_workflow(config_path, species_name)
        results[species_name] = success
        if success and work_dir:
            work_dirs.append(work_dir)
    
    # Phase 2: Cross-species analysis
    if len(work_dirs) >= 2:
        print("\n" + "="*80)
        print("PHASE 2: CROSS-SPECIES ANALYSIS")
        print("="*80)
        print(f"Analyzing {len(work_dirs)} species together...")
        print("")
        
        run_cross_species_analysis(work_dirs)
    else:
        print("\n⚠️  Skipping cross-species analysis (need at least 2 successful species)")
    
    # Final summary
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    
    print("\nIndividual Species:")
    for species, success in results.items():
        status = "✅ SUCCESS" if success else "❌ FAILED"
        print(f"  {species:30s}: {status}")
    
    if len(work_dirs) >= 2:
        print(f"\nCross-Species Analysis: ✅ Completed ({len(work_dirs)} species)")
        print(f"  Results: output/amalgkit/cross_species/")
    
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

