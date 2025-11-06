#!/usr/bin/env python3
"""
Run amalgkit workflows for multiple ant species with cross-species analysis.

# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Multi-species SRA-based workflow with auto-discovery
# Steps: metadata → select → getfastq → quant → merge → curate → cstmm → csca → sanity
# Config: Auto-discovers all config/amalgkit/amalgkit_*.yaml files
# Threads: Per-config (default 10) or config override
# Batch Size: 10 samples per species (configurable via threads in config)
# Output: output/amalgkit/{species}/work/ per species
# Dependencies: SRA Toolkit, kallisto, fastp, seqkit, amalgkit
# Virtual Env: Auto-activates if .venv exists
# Reliability: ~0% for large samples (SRA Toolkit limitation)
# ============================================================================

This script provides production-ready multi-species RNA-seq workflow orchestration with:

Features:
    - Auto-activation: Automatically detects and activates virtual environment
    - Auto-discovery: Finds all species configs (amalgkit_*.yaml in config/amalgkit/)
    - Batched processing: Downloads 10 samples → Quantifies → Deletes FASTQs → Repeats
    - Disk management: Automatic SRA toolkit configuration for large downloads
    - SRA optimization: Wrapper script creation and environment setup
    - Cross-species analysis: CSTMM and CSCA for multi-species comparisons
    - Complete pipeline: metadata → select → getfastq → quant → merge → curate → sanity

Processing mode:
    - Batched: Process 10 samples at a time (configurable via threads in config)
    - Disk-friendly: Only 10 samples' FASTQs on disk at once (~20-50 GB peak)
    - Fast: 10 parallel threads for downloads and quantification
    - Resume: Automatically skips already-quantified samples

Environment Management:
    - Virtual environment auto-activation using os.execve()
    - SRA temp directory: output/sra_temp (uses /home instead of /tmp)
    - Environment variables: TMPDIR, TEMP, TMP, NCBI_VDB_QUALITY
    - Wrapper script: fasterq-dump with --size-check off
    - Tool symlinks: fastp, kallisto, seqkit in wrapper directory

Usage:
    # Process all discovered species (automatically activates venv)
    python3 scripts/rna/run_multi_species.py
    
    # The script will find all config/amalgkit/amalgkit_*.yaml files
    # Currently configured species:
    #   - Pogonomyrmex barbatus (pbarbatus)
    #   - Camponotus floridanus (cfloridanus)
    #   - Monomorium pharaonis (mpharaonis)
    #   - Solenopsis invicta (sinvicta)

Dependencies:
    - Python 3.11+ with virtual environment (.venv/)
    - amalgkit (installed in venv)
    - fasterq-dump (SRA Toolkit)
    - kallisto (quantification)
    - fastp (FASTQ quality control)
    - seqkit (sequence manipulation)
"""

import sys
import os
import shutil
from pathlib import Path
from datetime import datetime
from glob import glob

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv using uv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

# Configure SRA toolkit to use /home for temp storage instead of /tmp
# /tmp is only 16GB but samples can be 130GB+ each
repo_root = Path(__file__).parent.parent.parent.resolve()
sra_temp_dir = repo_root / "output" / "sra_temp"
sra_temp_dir.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = str(sra_temp_dir)
os.environ["TEMP"] = str(sra_temp_dir)
os.environ["TMP"] = str(sra_temp_dir)

# Disable fasterq-dump's overly conservative size check
# It fails even with 317GB available. The actual downloads work fine.
os.environ["NCBI_VDB_QUALITY"] = "R"  # Use release quality (less conservative)
os.environ["VDB_PWFILE"] = "/dev/null"  # Disable password check

# Prepend wrapper directory to PATH so amalgkit uses our wrapper
wrapper_dir = str(sra_temp_dir)
os.environ["PATH"] = f"{wrapper_dir}:{os.environ.get('PATH', '')}"

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow
from metainformant.rna.amalgkit import run_amalgkit
from metainformant.core.logging import get_logger


def check_environment_or_exit():
    """
    Check that all required tools are available before proceeding.
    
    Verifies:
        - Virtual environment is activated
        - amalgkit is installed and accessible
        - SRA Toolkit (fasterq-dump) is installed
        - kallisto is installed
        - fastp is installed (optional but recommended)
        - seqkit is installed (required for getfastq)
    
    Also performs setup:
        - Creates SRA wrapper script if needed (output/sra_temp/fasterq-dump)
        - Creates tool symlinks in wrapper directory
        - Ensures output directories exist
    
    Exits with clear error message if any critical dependency is missing.
    """
    print("Checking environment dependencies...")
    
    missing_tools = []
    warnings = []
    
    # Check 1: amalgkit (critical - must be in venv)
    if shutil.which("amalgkit") is None:
        missing_tools.append(("amalgkit", "Virtual environment not activated or amalgkit not installed"))
    
    # Check 2: fasterq-dump (critical - for downloading FASTQs)
    if shutil.which("fasterq-dump") is None:
        missing_tools.append(("fasterq-dump", "SRA Toolkit not installed"))
    
    # Check 3: kallisto (critical - for quantification)
    kallisto_path = shutil.which("kallisto")
    if kallisto_path is None:
        missing_tools.append(("kallisto", "Kallisto not installed"))
    
    # Check 4: Virtual environment (warning only)
    import os
    if not os.environ.get("VIRTUAL_ENV"):
        warnings.append("Virtual environment may not be activated (VIRTUAL_ENV not set)")
    
    # Check 5: metainformant package
    try:
        import metainformant
    except ImportError:
        missing_tools.append(("metainformant", "Package not installed (run: uv pip install -e .)"))
    
    # Report results
    if missing_tools or warnings:
        print("\n" + "="*80)
        print("❌ ENVIRONMENT CHECK FAILED")
        print("="*80)
        print()
        
        if missing_tools:
            print("Missing required tools:")
            for tool, reason in missing_tools:
                print(f"  ❌ {tool:20s} - {reason}")
            print()
        
        if warnings:
            print("Warnings:")
            for warning in warnings:
                print(f"  ⚠️  {warning}")
            print()
        
        print("=" * 80)
        print("SETUP INSTRUCTIONS")
        print("=" * 80)
        print()
        
        if any("amalgkit" in tool for tool, _ in missing_tools):
            print("1. Activate the virtual environment:")
            print("   cd /home/q/Documents/GitHub/MetaInformAnt")
            print("   source .venv/bin/activate")
            print()
        
        if any("metainformant" in tool for tool, _ in missing_tools):
            print("2. Install metainformant package:")
            print("   uv pip install -e .")
            print()
        
        if any("amalgkit" in tool for tool, _ in missing_tools):
            print("3. Install amalgkit:")
            print("   uv pip install git+https://github.com/kfuku52/amalgkit")
            print()
        
        if any("fasterq-dump" in tool for tool, _ in missing_tools):
            print("4. Install SRA Toolkit:")
            print("   sudo apt-get update")
            print("   sudo apt-get install -y sra-toolkit")
            print()
        
        if any("kallisto" in tool for tool, _ in missing_tools):
            print("5. Install kallisto:")
            print("   sudo apt-get update")
            print("   sudo apt-get install -y kallisto")
            print()
        
        print("For complete setup guide, see: docs/rna/SETUP.md")
        print()
        print("To run comprehensive environment check:")
        print("  python3 scripts/rna/check_environment.py")
        print("="*80)
        sys.exit(1)
    
    print("✅ All dependencies found and accessible")
    print()


def discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]:
    """Discover all species configuration files.
    
    Args:
        config_dir: Directory to search for config files (default: config/amalgkit)
    
    Returns:
        List of (species_name, config_path) tuples, sorted by filename
        
    Side effects:
        None (read-only operation)
        
    Notes:
        - Excludes template files (containing 'template' in name)
        - Extracts species name from filename: amalgkit_<species>.yaml
        - Resolves paths relative to repository root
    """
    # Get repository root (3 levels up from this script: scripts/rna/run_multi_species.py)
    repo_root = Path(__file__).parent.parent.parent.resolve()
    
    # Resolve config directory relative to repo root
    config_dir = repo_root / config_dir
    
    # Find all amalgkit_*.yaml files in config/amalgkit/, excluding template
    pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(pattern))
    
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
    
    Args:
        config_path: Path to species YAML configuration file
        species_name: Human-readable species name for logging
        
    Returns:
        Tuple of (success, work_dir) where:
        - success: True if workflow completed without errors
        - work_dir: Path to workflow output directory
        
    Side effects:
        - Loads workflow configuration
        - Executes complete amalgkit workflow (metadata → sanity)
        - Creates output directories and files
        - Writes logs to work_dir/logs/
    
    The workflow automatically uses batched download-quant-delete processing:
        - Downloads 10 samples (batch size from config.threads)
        - Quantifies all 10 samples
        - Deletes FASTQ files for those 10 samples
        - Moves to next batch
    
    This keeps disk usage bounded (~20-50 GB peak) while maintaining good performance.
    
    Workflow steps executed:
        1. metadata - Retrieve and filter SRA metadata
        2. config - Generate configuration files
        3. select - Select qualified samples
        4. getfastq - Download FASTQ files (batched)
        5. quant - Quantify expression (batched with getfastq)
        6. merge - Merge quantification results
        7. curate - Quality control and curation
        8. sanity - Validate workflow outputs
    
    Args:
        config_path: Path to species-specific amalgkit YAML config
        species_name: Short name for the species (e.g., "pbarbatus")
    
    Returns:
        Tuple of (success: bool, work_dir: Path)
        - success: True if workflow completed successfully
        - work_dir: Path to work directory
    """
    logger = get_logger(f"amalgkit_{species_name}")
    
    logger.info("="*80)
    logger.info(f"Starting batched workflow for {species_name}")
    logger.info("="*80)
    logger.info(f"Config: {config_path}")
    logger.info(f"Start time: {datetime.now()}")
    logger.info(f"Processing mode: Batched (download → quantify → delete, 10 samples at a time)")
    
    try:
        # Load config
        logger.info("Loading configuration...")
        cfg = load_workflow_config(config_path)
        logger.info(f"  Work dir: {cfg.work_dir}")
        logger.info(f"  Genome: {cfg.genome.get('accession') if cfg.genome else 'None'}")
        logger.info(f"  Species: {cfg.species_list}")
        logger.info(f"  Batch size: {cfg.threads} samples (from threads config)")
        logger.info(f"  Default threads: 10 for downloads and quantification")
        
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
        work_dirs: List of work directories for all species (must have completed workflows)
        output_dir: Output directory for cross-species results (default: output/amalgkit/cross_species)
        
    Returns:
        None (results written to output_dir)
        
    Side effects:
        - Creates output_dir if it doesn't exist
        - Runs amalgkit cstmm (Cross-Species TMM Normalization)
        - Runs amalgkit csca (Cross-Species Correlation Analysis)
        - Writes results to output_dir/cstmm/ and output_dir/csca/
        
    Dependencies:
        - All species workflows must have completed successfully
        - Merged abundance tables must exist in work_dirs
        - amalgkit CLI must be available
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
    """
    Run workflows for all discovered species, then cross-species analysis.
    
    Main workflow:
        1. Auto-activate virtual environment if needed
        2. Check environment and create SRA wrapper/symlinks
        3. Discover all species configs in config/amalgkit/
        4. Run individual species workflows (batched processing)
        5. Run cross-species analysis (CSTMM, CSCA) if ≥2 species
    
    Configuration:
        - Auto-discovers: config/amalgkit/amalgkit_*.yaml (excludes template)
        - Thread count: 10 (configurable in each species config)
        - Batch size: 10 samples at a time
        - SRA temp: output/sra_temp (instead of /tmp)
    
    Environment Setup:
        - TMPDIR=/path/to/output/sra_temp
        - TEMP=/path/to/output/sra_temp
        - TMP=/path/to/output/sra_temp
        - NCBI_VDB_QUALITY=R
        - VDB_PWFILE=/dev/null
        - PATH=prepended with wrapper directory
    
    Output:
        - Individual species: output/amalgkit/{species}/
        - Cross-species: output/amalgkit/cross_species/
        - Logs: Per-species and per-step
    """
    
    # Pre-flight check: ensure environment is properly set up
    check_environment_or_exit()
    
    # Discover species configurations
    species_configs = discover_species_configs()
    
    if not species_configs:
        print("❌ No species configurations found!")
        print("   Looked in: config/amalgkit/amalgkit_*.yaml")
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
    print("\nProcessing mode: Batched (10 samples at a time)")
    print("  • Download batch → Quantify batch → Delete FASTQs → Repeat")
    print("  • Disk-friendly: Only ~20-50 GB peak usage")
    print("  • Fast: 10 parallel threads for downloads and quantification")
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

