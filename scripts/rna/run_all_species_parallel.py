#!/usr/bin/env python3
"""
Run all 24 species workflows in parallel with 24 threads total (1 thread per species).

This script runs the complete end-to-end workflow (metadata → sanity) for all species
simultaneously, with 1 thread allocated per species.

# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: All 24 species in parallel with 1 thread per species
# Steps: metadata → select → getfastq → quant → merge → curate → cstmm → csca → sanity
# Config: Auto-discovers all config/amalgkit/amalgkit_*.yaml files
# Threads: 1 thread per species (24 total)
# Batch Size: Per-species (configurable in each species config)
# Output: output/amalgkit/{species}/work/ per species
# Dependencies: SRA Toolkit, kallisto, fastp, seqkit, amalgkit
# Virtual Env: Auto-activates if .venv exists
# Parallelism: All species run simultaneously (multiprocessing)
# ============================================================================

Usage:
    # Run all 24 species in parallel with 1 thread each (24 total threads)
    export AK_THREADS=1
    python3 scripts/rna/run_all_species_parallel.py
    
    # Or specify threads explicitly
    python3 scripts/rna/run_all_species_parallel.py --threads-per-species 1
"""

import sys
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow
from metainformant.core.logging import get_logger

logger = get_logger("parallel_species")


def discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]:
    """Discover all species configuration files."""
    configs = sorted(config_dir.glob("amalgkit_*.yaml"))
    # Exclude template and test configs
    configs = [
        c for c in configs
        if "template" not in c.stem.lower() and "test" not in c.stem.lower()
    ]
    
    species_configs = []
    for config_path in configs:
        # Extract species name from filename
        species_name = config_path.stem.replace("amalgkit_", "")
        species_configs.append((species_name, config_path))
    
    return species_configs


def run_species_workflow_worker(species_name: str, config_path: str, threads_per_species: int) -> tuple[bool, str]:
    """Worker function to run workflow for a single species (runs in separate process).
    
    Args:
        species_name: Human-readable species name
        config_path: Path to species YAML configuration file (as string for multiprocessing)
        threads_per_species: Number of threads to use for this species
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    # Import in worker process (must be done here for multiprocessing)
    import sys
    from pathlib import Path
    import os
    
    # Add paths
    repo_root = Path(__file__).parent.parent.parent.parent.resolve()
    sys.path.insert(0, str(repo_root / "src"))
    
    from metainformant.core.logging import get_logger
    from metainformant.rna.workflow import load_workflow_config, execute_workflow
    
    worker_logger = get_logger(f"parallel_{species_name}")
    config_path = Path(config_path)
    
    try:
        # Set thread count via environment variable for this process
        os.environ["AK_THREADS"] = str(threads_per_species)
        
        worker_logger.info("="*80)
        worker_logger.info(f"Starting workflow for {species_name}")
        worker_logger.info(f"Config: {config_path}")
        worker_logger.info(f"Threads: {threads_per_species}")
        worker_logger.info(f"Start time: {datetime.now()}")
        worker_logger.info("="*80)
        
        # Load config
        cfg = load_workflow_config(config_path)
        worker_logger.info(f"Work dir: {cfg.work_dir}")
        worker_logger.info(f"Species: {cfg.species_list}")
        
        # Execute workflow
        worker_logger.info("Executing full amalgkit workflow...")
        return_codes = execute_workflow(cfg, check=False)
        
        # Check results
        success = all(code == 0 or code == 126 for code in return_codes)  # 126 = skipped due to deps
        failed_steps = [i for i, code in enumerate(return_codes) if code != 0 and code != 126]
        
        worker_logger.info("="*80)
        if success:
            worker_logger.info(f"✅ {species_name}: SUCCESS")
            message = f"Completed successfully"
        else:
            worker_logger.warning(f"⚠️  {species_name}: Some steps failed")
            message = f"Failed steps: {failed_steps}"
        worker_logger.info("="*80)
        
        return success, message
        
    except Exception as e:
        worker_logger.error(f"❌ {species_name}: Exception: {e}")
        return False, f"Exception: {e}"


def main():
    """Main entry point for parallel species workflow execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Run all species workflows in parallel with 1 thread per species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--threads-per-species",
        type=int,
        default=None,
        help="Threads per species (default: 1, or from AK_THREADS env var)"
    )
    parser.add_argument(
        "--max-species",
        type=int,
        default=None,
        help="Maximum number of species to process (default: all)"
    )
    
    args = parser.parse_args()
    
    # Pre-flight check (auto-setup enabled)
    check_environment_or_exit(auto_setup=True)
    
    # Discover species
    species_configs = discover_species_configs()
    
    if not species_configs:
        logger.error("No species configurations found!")
        logger.error("Looked in: config/amalgkit/amalgkit_*.yaml")
        sys.exit(1)
    
    # Limit species if requested
    if args.max_species:
        species_configs = species_configs[:args.max_species]
    
    # Determine threads per species
    threads_per_species = args.threads_per_species or int(os.environ.get("AK_THREADS", "1"))
    
    print("\n" + "="*80)
    print("PARALLEL MULTI-SPECIES AMALGKIT WORKFLOW")
    print("="*80)
    print(f"Date: {datetime.now()}")
    print(f"Species: {len(species_configs)}")
    print(f"Threads per species: {threads_per_species}")
    print(f"Total threads: {len(species_configs) * threads_per_species}")
    print(f"Execution mode: Parallel (all species simultaneously)")
    print("="*80)
    print()
    
    for name, config in species_configs:
        print(f"  - {name} ({config})")
    print()
    
    # Run all species in parallel
    results = {}
    
    with ProcessPoolExecutor(max_workers=len(species_configs)) as executor:
        # Submit all species workflows (convert Path to str for multiprocessing)
        future_to_species = {
            executor.submit(run_species_workflow_worker, name, str(config), threads_per_species): name
            for name, config in species_configs
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_species):
            species_name = future_to_species[future]
            try:
                success, message = future.result()
                results[species_name] = {"success": success, "message": message}
                status = "✅" if success else "❌"
                print(f"{status} {species_name}: {message}")
            except Exception as e:
                results[species_name] = {"success": False, "message": f"Exception: {e}"}
                print(f"❌ {species_name}: Exception: {e}")
    
    # Final summary
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    
    successful = sum(1 for r in results.values() if r["success"])
    failed = len(results) - successful
    
    print(f"\nResults: {successful} successful, {failed} failed")
    print()
    
    for species_name, result in sorted(results.items()):
        status = "✅ SUCCESS" if result["success"] else "❌ FAILED"
        print(f"  {species_name:30s}: {status} - {result['message']}")
    
    print("="*80)
    
    # Exit code
    if all(r["success"] for r in results.values()):
        print("\n✅ All species processed successfully!")
        sys.exit(0)
    else:
        print("\n⚠️  Some species had failures - check logs")
        sys.exit(1)


if __name__ == "__main__":
    main()

