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
import time
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
from metainformant.rna.progress_tracker import get_tracker
from metainformant.rna.monitoring import analyze_species_status
from metainformant.core.io import read_delimited

logger = get_logger("parallel_species")


def discover_species_configs(config_dir: Path = Path("config/amalgkit")) -> list[tuple[str, Path]]:
    """Discover all species configuration files.
    
    Args:
        config_dir: Directory containing species config files
        
    Returns:
        List of (species_name, config_path) tuples
    """
    logger.info(f"Discovering species configs in {config_dir}...")
    configs = sorted(config_dir.glob("amalgkit_*.yaml"))
    
    # Exclude template and test configs
    skipped = []
    valid_configs = []
    for config_path in configs:
        if "template" in config_path.stem.lower() or "test" in config_path.stem.lower():
            skipped.append(config_path.name)
            continue
        valid_configs.append(config_path)
    
    if skipped:
        logger.debug(f"Skipped configs: {', '.join(skipped)}")
    
    species_configs = []
    for config_path in valid_configs:
        # Extract species name from filename
        species_name = config_path.stem.replace("amalgkit_", "")
        species_configs.append((species_name, config_path))
    
    logger.info(f"Found {len(species_configs)} species configs")
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
        worker_logger.info("Workflow steps: metadata → select → getfastq → quant → merge → curate → sanity")
        
        import time
        workflow_start = time.time()
        return_codes = execute_workflow(cfg, check=False)
        workflow_time = time.time() - workflow_start
        
        # Check results
        success = all(code == 0 or code == 126 for code in return_codes)  # 126 = skipped due to deps
        failed_steps = [i for i, code in enumerate(return_codes) if code != 0 and code != 126]
        
        worker_logger.info("="*80)
        worker_logger.info(f"Workflow completed in {workflow_time:.1f}s")
        if success:
            worker_logger.info(f"✅ {species_name}: SUCCESS")
            message = f"Completed successfully in {workflow_time:.1f}s"
        else:
            worker_logger.warning(f"⚠️  {species_name}: Some steps failed")
            message = f"Failed steps: {failed_steps}"
        worker_logger.info("="*80)
        
        return success, message
        
    except Exception as e:
        worker_logger.error(f"❌ {species_name}: Exception: {e}")
        return False, f"Exception: {e}"


def initialize_tracker_for_all_species(species_configs: list[tuple[str, Path]], tracker) -> None:
    """Initialize progress tracker for all species.
    
    Args:
        species_configs: List of (species_name, config_path) tuples
        tracker: ProgressTracker instance
    """
    import time
    
    total_species = len(species_configs)
    total_samples = 0
    start_time = time.time()
    
    logger.info(f"Initializing progress tracker for {total_species} species...")
    
    for idx, (species_name, config_path) in enumerate(species_configs, 1):
        species_start = time.time()
        logger.info(f"  [{idx}/{total_species}] Initializing {species_name}...")
        
        try:
            cfg = load_workflow_config(config_path)
            metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
            if not metadata_file.exists():
                metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
            
            if metadata_file.exists():
                rows = list(read_delimited(metadata_file, delimiter="\t"))
                sample_ids = [row.get("run") for row in rows if row.get("run")]
                
                if sample_ids:
                    tracker.initialize_species(species_name, len(sample_ids), sample_ids)
                    total_samples += len(sample_ids)
                    
                    logger.info(f"    Found {len(sample_ids)} samples, analyzing current state...")
                    state_start = time.time()
                    
                    # Check current state and populate accordingly
                    status = analyze_species_status(config_path)
                    categories = status.get("categories", {})
                    
                    state_time = time.time() - state_start
                    
                    # Mark completed samples
                    completed = len(categories.get("quantified_and_deleted", []))
                    for sample_id in categories.get("quantified_and_deleted", []):
                        tracker.on_quant_complete(species_name, sample_id)
                        tracker.on_delete_complete(species_name, sample_id)
                    
                    # Mark samples that need delete
                    needs_delete = len(categories.get("quantified_not_deleted", []))
                    for sample_id in categories.get("quantified_not_deleted", []):
                        tracker.on_quant_complete(species_name, sample_id)
                    
                    # Mark downloading samples
                    downloading = len(categories.get("downloading", []))
                    for sample_id in categories.get("downloading", []):
                        tracker.on_download_start(species_name, sample_id)
                    
                    # Mark failed downloads
                    failed = len(categories.get("failed_download", []))
                    for sample_id in categories.get("failed_download", []):
                        tracker.on_download_failed(species_name, sample_id)
                    
                    species_time = time.time() - species_start
                    logger.info(
                        f"    ✓ {species_name}: {len(sample_ids)} samples "
                        f"(completed: {completed}, needs_delete: {needs_delete}, "
                        f"downloading: {downloading}, failed: {failed}) "
                        f"in {species_time:.1f}s (state analysis: {state_time:.1f}s)"
                    )
                else:
                    logger.warning(f"    ⚠️  {species_name}: No samples found in metadata")
            else:
                logger.warning(f"    ⚠️  {species_name}: Metadata file not found")
        except Exception as e:
            logger.error(f"    ❌ {species_name}: Failed to initialize - {e}")
    
    elapsed = time.time() - start_time
    logger.info(
        f"✓ Initialized {total_species} species, {total_samples} total samples in {elapsed:.1f}s"
    )


def dashboard_update_worker(tracker, interval: int = 30):
    """Background worker to periodically update dashboard.
    
    Args:
        tracker: ProgressTracker instance
        interval: Update interval in seconds
    """
    import time
    update_count = 0
    while True:
        try:
            time.sleep(interval)
            update_start = time.time()
            tracker.update_dashboard()
            update_time = time.time() - update_start
            update_count += 1
            
            if update_time > 1.0:
                logger.debug(f"Dashboard update #{update_count} took {update_time:.1f}s (slow)")
            else:
                logger.debug(f"Dashboard update #{update_count} completed in {update_time:.2f}s")
        except Exception as e:
            logger.warning(f"Dashboard update error: {e}")


def main():
    """Main entry point for parallel species workflow execution."""
    import argparse
    import threading
    
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
    
    # Initialize progress tracker
    tracker = get_tracker()
    initialize_tracker_for_all_species(species_configs, tracker)
    
    # Start dashboard update thread
    dashboard_thread = threading.Thread(
        target=dashboard_update_worker,
        args=(tracker, 30),
        daemon=True
    )
    dashboard_thread.start()
    logger.info("Started dashboard update thread (30s interval)")
    
    # Initial dashboard update
    tracker.update_dashboard()
    
    print("\n" + "="*80)
    print("PARALLEL MULTI-SPECIES AMALGKIT WORKFLOW")
    print("="*80)
    print(f"Date: {datetime.now()}")
    print(f"Species: {len(species_configs)}")
    print(f"Progress dashboard: output/progress_dashboard.txt")
    print("="*80)
    print(f"Threads per species: {threads_per_species}")
    print(f"Total threads: {len(species_configs) * threads_per_species}")
    print(f"Execution mode: Parallel (all species simultaneously)")
    print("="*80)
    print()
    
    for name, config in species_configs:
        print(f"  - {name} ({config})")
    print()
    
    # Run all species in parallel
    logger.info(f"Starting parallel execution with {len(species_configs)} workers...")
    logger.info("Submitting species workflows to executor...")
    for idx, (name, config) in enumerate(species_configs, 1):
        logger.info(f"  [{idx}/{len(species_configs)}] Submitting {name}...")
    
    results = {}
    
    logger.info("All species submitted, starting parallel execution...")
    with ProcessPoolExecutor(max_workers=len(species_configs)) as executor:
        # Submit all species workflows (convert Path to str for multiprocessing)
        future_to_species = {
            executor.submit(run_species_workflow_worker, name, str(config), threads_per_species): name
            for name, config in species_configs
        }
        
        # Collect results as they complete with periodic status updates
        import time
        last_status_update = time.time()
        completed_count = 0
        in_progress = len(species_configs)
        
        for future in as_completed(future_to_species):
            species_name = future_to_species[future]
            try:
                success, message = future.result()
                results[species_name] = {"success": success, "message": message}
                completed_count += 1
                in_progress -= 1
                status = "✅" if success else "❌"
                print(f"{status} {species_name}: {message}")
                logger.info(f"{status} {species_name}: {message} ({completed_count}/{len(species_configs)} completed)")
            except Exception as e:
                results[species_name] = {"success": False, "message": f"Exception: {e}"}
                completed_count += 1
                in_progress -= 1
                logger.error(f"❌ {species_name}: Exception: {e}")
                print(f"❌ {species_name}: Exception: {e}")
            
            # Periodic status updates every 60 seconds
            current_time = time.time()
            if current_time - last_status_update >= 60:
                remaining = len(species_configs) - completed_count
                logger.info(
                    f"Progress update: {completed_count} completed, {in_progress} in progress, "
                    f"{remaining} remaining"
                )
                last_status_update = current_time
    
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
    
    # Final dashboard update
    tracker.update_dashboard()
    logger.info(f"Progress dashboard: output/progress_dashboard.txt")
    
    # Exit code
    if all(r["success"] for r in results.values()):
        print("\n✅ All species processed successfully!")
        sys.exit(0)
    else:
        print("\n⚠️  Some species had failures - check logs")
        sys.exit(1)


if __name__ == "__main__":
    main()

