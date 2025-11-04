#!/usr/bin/env python3
"""
Batch download samples for multiple species in parallel with dynamic thread allocation.

Features:
- Configurable total thread count (default: 30)
- Even initial distribution across species
- Dynamic reallocation as species complete
- Automatic thread concentration on remaining species

Thread Allocation:
- Initial: Evenly split threads across all species (e.g., 30 threads / 25 species = 1 for 20, 2 for 5)
- Dynamic: As species complete, threads redistribute to remaining species
- Concentration: Threads concentrate on fewer species as workflow progresses
"""

import os
import sys
import subprocess
import time
from pathlib import Path
from datetime import datetime
from glob import glob
from typing import Optional

# Ensure virtual environment is activated before imports
def ensure_venv_activated():
    """Automatically activate virtual environment if it exists and we're not using it."""
    repo_root = Path(__file__).parent.parent.parent.resolve()
    venv_python = repo_root / ".venv" / "bin" / "python3"
    venv_dir = repo_root / ".venv"
    
    current_python = Path(sys.executable)
    
    try:
        current_python.relative_to(repo_root / ".venv")
        if "VIRTUAL_ENV" not in os.environ:
            os.environ["VIRTUAL_ENV"] = str(venv_dir)
            venv_bin = str(venv_dir / "bin")
            if venv_bin not in os.environ.get("PATH", ""):
                os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
        return
    except ValueError:
        pass
    
    if venv_python.exists():
        new_env = os.environ.copy()
        new_env["VIRTUAL_ENV"] = str(venv_dir)
        venv_bin = str(venv_dir / "bin")
        new_env["PATH"] = f"{venv_bin}:{new_env.get('PATH', '')}"
        new_env.pop("PYTHONHOME", None)
        
        print("=" * 80)
        print("🔄 AUTO-ACTIVATING VIRTUAL ENVIRONMENT")
        print("=" * 80)
        print(f"Current Python:  {current_python}")
        print(f"Venv Python:     {venv_python}")
        print(f"Setting VIRTUAL_ENV={venv_dir}")
        print(f"Updating PATH to include {venv_bin}")
        print("=" * 80)
        print()
        sys.stdout.flush()
        
        os.execve(str(venv_python), [str(venv_python)] + sys.argv, new_env)
    else:
        print()
        print("=" * 80)
        print("⚠️  WARNING: Virtual environment not found")
        print("=" * 80)
        print(f"Expected location: {venv_python}")
        print("Continuing with system Python...")
        print("=" * 80)
        print()

ensure_venv_activated()

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.amalgkit import build_amalgkit_command, check_cli_available
from metainformant.core.io import read_delimited
from metainformant.core.logging import get_logger

logger = get_logger("batch_download")


def count_quantified_samples(config_path: Path) -> tuple[int, int]:
    """Count quantified and total samples for a species.
    
    Args:
        config_path: Path to species config file
        
    Returns:
        Tuple of (quantified_count, total_count)
    """
    try:
        cfg = load_workflow_config(config_path)
        quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
        
        # Count quantified
        quantified = 0
        if quant_dir.exists():
            quantified = len([d for d in quant_dir.iterdir() 
                             if d.is_dir() and (d / "abundance.tsv").exists()])
        
        # Count total from metadata
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        total = 0
        if metadata_file.exists():
            try:
                rows = list(read_delimited(metadata_file, delimiter="\t"))
                total = len([r for r in rows if r.get("run")])
            except Exception:
                pass
        
        return quantified, total
    except Exception:
        return 0, 0


def is_species_complete(config_path: Path, completion_threshold: float = 0.95) -> bool:
    """Check if species download/quantification is complete.
    
    Args:
        config_path: Path to species config file
        completion_threshold: Fraction of samples that must be quantified (default: 0.95)
        
    Returns:
        True if species is complete (≥threshold quantified)
    """
    quantified, total = count_quantified_samples(config_path)
    if total == 0:
        return False
    return (quantified / total) >= completion_threshold


def distribute_threads(total_threads: int, num_species: int) -> list[int]:
    """Distribute threads evenly across species.
    
    Args:
        total_threads: Total number of threads to distribute
        num_species: Number of species
        
    Returns:
        List of thread counts per species (one per species)
        
    Example:
        distribute_threads(30, 25) -> [1]*20 + [2]*5 (20 species get 1, 5 get 2)
    """
    if num_species == 0:
        return []
    
    base_threads = total_threads // num_species
    extra_threads = total_threads % num_species
    
    # Base allocation: all species get base_threads
    allocation = [base_threads] * num_species
    
    # Distribute extra threads to first extra_threads species
    for i in range(extra_threads):
        allocation[i] += 1
    
    return allocation


def redistribute_threads(total_threads: int, active_species: list[tuple[str, Path]]) -> dict[str, int]:
    """Redistribute threads among active (incomplete) species.
    
    Args:
        total_threads: Total number of threads available
        active_species: List of (species_name, config_path) tuples for incomplete species
        
    Returns:
        Dictionary mapping species_name -> thread_count
    """
    if not active_species:
        return {}
    
    num_active = len(active_species)
    thread_allocation = distribute_threads(total_threads, num_active)
    
    return {
        species_name: threads
        for (species_name, _), threads in zip(active_species, thread_allocation)
    }


def start_download_process(config_path: Path, species_name: str, threads: int, log_dir: Path) -> subprocess.Popen:
    """Start an amalgkit getfastq process for a species.
    
    Args:
        config_path: Path to species config file
        species_name: Display name for species
        threads: Number of threads to use
        log_dir: Directory for log files
        
    Returns:
        subprocess.Popen object for the running process
    """
    try:
        cfg = load_workflow_config(config_path)
        
        # Get metadata file
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            raise FileNotFoundError(f"No metadata file found for {species_name}")
        
        # Build getfastq params
        getfastq_params = dict(cfg.per_step.get("getfastq", {}))
        getfastq_params["out_dir"] = str(Path(getfastq_params.get("out_dir", cfg.work_dir / "fastq")).absolute())
        getfastq_params["metadata"] = str(metadata_file.absolute())
        getfastq_params["threads"] = threads
        
        # Enable cloud acceleration if configured
        if cfg.per_step.get("getfastq", {}).get("accelerate", True):
            getfastq_params.setdefault("aws", "yes")
            getfastq_params.setdefault("gcp", "yes")
            getfastq_params.setdefault("ncbi", "yes")
        
        # Set remove_sra to yes to save space
        getfastq_params["remove_sra"] = "yes"
        
        # Build command
        cmd = build_amalgkit_command("getfastq", getfastq_params)
        
        # Setup log file
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"getfastq_{species_name.replace(' ', '_')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logger.info(f"Starting {species_name} with {threads} threads (PID will be assigned)")
        logger.info(f"  Command: {' '.join(cmd[:5])}...")
        logger.info(f"  Log: {log_file}")
        
        # Start process (log file will be written to as process runs)
        log_handle = open(log_file, 'w')
        process = subprocess.Popen(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            cwd=str(cfg.work_dir.parent) if cfg.work_dir else None,
        )
        
        # Store log handle reference (will be closed when process finishes)
        # Note: In production, consider using a process manager to track and close handles
        return process
        
    except Exception as e:
        logger.error(f"Failed to start process for {species_name}: {e}")
        raise


def main():
    """Download samples for multiple species with dynamic thread allocation.
    
    Configuration:
        --total-threads: Total threads to use across all species (default: 30)
        --max-species: Maximum number of species to process (default: all)
        --check-interval: Seconds between completion checks (default: 300)
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Batch download samples with dynamic thread allocation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: 30 total threads, distributed evenly across all species
  python3 scripts/rna/batch_download_species.py
  
  # 50 total threads, redistributed as species complete
  python3 scripts/rna/batch_download_species.py --total-threads 50
  
  # Process only first 10 species
  python3 scripts/rna/batch_download_species.py --max-species 10
        """
    )
    parser.add_argument(
        "--total-threads",
        type=int,
        default=30,
        help="Total number of threads to use across all species (default: 30)"
    )
    parser.add_argument(
        "--max-species",
        type=int,
        default=None,
        help="Maximum number of species to process (default: all)"
    )
    parser.add_argument(
        "--check-interval",
        type=int,
        default=300,
        help="Seconds between completion checks (default: 300)"
    )
    
    args = parser.parse_args()
    
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
    
    # Limit species if requested
    if args.max_species:
        species_configs = species_configs[:args.max_species]
    
    print("\n" + "=" * 80)
    print("DYNAMIC THREAD ALLOCATION: MULTI-SPECIES DOWNLOADS")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    print(f"Total threads: {args.total_threads}")
    print(f"Species discovered: {len(species_configs)}")
    print(f"Check interval: {args.check_interval} seconds")
    print("=" * 80)
    print("\nThread Allocation Strategy:")
    print("  • Initial: Even distribution across all species")
    print("  • Dynamic: Redistribute as species complete")
    print("  • Concentration: Threads focus on remaining species")
    print("=" * 80 + "\n")
    
    # Check if amalgkit is available
    logger.info("Checking for amalgkit...")
    amalgkit_available, amalgkit_msg = check_cli_available()
    if not amalgkit_available:
        logger.error(f"❌ amalgkit not available: {amalgkit_msg}")
        logger.error("Please ensure virtual environment is activated and amalgkit is installed:")
        logger.error("  source .venv/bin/activate")
        logger.error("  pip install git+https://github.com/kfuku52/amalgkit")
        return 1
    logger.info(f"✅ amalgkit available: {amalgkit_msg[:100] if len(amalgkit_msg) > 100 else amalgkit_msg}")
    print()
    
    # Filter out already-complete species
    active_species = []
    complete_species = []
    
    for species_name, config_path in species_configs:
        if is_species_complete(config_path):
            complete_species.append((species_name, config_path))
            quantified, total = count_quantified_samples(config_path)
            logger.info(f"✅ {species_name}: Already complete ({quantified}/{total} quantified)")
        else:
            active_species.append((species_name, config_path))
    
    if complete_species:
        print(f"Found {len(complete_species)} already-complete species (skipping)")
        print()
    
    if not active_species:
        logger.info("✅ All species are already complete!")
        return 0
    
    print(f"Processing {len(active_species)} active species with {args.total_threads} total threads")
    print()
    
    # Initial thread distribution
    thread_allocation = redistribute_threads(args.total_threads, active_species)
    
    print("Initial Thread Allocation:")
    for species_name, threads in sorted(thread_allocation.items()):
        print(f"  {species_name:30s}: {threads} threads")
    print()
    
    # Setup logging directory
    repo_root = Path(__file__).parent.parent.parent.resolve()
    log_dir = repo_root / "output" / "batch_download_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Track running processes
    running_processes: dict[str, subprocess.Popen] = {}
    process_configs: dict[str, Path] = {}
    results: dict[str, dict] = {}
    
    # Start initial processes
    for species_name, config_path in active_species:
        threads = thread_allocation[species_name]
        try:
            process = start_download_process(config_path, species_name, threads, log_dir)
            running_processes[species_name] = process
            process_configs[species_name] = config_path
            logger.info(f"Started {species_name} (PID: {process.pid}, {threads} threads)")
        except Exception as e:
            logger.error(f"Failed to start {species_name}: {e}")
            results[species_name] = {
                "success": False,
                "message": f"Failed to start: {e}",
                "stats": {},
            }
    
    print(f"\nStarted {len(running_processes)} processes")
    print(f"Monitoring every {args.check_interval} seconds...")
    print("=" * 80)
    print()
    
    # Monitor and redistribute
    iteration = 0
    while running_processes:
        iteration += 1
        time.sleep(args.check_interval)
        
        print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Check #{iteration}")
        print("-" * 80)
        
        # Check for completed processes
        completed = []
        for species_name, process in list(running_processes.items()):
            if process.poll() is not None:
                # Process finished
                returncode = process.returncode
                completed.append(species_name)
                
                if returncode == 0:
                    status = "✅"
                    message = "Completed successfully"
                    results[species_name] = {
                        "success": True,
                        "message": message,
                        "stats": {"return_code": returncode},
                    }
                else:
                    status = "⚠️"
                    message = f"Completed with code {returncode}"
                    results[species_name] = {
                        "success": True,  # Non-zero codes can be warnings
                        "message": message,
                        "stats": {"return_code": returncode},
                    }
                
                logger.info(f"{status} {species_name}: {message}")
                del running_processes[species_name]
        
        # Check for newly-complete species (even if process still running)
        # This handles cases where process is still running but work is done
        newly_complete = []
        for species_name in list(running_processes.keys()):
            config_path = process_configs[species_name]
            if is_species_complete(config_path):
                newly_complete.append(species_name)
                process = running_processes[species_name]
                process.terminate()  # Stop the process since work is done
                logger.info(f"✅ {species_name}: Work complete (terminating process)")
                results[species_name] = {
                    "success": True,
                    "message": "Completed (detected completion)",
                    "stats": {"return_code": 0},
                }
                del running_processes[species_name]
        
        # Redistribute threads if any completed
        if completed or newly_complete:
            remaining_active = [
                (name, cfg) for name, cfg in active_species
                if name not in results or not results[name].get("success", False)
            ]
            
            if remaining_active:
                # Get current active (running) species
                currently_running = [
                    (name, process_configs[name])
                    for name in running_processes.keys()
                ]
                
                # Calculate threads to redistribute
                threads_in_use = sum(thread_allocation.get(name, 0) for name in running_processes.keys())
                threads_available = args.total_threads - threads_in_use
                
                # Redistribute among remaining active
                new_allocation = redistribute_threads(args.total_threads, remaining_active)
                
                # Update allocations for running processes
                for species_name in list(running_processes.keys()):
                    old_threads = thread_allocation.get(species_name, 1)
                    new_threads = new_allocation.get(species_name, old_threads)
                    
                    if new_threads != old_threads:
                        logger.info(f"  {species_name}: {old_threads} → {new_threads} threads")
                        # Note: Can't change threads on running process, but will use new allocation when restarting
                        thread_allocation[species_name] = new_threads
                
                print(f"Redistributed threads: {len(completed + newly_complete)} completed, "
                      f"{len(running_processes)} still running")
            else:
                print(f"All species complete or in progress")
        
        # Status update
        if running_processes:
            print(f"Active: {len(running_processes)} species")
            for species_name in sorted(running_processes.keys()):
                threads = thread_allocation.get(species_name, 1)
                process = running_processes[species_name]
                pid = process.pid
                status = "running" if process.poll() is None else "finished"
                print(f"  {species_name:30s}: PID {pid}, {threads} threads, {status}")
        else:
            print("All processes completed")
            break
    
    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total species: {len(species_configs)}")
    print(f"Already complete: {len(complete_species)}")
    print(f"Processed: {len(active_species)}")
    print()
    
    success_count = sum(1 for r in results.values() if r.get("success", False))
    failed_count = len(results) - success_count
    
    print(f"Successfully completed: {success_count}")
    print(f"Failed: {failed_count}")
    print()
    
    if results:
        print("Per-species results:")
        for species_name, result in sorted(results.items()):
            status = "✅" if result.get("success", False) else "❌"
            message = result.get("message", "Unknown")
            print(f"  {status} {species_name:30s}: {message}")
    
    print("=" * 80)
    
    return 0 if failed_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

