#!/usr/bin/env python3
"""
Batch download samples for multiple species in parallel.

Downloads samples for 4 species simultaneously using amalgkit getfastq.
Each species runs in its own process/thread to maximize throughput.
"""

import os
import sys
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob

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
from metainformant.rna.amalgkit import run_amalgkit, check_cli_available
from metainformant.core.logging import get_logger

logger = get_logger("batch_download")


def download_species_samples(config_path: Path, species_name: str, threads_override: int | None = None) -> tuple[bool, str, dict]:
    """Download samples for a single species using amalgkit getfastq.
    
    Args:
        config_path: Path to species config file
        species_name: Display name for species
        threads_override: Optional override for thread count (defaults to config value)
    
    Returns:
        Tuple of (success: bool, message: str, stats: dict)
    """
    species_logger = get_logger(f"download_{species_name.replace(' ', '_')}")
    
    try:
        species_logger.info(f"Starting download for {species_name}")
        species_logger.info(f"Config: {config_path}")
        
        # Load config
        cfg = load_workflow_config(config_path)
        
        # Get metadata file
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            return False, f"No metadata file found", {}
        
        # Get getfastq params
        getfastq_params = dict(cfg.per_step.get("getfastq", {}))
        getfastq_params["out_dir"] = str(Path(getfastq_params.get("out_dir", cfg.work_dir / "fastq")).absolute())
        getfastq_params["metadata"] = str(metadata_file.absolute())
        getfastq_params["threads"] = threads_override if threads_override is not None else cfg.threads
        
        # Enable cloud acceleration if configured
        if cfg.per_step.get("getfastq", {}).get("accelerate", True):
            getfastq_params.setdefault("aws", "yes")
            getfastq_params.setdefault("gcp", "yes")
            getfastq_params.setdefault("ncbi", "yes")
        
        # Set remove_sra to yes to save space (delete SRA after extraction)
        getfastq_params["remove_sra"] = "yes"
        
        species_logger.info(f"Running amalgkit getfastq for {species_name}...")
        species_logger.info(f"  Output directory: {getfastq_params['out_dir']}")
        species_logger.info(f"  Threads: {getfastq_params['threads']} {'(overridden)' if threads_override else ''}")
        
        # Run getfastq
        result = run_amalgkit(
            "getfastq",
            getfastq_params,
            work_dir=None,
            log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
            step_name=f"getfastq_{species_name.replace(' ', '_')}",
            check=False,
        )
        
        stats = {
            "return_code": result.returncode,
            "species": species_name,
            "config": str(config_path),
        }
        
        if result.returncode == 0:
            species_logger.info(f"✅ Successfully completed download for {species_name}")
            return True, f"Download completed successfully", stats
        else:
            species_logger.warning(f"⚠️  Download completed with warnings (code {result.returncode}) for {species_name}")
            return True, f"Download completed with warnings (code {result.returncode})", stats
        
    except Exception as e:
        species_logger.error(f"❌ Error downloading {species_name}: {e}", exc_info=True)
        return False, str(e), {"species": species_name, "error": str(e)}


def main():
    """Download samples for multiple species in parallel batches.
    
    Configuration:
        --species-count: Number of species to download in parallel (default: 3)
        --threads-per-species: Threads per species download (default: 10)
        --max-species: Maximum number of species to process (default: all)
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Batch download samples for multiple species in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: 3 species, 10 threads each (30 total downloads)
  python3 scripts/rna/batch_download_species.py
  
  # 4 species, 12 threads each (48 total downloads)
  python3 scripts/rna/batch_download_species.py --species-count 4 --threads-per-species 12
  
  # 2 species, 8 threads each (16 total downloads)
  python3 scripts/rna/batch_download_species.py --species-count 2 --threads-per-species 8
        """
    )
    parser.add_argument(
        "--species-count",
        type=int,
        default=3,
        help="Number of species to download in parallel (default: 3)"
    )
    parser.add_argument(
        "--threads-per-species",
        type=int,
        default=10,
        help="Number of threads per species download (default: 10)"
    )
    parser.add_argument(
        "--max-species",
        type=int,
        default=None,
        help="Maximum number of species to process (default: all)"
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
    print("BATCH DOWNLOAD: MULTI-SPECIES PARALLEL DOWNLOADS")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    print(f"Species discovered: {len(species_configs)}")
    print(f"Parallel downloads: {args.species_count} species at once")
    print(f"Threads per species: {args.threads_per_species}")
    print(f"Total concurrent downloads: {args.species_count * args.threads_per_species}")
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
    
    # Process species in batches
    batch_size = args.species_count
    threads_per_species = args.threads_per_species
    total_success = 0
    total_failed = 0
    results = {}
    
    for batch_start in range(0, len(species_configs), batch_size):
        batch = species_configs[batch_start:batch_start + batch_size]
        batch_num = (batch_start // batch_size) + 1
        total_batches = (len(species_configs) + batch_size - 1) // batch_size
        
        print("=" * 80)
        print(f"BATCH {batch_num}/{total_batches}: Downloading {len(batch)} species in parallel")
        print("=" * 80)
        for name, config in batch:
            print(f"  - {name}")
        print("=" * 80)
        print()
        
        # Run downloads in parallel
        with ThreadPoolExecutor(max_workers=len(batch)) as executor:
            futures = {
                executor.submit(download_species_samples, config_path, species_name, threads_per_species): (species_name, config_path)
                for species_name, config_path in batch
            }
            
            # Process results as they complete
            for future in as_completed(futures):
                species_name, config_path = futures[future]
                try:
                    success, message, stats = future.result()
                    
                    if success:
                        total_success += 1
                        status = "✅"
                        logger.info(f"{status} {species_name}: {message}")
                    else:
                        total_failed += 1
                        status = "❌"
                        logger.error(f"{status} {species_name}: {message}")
                    
                    results[species_name] = {
                        "success": success,
                        "message": message,
                        "stats": stats,
                    }
                    
                except Exception as e:
                    total_failed += 1
                    logger.error(f"❌ {species_name}: Exception - {e}", exc_info=True)
                    results[species_name] = {
                        "success": False,
                        "message": f"Exception: {e}",
                        "stats": {},
                    }
        
        print()
        print(f"Batch {batch_num} complete: {len([r for r in results.values() if r['success']])} succeeded, "
              f"{len([r for r in results.values() if not r['success']])} failed")
        print()
    
    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total species: {len(species_configs)}")
    print(f"Successfully downloaded: {total_success}")
    print(f"Failed: {total_failed}")
    print()
    
    if results:
        print("Per-species results:")
        for species_name, result in sorted(results.items()):
            status = "✅" if result["success"] else "❌"
            print(f"  {status} {species_name:30s}: {result['message']}")
    
    print("=" * 80)
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

