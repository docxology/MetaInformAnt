#!/usr/bin/env python3
"""
Quantify samples that have already been downloaded but not yet quantified.
This script finds FASTQ files and runs quantification, then deletes the FASTQs.
"""

import os
import sys
from pathlib import Path
from datetime import datetime

# Ensure virtual environment is activated before imports
def ensure_venv_activated():
    """Automatically activate virtual environment if it exists and we're not using it."""
    repo_root = Path(__file__).parent.parent.parent.resolve()
    venv_python = repo_root / ".venv" / "bin" / "python3"
    venv_dir = repo_root / ".venv"
    
    # Check if we're already running with venv Python
    current_python = Path(sys.executable)
    
    # Check if current Python is inside .venv directory
    try:
        current_python.relative_to(repo_root / ".venv")
        # We're already using venv Python - ensure environment variables are set
        if "VIRTUAL_ENV" not in os.environ:
            os.environ["VIRTUAL_ENV"] = str(venv_dir)
            venv_bin = str(venv_dir / "bin")
            if venv_bin not in os.environ.get("PATH", ""):
                os.environ["PATH"] = f"{venv_bin}:{os.environ.get('PATH', '')}"
        return
    except ValueError:
        # Not using venv Python - need to switch
        pass
    
    if venv_python.exists():
        # Set up environment variables BEFORE re-exec
        new_env = os.environ.copy()
        new_env["VIRTUAL_ENV"] = str(venv_dir)
        venv_bin = str(venv_dir / "bin")
        new_env["PATH"] = f"{venv_bin}:{new_env.get('PATH', '')}"
        new_env.pop("PYTHONHOME", None)
        
        # Re-exec this script using venv Python with proper environment
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
        print("Continuing with system Python (amalgkit should be installed)...")
        print("=" * 80)
        print()

# Activate venv BEFORE any other imports
ensure_venv_activated()

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.amalgkit import check_cli_available
from metainformant.rna.steps import quantify_sample, delete_sample_fastqs
from metainformant.rna import find_unquantified_samples
from metainformant.core.logging import get_logger
from metainformant.core.io import read_delimited

logger = get_logger("quant_downloaded")


def quantify_samples(config_path: Path, species_name: str):
    """Quantify all downloaded but unquantified samples for a species."""
    logger.info("=" * 80)
    logger.info(f"QUANTIFYING DOWNLOADED SAMPLES: {species_name}")
    logger.info("=" * 80)
    
    # Load config
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Find unquantified samples using metainformant function
    unquantified = find_unquantified_samples(config_path)
    
    if not unquantified:
        logger.info(f"✅ No downloaded samples need quantification for {species_name}")
        return 0, 0
    
    logger.info(f"Found {len(unquantified)} samples with FASTQs needing quantification:")
    for sample_id in unquantified[:10]:
        logger.info(f"  - {sample_id}")
    if len(unquantified) > 10:
        logger.info(f"  ... and {len(unquantified) - 10} more")
    
    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        logger.error(f"❌ No metadata file found in {cfg.work_dir / 'metadata'}")
        return 0, len(unquantified)
    
    rows = list(read_delimited(metadata_file, delimiter="\t"))
    
    # Get quant params from config
    quant_params = dict(cfg.per_step.get("quant", {}))
    quant_params["out_dir"] = str(quant_dir.absolute())
    quant_params["threads"] = cfg.threads
    
    # Inject index_dir if needed (amalgkit quant uses --index_dir, not --genome_dir)
    if "index_dir" not in quant_params and "index-dir" not in quant_params:
        # Try to find index directory - default is out_dir/index/
        index_dir = quant_dir.parent / "work" / "index"
        if not index_dir.exists():
            # Fallback: check if genome_dir has index files
            genome_dir = Path(cfg.genome.get("dest_dir", cfg.work_dir.parent / "genome"))
            if genome_dir.exists():
                # Check if there's an index subdirectory
                potential_index = genome_dir / "index"
                if potential_index.exists():
                    index_dir = potential_index
        if index_dir.exists():
            quant_params["index_dir"] = str(index_dir.absolute())
        else:
            # If no index, try to use fasta_dir with build_index
            fasta_dir = quant_dir.parent / "work" / "fasta"
            if not fasta_dir.exists():
                fasta_dir = quant_dir.parent / "fasta"
            if fasta_dir.exists():
                quant_params["fasta_dir"] = str(fasta_dir.absolute())
                quant_params["build_index"] = True
    
    success_count = 0
    failed_count = 0
    
    # Process each sample using metainformant functions
    for idx, sample_id in enumerate(unquantified, 1):
        logger.info(f"\n[{idx}/{len(unquantified)}] Processing {sample_id}")
        
        try:
            # Get metadata rows for this sample
            sample_rows = [row for row in rows if row.get("run") == sample_id]
            
            if not sample_rows:
                logger.warning(f"  ⚠️  {sample_id} not in metadata, skipping")
                failed_count += 1
                continue
            
            # Use metainformant quantify_sample function
            logger.info(f"  🔬 Quantifying {sample_id}...")
            success, message, abundance_file = quantify_sample(
                sample_id=sample_id,
                metadata_rows=sample_rows,
                quant_params=quant_params,
                log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
                step_name=f"quant_{sample_id}",
            )
            
            if success:
                logger.info(f"  ✅ {message}")
                success_count += 1
            else:
                logger.error(f"  ❌ {message}")
                failed_count += 1
            
            # Always delete FASTQs to free space (even if quant failed)
            logger.info(f"  🗑️  Deleting FASTQs for {sample_id}...")
            delete_sample_fastqs(sample_id, fastq_dir)
        
        except Exception as e:
            logger.error(f"  ❌ Error processing {sample_id}: {e}", exc_info=True)
            failed_count += 1
            # Attempt cleanup
            try:
                delete_sample_fastqs(sample_id, fastq_dir)
            except Exception:
                pass
    
    logger.info("\n" + "=" * 80)
    logger.info(f"QUANTIFICATION SUMMARY: {species_name}")
    logger.info(f"  Total samples: {len(unquantified)}")
    logger.info(f"  Success: {success_count}")
    logger.info(f"  Failed: {failed_count}")
    logger.info("=" * 80)
    
    return success_count, failed_count


def main():
    """Quantify all downloaded samples across all species."""
    from glob import glob
    
    repo_root = Path(__file__).parent.parent.parent.resolve()
    config_dir = repo_root / "config" / "amalgkit"
    
    # Try new location first, fall back to old location
    if not config_dir.exists() or not list(config_dir.glob("amalgkit_*.yaml")):
        config_dir = repo_root / "config"
    
    # Discover all config files
    config_pattern = str(config_dir / "amalgkit_*.yaml")
    config_files = sorted(glob(config_pattern))
    
    species_configs = []
    for config_file in config_files:
        path = Path(config_file)
        # Skip template
        if "template" in path.stem.lower():
            continue
        
        # Extract species name from filename: amalgkit_<species>.yaml
        # Convert to display name: amalgkit_cfloridanus -> C. floridanus
        species_code = path.stem.replace("amalgkit_", "")
        display_name = species_code.replace("_", " ").title()
        species_configs.append((display_name, path))
    
    if not species_configs:
        logger.warning("⚠️  No species configs found")
        return 1
    
    print("\n" + "=" * 80)
    print("QUANTIFY DOWNLOADED SAMPLES")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
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
    
    total_success = 0
    total_failed = 0
    
    for species_name, config_path in species_configs:
        if not config_path.exists():
            logger.warning(f"⚠️  Config not found: {config_path}")
            continue
        
        success, failed = quantify_samples(config_path, species_name)
        total_success += success
        total_failed += failed
    
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total quantified: {total_success}")
    print(f"Total failed: {total_failed}")
    print("=" * 80)
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

