#!/usr/bin/env python3
"""
Comprehensive script to process SRA files: convert to FASTQ → quantify → delete.

This script:
1. Finds SRA files that haven't been quantified
2. Converts SRA to FASTQ using amalgkit getfastq
3. Quantifies the FASTQ files using amalgkit quant
4. Deletes both SRA and FASTQ files after successful quantification
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import shutil

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
from metainformant.rna.amalgkit import check_cli_available
from metainformant.rna.steps import process_sample_pipeline, delete_sample_fastqs
from metainformant.core.logging import get_logger
from glob import glob

logger = get_logger("process_sra")


def find_sra_samples(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples with SRA files that haven't been quantified."""
    samples_with_sra = set()
    
    if not fastq_dir.exists():
        return []
    
    # Check for SRA files in getfastq subdirectory
    for sra_file in fastq_dir.glob("getfastq/*/*.sra"):
        sample_id = sra_file.parent.name
        samples_with_sra.add(sample_id)
    
    # Also check direct structure
    for sra_file in fastq_dir.glob("*/*.sra"):
        sample_id = sra_file.parent.name
        samples_with_sra.add(sample_id)
    
    # Filter to samples not yet quantified
    unquantified = []
    for sample_id in sorted(samples_with_sra):
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if not abundance_file.exists():
            unquantified.append(sample_id)
    
    return unquantified


# delete_sample_files and process_sra_sample now use metainformant.rna.steps functions


def process_species(config_path: Path, species_name: str):
    """Process all SRA samples for a species."""
    logger.info("=" * 80)
    logger.info(f"PROCESSING SRA SAMPLES: {species_name}")
    logger.info("=" * 80)
    
    # Load config
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Find SRA samples
    sra_samples = find_sra_samples(fastq_dir, quant_dir)
    
    if not sra_samples:
        logger.info(f"✅ No SRA samples need processing for {species_name}")
        return 0, 0
    
    logger.info(f"Found {len(sra_samples)} SRA samples needing processing:")
    for sample_id in sra_samples[:10]:
        logger.info(f"  - {sample_id}")
    if len(sra_samples) > 10:
        logger.info(f"  ... and {len(sra_samples) - 10} more")
    
    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
    
    if not metadata_file.exists():
        logger.error(f"❌ No metadata file found in {cfg.work_dir / 'metadata'}")
        return 0, len(sra_samples)
    
    rows = list(read_delimited(metadata_file, delimiter="\t"))
    
    success_count = 0
    failed_count = 0
    
    # Process each sample
    for idx, sample_id in enumerate(sra_samples, 1):
        logger.info(f"\n[{idx}/{len(sra_samples)}] Processing {sample_id}")
        
        success, msg = process_sra_sample(
            sample_id,
            config_path,
            species_name,
            fastq_dir,
            quant_dir,
            metadata_file,
            rows,
        )
        
        if success:
            success_count += 1
        else:
            logger.error(f"  ❌ Failed: {msg}")
            failed_count += 1
            # Still try to delete to free space (might have partial extraction)
            try:
                delete_sample_files(sample_id, fastq_dir)
            except Exception:
                pass
    
    logger.info("\n" + "=" * 80)
    logger.info(f"PROCESSING SUMMARY: {species_name}")
    logger.info(f"  Total samples: {len(sra_samples)}")
    logger.info(f"  Success: {success_count}")
    logger.info(f"  Failed: {failed_count}")
    logger.info("=" * 80)
    
    return success_count, failed_count


def main():
    """Process SRA samples across all species."""
    from glob import glob
    
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
        logger.warning("⚠️  No species configs found")
        return 1
    
    print("\n" + "=" * 80)
    print("PROCESS SRA SAMPLES: CONVERT → QUANTIFY → DELETE")
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
        
        success, failed = process_species(config_path, species_name)
        total_success += success
        total_failed += failed
    
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total processed: {total_success}")
    print(f"Total failed: {total_failed}")
    print("=" * 80)
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())


