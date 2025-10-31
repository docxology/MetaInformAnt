#!/usr/bin/env python3
"""
Quantify samples that have already been downloaded but not yet quantified.
This script finds FASTQ files and runs quantification, then deletes the FASTQs.
"""

import sys
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.amalgkit import run_amalgkit
from metainformant.core.logging import get_logger
from metainformant.core.io import read_delimited, write_delimited
import shutil

logger = get_logger("quant_downloaded")


def find_unquantified_samples(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples with FASTQs that haven't been quantified."""
    samples_with_fastq = set()
    
    # Find all samples with FASTQ files
    if not fastq_dir.exists():
        return []
    
    for fastq_file in fastq_dir.glob("getfastq/*/*.fastq*"):
        sample_id = fastq_file.parent.name
        samples_with_fastq.add(sample_id)
    
    # Filter to samples not yet quantified
    unquantified = []
    for sample_id in sorted(samples_with_fastq):
        abundance_file = quant_dir / sample_id / "abundance.tsv"
        if not abundance_file.exists():
            unquantified.append(sample_id)
    
    return unquantified


def delete_sample_fastqs(sample_id: str, fastq_dir: Path):
    """Delete FASTQ files for a sample."""
    sample_dir = fastq_dir / "getfastq" / sample_id
    if sample_dir.exists():
        try:
            shutil.rmtree(sample_dir)
            logger.info(f"  🗑️  Deleted FASTQs for {sample_id}")
        except Exception as e:
            logger.warning(f"  ⚠️  Failed to delete {sample_dir}: {e}")


def quantify_samples(config_path: Path, species_name: str):
    """Quantify all downloaded but unquantified samples for a species."""
    logger.info("=" * 80)
    logger.info(f"QUANTIFYING DOWNLOADED SAMPLES: {species_name}")
    logger.info("=" * 80)
    
    # Load config
    cfg = load_workflow_config(config_path)
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))
    
    # Find unquantified samples
    unquantified = find_unquantified_samples(fastq_dir, quant_dir)
    
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
    
    # Inject genome_dir if needed
    if "index" not in quant_params and "genome_dir" not in quant_params:
        genome_dir = Path(cfg.genome.get("dest_dir", cfg.work_dir.parent / "genome"))
        if genome_dir.exists():
            quant_params["genome_dir"] = str(genome_dir)
    
    success_count = 0
    failed_count = 0
    
    # Process each sample
    for idx, sample_id in enumerate(unquantified, 1):
        logger.info(f"\n[{idx}/{len(unquantified)}] Processing {sample_id}")
        
        try:
            # Create temp metadata with just this sample
            sample_rows = [row for row in rows if row.get("run") == sample_id]
            
            if not sample_rows:
                logger.warning(f"  ⚠️  {sample_id} not in metadata, skipping")
                failed_count += 1
                continue
            
            temp_metadata = cfg.work_dir / f"metadata.quant.{sample_id}.tsv"
            write_delimited(sample_rows, temp_metadata, delimiter="\t")
            
            # Run quantification
            logger.info(f"  🔬 Quantifying {sample_id}...")
            quant_params_single = quant_params.copy()
            quant_params_single["metadata"] = str(temp_metadata.absolute())
            
            result = run_amalgkit(
                "quant",
                quant_params_single,
                work_dir=None,
                log_dir=cfg.log_dir or (cfg.work_dir / "logs"),
                step_name=f"quant_{sample_id}",
                check=False,
            )
            
            # Clean up temp metadata
            try:
                temp_metadata.unlink()
            except Exception:
                pass
            
            if result.returncode == 0:
                logger.info(f"  ✅ Quantified {sample_id}")
                success_count += 1
                
                # Delete FASTQs
                logger.info(f"  🗑️  Deleting FASTQs for {sample_id}...")
                delete_sample_fastqs(sample_id, fastq_dir)
            else:
                logger.error(f"  ❌ Quantification failed for {sample_id} (code {result.returncode})")
                failed_count += 1
                # Still delete FASTQs to free space
                delete_sample_fastqs(sample_id, fastq_dir)
        
        except Exception as e:
            logger.error(f"  ❌ Error processing {sample_id}: {e}")
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
    config_dir = Path("config/amalgkit")
    
    # Try new location first, fall back to old location
    if not config_dir.exists() or not list(config_dir.glob("amalgkit_*.yaml")):
        config_dir = Path("config")
    
    species_configs = [
        ("C. floridanus", config_dir / "amalgkit_cfloridanus.yaml"),
        ("M. pharaonis", config_dir / "amalgkit_mpharaonis.yaml"),
        ("P. barbatus", config_dir / "amalgkit_pbarbatus.yaml"),
        ("S. invicta", config_dir / "amalgkit_sinvicta.yaml"),
    ]
    
    print("\n" + "=" * 80)
    print("QUANTIFY DOWNLOADED SAMPLES")
    print("=" * 80)
    print(f"Date: {datetime.now()}")
    print("=" * 80 + "\n")
    
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

