#!/usr/bin/env python3
"""
Manually quantify completed downloads and cleanup FASTQs.
Runs sequentially with proper threading to avoid overwhelming the system.
"""

import sys
import subprocess
import shutil
from pathlib import Path
import yaml
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)
logger = logging.getLogger(__name__)


def load_config(config_path: Path) -> dict:
    """Load amalgkit config."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_kallisto_index(species_dir: Path) -> Path:
    """Get kallisto index path."""
    index_dir = species_dir / "work" / "index"
    if index_dir.exists():
        indices = list(index_dir.glob("*.idx"))
        if indices:
            return indices[0]
    return None


def find_complete_unquantified(species_dir: Path) -> list:
    """Find samples with complete FASTQs but no quantification."""
    fastq_dir = species_dir / "fastq"
    quant_dir = species_dir / "quant"
    
    if not fastq_dir.exists():
        return []
    
    ready = []
    for sample_dir in fastq_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        
        # Check if already quantified
        if (quant_dir / sample_id / "abundance.tsv").exists():
            continue
        
        # Check for complete FASTQ files (>1KB, not empty)
        fastq_files = []
        for fq in sample_dir.glob("*.fastq.gz"):
            if fq.stat().st_size > 1024:  # >1KB
                fastq_files.append(fq)
        
        if fastq_files:
            ready.append((sample_id, sample_dir, fastq_files))
    
    return ready


def quantify_sample(sample_id: str, fastq_files: list, kallisto_index: Path, 
                   output_dir: Path, species: str) -> bool:
    """Run kallisto quantification."""
    try:
        quant_output = output_dir / sample_id
        quant_output.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "kallisto", "quant",
            "-i", str(kallisto_index),
            "-o", str(quant_output),
            "-t", "12"  # Use 12 threads for speed
        ]
        
        # Single-end or paired-end
        if len(fastq_files) == 1:
            cmd.extend(["--single", "-l", "200", "-s", "20"])
        
        cmd.extend([str(f) for f in sorted(fastq_files)])
        
        logger.info(f"  Quantifying {species}/{sample_id} ...")
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0 and (quant_output / "abundance.tsv").exists():
            logger.info(f"  ✅ {species}/{sample_id}: Quantified")
            return True
        else:
            logger.error(f"  ❌ {species}/{sample_id}: Failed")
            return False
            
    except Exception as e:
        logger.error(f"  ❌ {species}/{sample_id}: {e}")
        return False


def cleanup_sample(sample_dir: Path, sample_id: str, species: str):
    """Delete FASTQs after successful quantification."""
    try:
        size_mb = sum(f.stat().st_size for f in sample_dir.rglob("*") if f.is_file()) / (1024*1024)
        shutil.rmtree(sample_dir)
        logger.info(f"  🗑️  {species}/{sample_id}: Deleted {size_mb:.1f} MB")
    except Exception as e:
        logger.error(f"  ❌ {species}/{sample_id}: Cleanup failed - {e}")


def process_species(species_name: str, config_path: Path):
    """Process one species."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing {species_name}")
    logger.info(f"{'='*80}")
    
    repo_root = Path(__file__).parent.parent.parent
    species_dir = repo_root / "output" / "amalgkit" / species_name
    
    kallisto_index = get_kallisto_index(species_dir)
    if not kallisto_index:
        logger.error(f"❌ No kallisto index for {species_name}")
        return
    
    samples = find_complete_unquantified(species_dir)
    if not samples:
        logger.info(f"✅ No complete samples need quantification")
        return
    
    logger.info(f"Found {len(samples)} samples ready for quantification")
    
    quantified = 0
    cleaned = 0
    
    # Process sequentially to avoid overwhelming system
    for sample_id, sample_dir, fastq_files in samples:
        if quantify_sample(sample_id, fastq_files, kallisto_index, 
                          species_dir / "quant", species_name):
            quantified += 1
            cleanup_sample(sample_dir, sample_id, species_name)
            cleaned += 1
    
    logger.info(f"\n{species_name} Summary: {quantified} quantified, {cleaned} cleaned")


def main():
    logger.info("="*80)
    logger.info("MANUAL QUANTIFICATION AND CLEANUP")
    logger.info("="*80)
    
    repo_root = Path(__file__).parent.parent.parent
    config_dir = repo_root / "config" / "amalgkit"
    
    species_configs = [
        ("cfloridanus", config_dir / "amalgkit_cfloridanus.yaml"),
        ("pbarbatus", config_dir / "amalgkit_pbarbatus.yaml"),
        ("mpharaonis", config_dir / "amalgkit_mpharaonis.yaml"),
        ("sinvicta", config_dir / "amalgkit_sinvicta.yaml")
    ]
    
    for species_name, config_path in species_configs:
        if config_path.exists():
            process_species(species_name, config_path)
    
    logger.info("\n" + "="*80)
    logger.info("COMPLETE")
    logger.info("="*80)


if __name__ == "__main__":
    main()

