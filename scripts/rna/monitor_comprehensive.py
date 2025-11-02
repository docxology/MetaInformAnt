#!/usr/bin/env python3
"""
Monitor comprehensive multi-species RNA-seq workflow progress.

Tracks download, quantification, and completion status across all 4 species.
"""

import sys
from pathlib import Path
from datetime import datetime
import time

# Species configurations
SPECIES = {
    'cfloridanus': {'name': 'C. floridanus', 'total': 307, 'batches': 25},
    'pbarbatus': {'name': 'P. barbatus', 'total': 83, 'batches': 7},
    'mpharaonis': {'name': 'M. pharaonis', 'total': 100, 'batches': 9},
    'sinvicta': {'name': 'S. invicta', 'total': 354, 'batches': 30},
}

def get_quantified_count(species):
    """Count quantified samples (abundance.tsv files)."""
    quant_dir = Path(f"output/amalgkit/{species}/quant")
    if not quant_dir.exists():
        return 0
    return len(list(quant_dir.glob("*/abundance.tsv")))

def get_downloading_count(species):
    """Count samples currently being downloaded."""
    fastq_dir = Path(f"output/amalgkit/{species}/fastq")
    if not fastq_dir.exists():
        return 0
    return len(list(fastq_dir.glob("SRR*")))

def get_fastq_size(species):
    """Get total size of FASTQ directory."""
    import subprocess
    fastq_dir = Path(f"output/amalgkit/{species}/fastq")
    if not fastq_dir.exists():
        return "0"
    result = subprocess.run(
        ["du", "-sh", str(fastq_dir)],
        capture_output=True,
        text=True
    )
    if result.returncode == 0:
        return result.stdout.split()[0]
    return "?"

def get_current_batch(species):
    """Parse log to find current batch number."""
    log_files = list(Path("output").glob(f"workflow_{species}_*.log"))
    if not log_files:
        return "?"
    
    log_file = sorted(log_files)[-1]  # Most recent
    try:
        with open(log_file) as f:
            for line in f:
                if "📦 BATCH" in line:
                    # Extract "BATCH X/Y"
                    parts = line.split("BATCH")[1].split(":")[0].strip()
                    return parts
    except:
        pass
    return "?"

def main():
    print("\n" + "="*80)
    print("  MULTI-SPECIES RNA-SEQ WORKFLOW MONITOR")
    print("="*80)
    print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    total_samples = 0
    total_quantified = 0
    
    for species_id, info in SPECIES.items():
        quantified = get_quantified_count(species_id)
        downloading = get_downloading_count(species_id)
        fastq_size = get_fastq_size(species_id)
        current_batch = get_current_batch(species_id)
        
        total_samples += info['total']
        total_quantified += quantified
        
        percent = (quantified * 100) // info['total']
        
        print(f"📊 {info['name']} ({species_id})")
        print(f"   Progress: {quantified}/{info['total']} ({percent}%)")
        print(f"   Batch: {current_batch}/{info['batches']}")
        print(f"   Downloading: {downloading} samples ({fastq_size})")
        print()
    
    print("="*80)
    overall_percent = (total_quantified * 100) // total_samples
    print(f"🎯 OVERALL: {total_quantified}/{total_samples} samples ({overall_percent}%)")
    print("="*80)
    
    # Estimate remaining time
    remaining = total_samples - total_quantified
    if remaining > 0:
        est_minutes = remaining * 7.5  # 7.5 min/sample
        est_hours = est_minutes / 60
        print(f"\n⏳ Estimated remaining: {est_hours:.1f} hours ({remaining} samples)")
    else:
        print("\n✅ ALL SAMPLES COMPLETE!")
    
    print()

if __name__ == '__main__':
    main()

