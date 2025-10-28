#!/usr/bin/env python3
"""
Batch Download + Quantification for P. barbatus RNA-seq samples.

This script processes remaining samples in batches to manage disk space.
Each batch: Download -> Quantify -> Delete FASTQs -> Repeat

Usage:
    cd /Users/4d/Documents/GitHub/metainformant
    PYTHONPATH=src python3 output/amalgkit/pbarbatus/batch_process.py
"""

from pathlib import Path
import sys
import subprocess
import time
from datetime import datetime

def run_kallisto_quant(sample_id, base_dir):
    """Run kallisto quantification on a single sample."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "getfastq" / sample_id
    index_file = work_dir / "index" / "Pogonomyrmex_barbatus.idx"
    output_dir = work_dir / "quant" / sample_id
    
    # Check if FASTQs exist
    fastq_1 = fastq_dir / f"{sample_id}_1.fastq.gz"
    fastq_2 = fastq_dir / f"{sample_id}_2.fastq.gz"
    
    if not (fastq_1.exists() and fastq_2.exists()):
        print(f"   ⚠️  FASTQs not found for {sample_id}")
        return False
    
    # Run kallisto
    cmd = [
        "kallisto", "quant",
        "-i", str(index_file),
        "-o", str(output_dir),
        "-t", "6",
        str(fastq_1),
        str(fastq_2)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            # Check abundance.tsv was created
            abundance = output_dir / "abundance.tsv"
            if abundance.exists():
                return True
        else:
            print(f"   ⚠️  Kallisto failed for {sample_id}: {result.stderr[:200]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"   ⚠️  Kallisto timeout for {sample_id}")
        return False
    except Exception as e:
        print(f"   ⚠️  Error quantifying {sample_id}: {e}")
        return False
    
    return False

def cleanup_fastq(sample_id, base_dir):
    """Delete FASTQ files for a sample after successful quantification."""
    import shutil
    
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "getfastq" / sample_id
    
    if fastq_dir.exists():
        try:
            size_before = sum(f.stat().st_size for f in fastq_dir.rglob("*") if f.is_file())
            shutil.rmtree(fastq_dir)
            size_gb = size_before / (1024**3)
            print(f"   ✅ Cleaned up {size_gb:.2f}GB")
            return True
        except Exception as e:
            print(f"   ⚠️  Cleanup failed: {e}")
            return False
    return True

def process_batch(batch_num, samples, base_dir, log_file):
    """Process a batch of samples."""
    print("\n" + "=" * 80)
    print(f"BATCH {batch_num}: {len(samples)} samples")
    print("=" * 80)
    
    with open(log_file, 'a') as log:
        log.write(f"\n\n{'=' * 80}\n")
        log.write(f"Batch {batch_num} started at {datetime.now()}\n")
        log.write(f"Samples: {samples}\n")
        log.write(f"{'=' * 80}\n")
    
    downloaded = 0
    quantified = 0
    cleaned = 0
    
    # Download all samples in batch
    print(f"\n📥 Step 1: Downloading {len(samples)} samples...")
    for i, sample in enumerate(samples, 1):
        print(f"   [{i}/{len(samples)}] Downloading {sample}...")
        
        # Use amalgkit getfastq for this specific sample
        # (This would be the actual download command)
        # For now, we'll check if already downloaded
        
        work_dir = base_dir / "work"
        fastq_dir = work_dir / "getfastq" / sample
        fastq_1 = fastq_dir / f"{sample}_1.fastq.gz"
        
        if fastq_1.exists():
            print(f"      ✅ Already downloaded")
            downloaded += 1
        else:
            print(f"      📥 Need to download (not yet implemented)")
            # Here you would call amalgkit.getfastq with batch parameter
    
    # Quantify all samples in batch
    print(f"\n📈 Step 2: Quantifying {len(samples)} samples...")
    for i, sample in enumerate(samples, 1):
        print(f"   [{i}/{len(samples)}] Quantifying {sample}...")
        
        if run_kallisto_quant(sample, base_dir):
            print(f"      ✅ Quantified")
            quantified += 1
            
            # Cleanup FASTQ immediately
            if cleanup_fastq(sample, base_dir):
                cleaned += 1
        else:
            print(f"      ❌ Failed")
    
    # Log results
    with open(log_file, 'a') as log:
        log.write(f"\nBatch {batch_num} completed at {datetime.now()}\n")
        log.write(f"Downloaded: {downloaded}/{len(samples)}\n")
        log.write(f"Quantified: {quantified}/{len(samples)}\n")
        log.write(f"Cleaned up: {cleaned}/{len(samples)}\n")
    
    print(f"\n✅ Batch {batch_num} complete:")
    print(f"   Downloaded: {downloaded}/{len(samples)}")
    print(f"   Quantified: {quantified}/{len(samples)}")
    print(f"   Cleaned up: {cleaned}/{len(samples)}")
    
    return quantified

def main():
    """Main batch processing function."""
    base_dir = Path("output/amalgkit/pbarbatus")
    
    # Load remaining samples
    remaining_file = base_dir / "remaining_samples.txt"
    if not remaining_file.exists():
        print("❌ remaining_samples.txt not found!")
        print("   Run the setup script first")
        return 1
    
    with open(remaining_file) as f:
        remaining = [line.strip() for line in f if line.strip()]
    
    print("=" * 80)
    print("BATCH PROCESSING START")
    print("=" * 80)
    print(f"\nTotal samples to process: {len(remaining)}")
    print(f"Batch size: 10")
    print(f"Estimated time: {len(remaining) * 0.15:.1f} - {len(remaining) * 0.25:.1f} hours")
    
    # Create log file
    log_file = base_dir / f"batch_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Process in batches
    batch_size = 10
    total_quantified = 0
    
    for i in range(0, len(remaining), batch_size):
        batch_samples = remaining[i:i+batch_size]
        batch_num = (i // batch_size) + 1
        
        quantified = process_batch(batch_num, batch_samples, base_dir, log_file)
        total_quantified += quantified
        
        print(f"\n📊 Progress: {total_quantified}/{len(remaining)} samples quantified")
        
        # Brief pause between batches
        if i + batch_size < len(remaining):
            print("\n⏸️  5 second pause before next batch...")
            time.sleep(5)
    
    print("\n" + "=" * 80)
    print("✅ BATCH PROCESSING COMPLETE")
    print("=" * 80)
    print(f"\nTotal quantified: {total_quantified}/{len(remaining)}")
    print(f"Success rate: {total_quantified/len(remaining)*100:.1f}%")
    print(f"\nLog file: {log_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

