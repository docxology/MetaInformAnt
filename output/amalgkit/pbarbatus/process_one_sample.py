#!/usr/bin/env python3
"""Download and quantify a single sample."""

import sys
import subprocess
from pathlib import Path
import shutil

def process_sample(sample_id, base_dir):
    """Download, quantify, and cleanup one sample."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "fastq" / sample_id  # amalgkit puts in fastq/{sample}
    index_file = work_dir / "index" / "Pogonomyrmex_barbatus.idx"
    quant_out = work_dir / "quant" / sample_id
    
    print(f"\n{'=' * 80}")
    print(f"PROCESSING: {sample_id}")
    print(f"{'=' * 80}")
    
    # Step 1: Download FASTQ using sra-tools directly
    print(f"\n📥 Step 1: Downloading FASTQ...")
    
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    # Use prefetch + fasterq-dump (faster than amalgkit for single samples)
    try:
        # Prefetch SRA file
        print(f"   Downloading SRA...")
        subprocess.run(["prefetch", sample_id], check=True, cwd=str(fastq_dir))
        
        # Convert to FASTQ
        print(f"   Converting to FASTQ...")
        subprocess.run([
            "fasterq-dump",
            sample_id,
            "-O", str(fastq_dir),
            "-e", "6",  # 6 threads
            "-p"  # Show progress
        ], check=True)
        
        # Compress FASTQs
        print(f"   Compressing...")
        for fastq in fastq_dir.glob(f"{sample_id}*.fastq"):
            subprocess.run(["gzip", str(fastq)], check=True)
        
        print(f"   ✅ Download complete")
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        return False
    
    # Step 2: Quantify with kallisto
    print(f"\n📈 Step 2: Quantifying with kallisto...")
    
    fastq_1 = fastq_dir / f"{sample_id}_1.fastq.gz"
    fastq_2 = fastq_dir / f"{sample_id}_2.fastq.gz"
    
    if not (fastq_1.exists() and fastq_2.exists()):
        print(f"   ❌ FASTQs not found")
        return False
    
    try:
        result = subprocess.run([
            "kallisto", "quant",
            "-i", str(index_file),
            "-o", str(quant_out),
            "-t", "6",
            str(fastq_1),
            str(fastq_2)
        ], capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0 and (quant_out / "abundance.tsv").exists():
            print(f"   ✅ Quantification complete")
            
            # Show results
            with open(quant_out / "abundance.tsv") as f:
                lines = f.readlines()
                n_transcripts = len(lines) - 1
                # Count expressed
                expressed = sum(1 for line in lines[1:] if float(line.split('\t')[4]) > 0)
            print(f"      • {n_transcripts} transcripts")
            print(f"      • {expressed} expressed ({expressed/n_transcripts*100:.1f}%)")
        else:
            print(f"   ❌ Quantification failed")
            print(f"   {result.stderr[:200]}")
            return False
    except Exception as e:
        print(f"   ❌ Quantification error: {e}")
        return False
    
    # Step 3: Cleanup FASTQs
    print(f"\n🗑️  Step 3: Cleaning up FASTQs...")
    
    try:
        size_before = sum(f.stat().st_size for f in fastq_dir.rglob("*") if f.is_file())
        shutil.rmtree(fastq_dir)
        size_gb = size_before / (1024**3)
        print(f"   ✅ Freed {size_gb:.2f}GB")
    except Exception as e:
        print(f"   ⚠️  Cleanup warning: {e}")
    
    print(f"\n✅ {sample_id} COMPLETE\n")
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 process_one_sample.py SRR_ID")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    base_dir = Path("output/amalgkit/pbarbatus")
    
    success = process_sample(sample_id, base_dir)
    sys.exit(0 if success else 1)
