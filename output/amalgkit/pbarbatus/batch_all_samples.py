#!/usr/bin/env python3
"""
Batch processor for all remaining P. barbatus samples.

Processes samples one at a time to minimize disk space usage:
- Download FASTQ
- Quantify with kallisto
- Delete FASTQ immediately
- Repeat for next sample

Usage:
    cd /Users/4d/Documents/GitHub/metainformant
    python3 output/amalgkit/pbarbatus/batch_all_samples.py
"""

import sys
import subprocess
from pathlib import Path
import shutil
import time
from datetime import datetime
import json
import signal
import psutil

def kill_blocking_processes(sample_id, log_entry):
    """Kill any blocking prefetch/fasterq/kallisto processes."""
    killed = []
    
    try:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                cmdline = ' '.join(proc.info['cmdline'] or [])
                
                # Check for blocking processes related to this sample or general
                if (proc.info['name'] in ['prefetch', 'fasterq-dump', 'kallisto'] and
                    (sample_id in cmdline or 'SRR' in cmdline)):
                    
                    print(f"   🔪 Killing blocking process: {proc.info['name']} (PID {proc.info['pid']})")
                    proc.kill()
                    proc.wait(timeout=5)
                    killed.append({
                        'pid': proc.info['pid'],
                        'name': proc.info['name'],
                        'cmdline': cmdline[:100]
                    })
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.TimeoutExpired):
                continue
    except Exception as e:
        print(f"   ⚠️  Warning during process cleanup: {e}")
    
    if killed:
        log_entry['killed_processes'] = killed
        print(f"   ✅ Killed {len(killed)} blocking process(es)")
    
    return killed

def process_sample(sample_id, base_dir, log_data):
    """Download, quantify, and cleanup one sample."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "fastq" / sample_id
    index_file = work_dir / "index" / "Pogonomyrmex_barbatus.idx"
    quant_out = work_dir / "quant" / sample_id
    
    start_time = time.time()
    sample_log = {
        "sample_id": sample_id,
        "start_time": datetime.now().isoformat(),
        "status": "started"
    }
    
    print(f"\n{'=' * 80}")
    print(f"PROCESSING: {sample_id}")
    print(f"{'=' * 80}")
    
    # Kill any blocking processes first
    print(f"\n🔍 Checking for blocking processes...")
    killed = kill_blocking_processes(sample_id, sample_log)
    if killed:
        time.sleep(2)  # Give processes time to fully terminate
    
    # Check if already quantified
    if (quant_out / "abundance.tsv").exists():
        print(f"   ⏭️  Already quantified, skipping")
        sample_log["status"] = "skipped"
        sample_log["reason"] = "already_quantified"
        return True, sample_log
    
    # Step 1: Download FASTQ
    print(f"\n📥 Step 1: Downloading FASTQ...")
    
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Prefetch
        print(f"   Downloading SRA...")
        subprocess.run(
            ["prefetch", sample_id],
            check=True,
            cwd=str(fastq_dir),
            timeout=600,
            capture_output=True
        )
        
        # Convert to FASTQ
        print(f"   Converting to FASTQ...")
        subprocess.run([
            "fasterq-dump",
            sample_id,
            "-O", str(fastq_dir),
            "-e", "6",
            "-p"
        ], check=True, timeout=1200)
        
        # Compress
        print(f"   Compressing...")
        for fastq in fastq_dir.glob(f"{sample_id}*.fastq"):
            subprocess.run(["gzip", str(fastq)], check=True, timeout=300)
        
        sample_log["download"] = "success"
        print(f"   ✅ Download complete")
    except subprocess.TimeoutExpired:
        print(f"   ❌ Download timeout")
        sample_log["download"] = "timeout"
        sample_log["status"] = "failed"
        return False, sample_log
    except Exception as e:
        print(f"   ❌ Download failed: {e}")
        sample_log["download"] = f"failed: {e}"
        sample_log["status"] = "failed"
        return False, sample_log
    
    # Step 2: Quantify
    print(f"\n📈 Step 2: Quantifying...")
    
    fastq_1 = fastq_dir / f"{sample_id}_1.fastq.gz"
    fastq_2 = fastq_dir / f"{sample_id}_2.fastq.gz"
    
    if not (fastq_1.exists() and fastq_2.exists()):
        print(f"   ❌ FASTQs not found")
        sample_log["quant"] = "fastq_not_found"
        sample_log["status"] = "failed"
        return False, sample_log
    
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
            # Parse results
            with open(quant_out / "abundance.tsv") as f:
                lines = f.readlines()
                n_transcripts = len(lines) - 1
                expressed = sum(1 for line in lines[1:] if float(line.split('\t')[4]) > 0)
            
            sample_log["quant"] = "success"
            sample_log["transcripts"] = n_transcripts
            sample_log["expressed"] = expressed
            sample_log["expression_rate"] = f"{expressed/n_transcripts*100:.1f}%"
            
            print(f"   ✅ Quantification complete")
            print(f"      • {n_transcripts} transcripts")
            print(f"      • {expressed} expressed ({expressed/n_transcripts*100:.1f}%)")
            
            # Parse mapping rate from stderr
            if "reads pseudoaligned" in result.stderr:
                for line in result.stderr.split('\n'):
                    if "reads pseudoaligned" in line:
                        sample_log["mapping_info"] = line.strip()
        else:
            print(f"   ❌ Quantification failed")
            sample_log["quant"] = f"failed: {result.stderr[:100]}"
            sample_log["status"] = "failed"
            return False, sample_log
    except subprocess.TimeoutExpired:
        print(f"   ❌ Quantification timeout")
        sample_log["quant"] = "timeout"
        sample_log["status"] = "failed"
        return False, sample_log
    except Exception as e:
        print(f"   ❌ Quantification error: {e}")
        sample_log["quant"] = f"error: {e}"
        sample_log["status"] = "failed"
        return False, sample_log
    
    # Step 3: Cleanup
    print(f"\n🗑️  Step 3: Cleaning up...")
    
    try:
        size_before = sum(f.stat().st_size for f in fastq_dir.rglob("*") if f.is_file())
        shutil.rmtree(fastq_dir)
        size_gb = size_before / (1024**3)
        sample_log["cleanup"] = f"{size_gb:.2f}GB freed"
        print(f"   ✅ Freed {size_gb:.2f}GB")
    except Exception as e:
        print(f"   ⚠️  Cleanup warning: {e}")
        sample_log["cleanup"] = f"warning: {e}"
    
    elapsed = time.time() - start_time
    sample_log["elapsed_minutes"] = f"{elapsed/60:.1f}"
    sample_log["status"] = "success"
    sample_log["end_time"] = datetime.now().isoformat()
    
    print(f"\n✅ {sample_id} COMPLETE ({elapsed/60:.1f} minutes)\n")
    return True, sample_log

def main():
    """Process all remaining samples."""
    base_dir = Path("output/amalgkit/pbarbatus")
    
    # Load remaining samples
    remaining_file = base_dir / "remaining_samples.txt"
    if not remaining_file.exists():
        print("❌ remaining_samples.txt not found")
        return 1
    
    with open(remaining_file) as f:
        remaining = [line.strip() for line in f if line.strip()]
    
    print("=" * 80)
    print("BATCH PROCESSING: ALL REMAINING SAMPLES")
    print("=" * 80)
    print(f"\nTotal samples: {len(remaining)}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Estimated total time: {len(remaining) * 0.15:.1f} - {len(remaining) * 0.25:.1f} hours")
    
    # Log file
    log_file = base_dir / f"batch_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    log_data = {
        "start_time": datetime.now().isoformat(),
        "total_samples": len(remaining),
        "samples": []
    }
    
    # Process each sample
    success_count = 0
    failed_count = 0
    
    for i, sample_id in enumerate(remaining, 1):
        print(f"\n[{i}/{len(remaining)}] {sample_id}")
        
        success, sample_log = process_sample(sample_id, base_dir, log_data)
        log_data["samples"].append(sample_log)
        
        if success:
            success_count += 1
        else:
            failed_count += 1
        
        # Save log after each sample
        with open(log_file, 'w') as f:
            json.dump(log_data, f, indent=2)
        
        # Progress update
        print(f"\n📊 Progress: {success_count}/{len(remaining)} successful, {failed_count} failed")
        
        # Brief pause between samples
        if i < len(remaining):
            time.sleep(2)
    
    # Final summary
    log_data["end_time"] = datetime.now().isoformat()
    log_data["success_count"] = success_count
    log_data["failed_count"] = failed_count
    
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    
    print("\n" + "=" * 80)
    print("✅ BATCH PROCESSING COMPLETE")
    print("=" * 80)
    print(f"\nTotal processed: {len(remaining)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {failed_count}")
    print(f"Success rate: {success_count/len(remaining)*100:.1f}%")
    print(f"\nLog file: {log_file}")
    
    return 0 if failed_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())

