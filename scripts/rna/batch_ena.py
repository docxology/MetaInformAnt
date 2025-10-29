#!/usr/bin/env python3
"""
Fast ENA-based parallel downloader for P. barbatus samples.

Uses ENA (European Nucleotide Archive) which provides:
- Direct FASTQ.gz downloads (no SRA conversion!)
- Often 5-10X faster than NCBI
- Complete mirror of SRA data

This replaces the slow NCBI prefetch + fasterq-dump workflow.
"""

import sys
import subprocess
from pathlib import Path
import shutil
import time
from datetime import datetime
import json
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
MAX_CONCURRENT_DOWNLOADS = 5  # Parallel sample downloads (user preference)
MAX_CONCURRENT_QUANTS = 3     # Parallel quantifications (user preference)
KALLISTO_THREADS = 3          # Threads per kallisto process
DOWNLOAD_TIMEOUT = 3600       # 1 hour per FASTQ file (more realistic for ENA)
ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

def get_ena_fastq_urls(sample_id):
    """Query ENA API for direct FASTQ URLs."""
    url = f"{ENA_API_URL}?accession={sample_id}&result=read_run&fields=fastq_ftp&format=tsv"
    
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
        
        lines = content.strip().split('\n')
        if len(lines) > 1:
            # Parse TSV
            header = lines[0].split('\t')
            data = lines[1].split('\t')
            
            if 'fastq_ftp' in header:
                idx = header.index('fastq_ftp')
                if len(data) > idx and data[idx]:
                    ftp_urls = data[idx].split(';')
                    # Convert FTP to HTTP (usually faster)
                    http_urls = [f"http://{url}" for url in ftp_urls]
                    return http_urls
        
        return None
        
    except Exception as e:
        print(f"   ⚠️  ENA API error for {sample_id}: {e}")
        return None

def download_fastq_file(url, dest_path, sample_id):
    """Download a single FASTQ file with progress."""
    filename = url.split('/')[-1]
    
    try:
        # Use wget for better progress and resume capability
        result = subprocess.run([
            "wget",
            "-c",  # Continue partial downloads
            "-t", "3",  # 3 retries
            "-T", "30",  # 30s timeout per connection
            "-O", str(dest_path),
            url
        ], capture_output=True, text=True, timeout=DOWNLOAD_TIMEOUT)
        
        if result.returncode == 0 and dest_path.exists():
            size_mb = dest_path.stat().st_size / (1024**2)
            return True, size_mb
        else:
            return False, f"wget failed: {result.stderr[:200]}"
            
    except subprocess.TimeoutExpired:
        return False, "Download timeout"
    except Exception as e:
        return False, str(e)

def download_sample_ena(sample_id, base_dir):
    """Download FASTQ for one sample from ENA."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "fastq" / sample_id
    
    log_entry = {
        "sample_id": sample_id,
        "method": "ena",
        "start_time": datetime.now().isoformat()
    }
    
    print(f"📥 [{sample_id}] Starting ENA download...")
    
    # Create directory
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Step 1: Get ENA URLs
        print(f"   [{sample_id}] Querying ENA...")
        urls = get_ena_fastq_urls(sample_id)
        
        if not urls:
            log_entry["status"] = "failed"
            log_entry["error"] = "No ENA URLs found"
            print(f"   ❌ [{sample_id}] Not available on ENA")
            return False, log_entry
        
        log_entry["urls"] = urls
        log_entry["file_count"] = len(urls)
        print(f"   [{sample_id}] Found {len(urls)} FASTQ files on ENA")
        
        # Step 2: Download files
        downloaded_files = []
        total_size_mb = 0
        
        for i, url in enumerate(urls, 1):
            filename = url.split('/')[-1]
            dest_path = fastq_dir / filename
            
            print(f"   [{sample_id}] Downloading {filename} ({i}/{len(urls)})...")
            success, result = download_fastq_file(url, dest_path, sample_id)
            
            if success:
                size_mb = result
                total_size_mb += size_mb
                downloaded_files.append(filename)
                print(f"   ✅ [{sample_id}] {filename}: {size_mb:.1f}MB")
            else:
                log_entry["status"] = "failed"
                log_entry["error"] = f"Download failed for {filename}: {result}"
                print(f"   ❌ [{sample_id}] {filename} failed: {result}")
                return False, log_entry
        
        # Verify all files
        fastq_1 = fastq_dir / f"{sample_id}_1.fastq.gz"
        fastq_2 = fastq_dir / f"{sample_id}_2.fastq.gz"
        
        if fastq_1.exists() and fastq_2.exists():
            log_entry["status"] = "success"
            log_entry["downloaded_files"] = downloaded_files
            log_entry["total_size_mb"] = total_size_mb
            log_entry["end_time"] = datetime.now().isoformat()
            print(f"✅ [{sample_id}] ENA download complete: {total_size_mb:.1f}MB total")
            return True, log_entry
        else:
            log_entry["status"] = "failed"
            log_entry["error"] = "Expected files not found after download"
            return False, log_entry
            
    except Exception as e:
        log_entry["status"] = "failed"
        log_entry["error"] = str(e)
        log_entry["end_time"] = datetime.now().isoformat()
        print(f"❌ [{sample_id}] ENA download error: {e}")
        return False, log_entry

def quantify_sample(sample_id, base_dir):
    """Quantify one sample with kallisto."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "fastq" / sample_id
    index_file = work_dir / "index" / "Pogonomyrmex_barbatus.idx"
    quant_out = work_dir / "quant" / sample_id
    
    log_entry = {
        "sample_id": sample_id,
        "stage": "quantify",
        "start_time": datetime.now().isoformat()
    }
    
    print(f"📈 [{sample_id}] Starting quantification...")
    
    fastq_1 = fastq_dir / f"{sample_id}_1.fastq.gz"
    fastq_2 = fastq_dir / f"{sample_id}_2.fastq.gz"
    
    if not (fastq_1.exists() and fastq_2.exists()):
        log_entry["status"] = "failed"
        log_entry["error"] = "FASTQs not found"
        return False, log_entry
    
    try:
        result = subprocess.run([
            "kallisto", "quant",
            "-i", str(index_file),
            "-o", str(quant_out),
            "-t", str(KALLISTO_THREADS),  # Configurable threads per sample
            str(fastq_1),
            str(fastq_2)
        ], capture_output=True, text=True, timeout=1200)  # 20 min timeout
        
        if result.returncode == 0 and (quant_out / "abundance.tsv").exists():
            # Parse results
            with open(quant_out / "abundance.tsv") as f:
                lines = f.readlines()
                n_transcripts = len(lines) - 1
                expressed = sum(1 for line in lines[1:] if float(line.split('\t')[4]) > 0)
            
            log_entry["status"] = "success"
            log_entry["transcripts"] = n_transcripts
            log_entry["expressed"] = expressed
            log_entry["expression_rate"] = f"{expressed/n_transcripts*100:.1f}%"
            log_entry["end_time"] = datetime.now().isoformat()
            
            print(f"✅ [{sample_id}] Quantification complete: {expressed}/{n_transcripts} expressed")
            return True, log_entry
        else:
            log_entry["status"] = "failed"
            log_entry["error"] = result.stderr[:200]
            return False, log_entry
            
    except Exception as e:
        log_entry["status"] = "failed"
        log_entry["error"] = str(e)
        return False, log_entry

def cleanup_sample(sample_id, base_dir):
    """Delete FASTQ files after successful quantification."""
    work_dir = base_dir / "work"
    fastq_dir = work_dir / "fastq" / sample_id
    
    if not fastq_dir.exists():
        return True
    
    try:
        size_before = sum(f.stat().st_size for f in fastq_dir.rglob("*") if f.is_file())
        shutil.rmtree(fastq_dir)
        size_gb = size_before / (1024**3)
        print(f"🗑️  [{sample_id}] Freed {size_gb:.2f}GB")
        return True
    except Exception as e:
        print(f"⚠️  [{sample_id}] Cleanup warning: {e}")
        return False

def process_sample_pipeline(sample_id, base_dir):
    """Complete pipeline: ENA download -> quantify -> cleanup."""
    quant_dir = base_dir / "work" / "quant" / sample_id
    
    # Skip if already done
    if (quant_dir / "abundance.tsv").exists():
        print(f"⏭️  [{sample_id}] Already quantified, skipping")
        return True, {"sample_id": sample_id, "status": "skipped"}
    
    start_time = time.time()
    
    # Download from ENA
    download_success, download_log = download_sample_ena(sample_id, base_dir)
    if not download_success:
        return False, download_log
    
    # Quantify
    quant_success, quant_log = quantify_sample(sample_id, base_dir)
    if not quant_success:
        cleanup_sample(sample_id, base_dir)
        return False, quant_log
    
    # Cleanup
    cleanup_sample(sample_id, base_dir)
    
    elapsed = time.time() - start_time
    combined_log = {
        "sample_id": sample_id,
        "status": "success",
        "method": "ena",
        "elapsed_minutes": f"{elapsed/60:.1f}",
        "download": download_log,
        "quantify": quant_log
    }
    
    return True, combined_log

def main():
    """Process all samples using ENA downloads."""
    base_dir = Path("output/amalgkit/pbarbatus")
    
    # Load remaining samples
    remaining_file = base_dir / "remaining_samples.txt"
    if not remaining_file.exists():
        print("❌ remaining_samples.txt not found")
        return 1
    
    with open(remaining_file) as f:
        remaining = [line.strip() for line in f if line.strip()]
    
    # Filter out already completed
    quant_dir = base_dir / "work" / "quant"
    completed = {d.name for d in quant_dir.iterdir() 
                 if d.is_dir() and (d / "abundance.tsv").exists()}
    to_process = [s for s in remaining if s not in completed]
    
    print("=" * 80)
    print("ENA PARALLEL BATCH PROCESSING")
    print("=" * 80)
    print(f"\nTotal samples: {len(remaining)}")
    print(f"Already completed: {len(completed)}")
    print(f"To process: {len(to_process)}")
    print(f"\nConfiguration:")
    print(f"  • Method: ENA direct FASTQ download")
    print(f"  • Concurrent downloads: {MAX_CONCURRENT_DOWNLOADS}")
    print(f"  • Concurrent quantifications: {MAX_CONCURRENT_QUANTS}")
    print(f"  • Kallisto threads: {KALLISTO_THREADS} per sample")
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Log file
    log_file = base_dir / f"batch_ena_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    log_data = {
        "method": "ena",
        "start_time": datetime.now().isoformat(),
        "total_samples": len(to_process),
        "max_concurrent_downloads": MAX_CONCURRENT_DOWNLOADS,
        "max_concurrent_quants": MAX_CONCURRENT_QUANTS,
        "kallisto_threads": KALLISTO_THREADS,
        "samples": []
    }
    
    # Process with thread pool
    success_count = 0
    failed_count = 0
    
    print(f"\n🚀 Starting ENA parallel processing...")
    print("=" * 80)
    
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_DOWNLOADS) as executor:
        # Submit all samples
        future_to_sample = {
            executor.submit(process_sample_pipeline, sample_id, base_dir): sample_id
            for sample_id in to_process
        }
        
        # Process as they complete
        for i, future in enumerate(as_completed(future_to_sample), 1):
            sample_id = future_to_sample[future]
            
            try:
                success, sample_log = future.result()
                log_data["samples"].append(sample_log)
                
                if success:
                    success_count += 1
                else:
                    failed_count += 1
                
                # Save log after each completion
                with open(log_file, 'w') as f:
                    json.dump(log_data, f, indent=2)
                
                # Progress update
                print(f"\n📊 Progress: {i}/{len(to_process)} | Success: {success_count} | Failed: {failed_count}")
                
            except Exception as e:
                print(f"❌ [{sample_id}] Unexpected error: {e}")
                failed_count += 1
    
    # Final summary
    log_data["end_time"] = datetime.now().isoformat()
    log_data["success_count"] = success_count
    log_data["failed_count"] = failed_count
    
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    
    print("\n" + "=" * 80)
    print("✅ ENA PARALLEL PROCESSING COMPLETE")
    print("=" * 80)
    print(f"\nTotal processed: {len(to_process)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {failed_count}")
    print(f"Success rate: {success_count/len(to_process)*100:.1f}%")
    print(f"\nLog file: {log_file}")
    
    return 0 if failed_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())

