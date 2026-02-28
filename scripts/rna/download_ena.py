#!/usr/bin/env python3
import os
import sys
import json
import urllib.request
import urllib.parse
import subprocess
import argparse
from pathlib import Path

def get_ena_metadata(srr_id):
    """Query ENA API for fastq FTP links and MD5s."""
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srr_id}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=json"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as response:
            data = json.loads(response.read().decode('utf-8'))
            if data and isinstance(data, list) and len(data) > 0:
                return data[0]
    except Exception as e:
        print(f"[{srr_id}] ENA API Error: {e}", file=sys.stderr)
    return None

def download_file(url, out_path):
    """Download file using wget (robust for large files)."""
    # Ensure URL has protocol
    if not url.startswith('http') and not url.startswith('ftp'):
        url = 'https://' + url

    print(f"Downloading {url} to {out_path}...")
    cmd = ['wget', '-q', '--show-progress', '-c', '-t', '5', '-T', '300', '-O', str(out_path), url]
    
    # Run wget and stream output
    process = subprocess.Popen(cmd)
    process.communicate()
    
    return process.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0

def check_md5(file_path, expected_md5):
    """Verify md5sum of a file."""
    try:
        # Use system md5sum for speed
        result = subprocess.run(['md5sum', str(file_path)], capture_output=True, text=True, check=True)
        actual_md5 = result.stdout.split()[0]
        return actual_md5 == expected_md5
    except subprocess.CalledProcessError:
        return False

def main():
    parser = argparse.ArgumentParser(description="Download fastq files directly from ENA")
    parser.add_argument("srr_id", help="SRR/ERR/DRR Accession (e.g., SRR37117694)")
    parser.add_argument("out_dir", help="Output directory for fastq files")
    args = parser.parse_args()

    srr_id = args.srr_id
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{srr_id}] Querying ENA API...")
    meta = get_ena_metadata(srr_id)

    if not meta or not meta.get('fastq_ftp'):
        print(f"[{srr_id}] Not found on ENA or no fastq FTP links available.")
        sys.exit(1)

    ftps = meta['fastq_ftp'].split(';')
    md5s = meta['fastq_md5'].split(';')
    
    if len(ftps) != len(md5s):
        print(f"[{srr_id}] Mismatch in FTP links and MD5 counts.")
        sys.exit(1)

    all_success = True
    
    # ENA returns files like: SRR..._1.fastq.gz, SRR..._2.fastq.gz
    for i, (ftp_url, expected_md5) in enumerate(zip(ftps, md5s)):
        # Determine proper amalgkit expected name
        if len(ftps) == 1:
            filename = f"{srr_id}.fastq.gz"
        elif len(ftps) == 2:
            filename = f"{srr_id}_{i+1}.fastq.gz"
        else:
            # Fallback for complex layouts (e.g. 3 files: normal, _1, _2)
            orig_name = ftp_url.split('/')[-1]
            filename = orig_name
            
        out_path = out_dir / filename
        
        # Check if already downloaded and valid
        if out_path.exists():
            print(f"[{srr_id}] File {filename} exists. Checking MD5...")
            if check_md5(out_path, expected_md5):
                print(f"[{srr_id}] MD5 matches! Skipping download.")
                continue
            else:
                print(f"[{srr_id}] MD5 mismatch. Re-downloading...")
                out_path.unlink() # Delete corrupted file

        # Download
        success = download_file(ftp_url, out_path)
        if not success:
            print(f"[{srr_id}] Download failed for {ftp_url}")
            all_success = False
            continue
            
        # Verify
        print(f"[{srr_id}] Verifying MD5 for {filename}...")
        if check_md5(out_path, expected_md5):
            print(f"[{srr_id}] MD5 verification successful.")
        else:
            print(f"[{srr_id}] MD5 verification failed after download.")
            out_path.unlink() # Delete corrupted file
            all_success = False

    if all_success:
        print(f"[{srr_id}] Successfully downloaded all fastq files from ENA.")
        sys.exit(0)
    else:
        print(f"[{srr_id}] Failed to complete ENA download.")
        sys.exit(1)

if __name__ == "__main__":
    main()
