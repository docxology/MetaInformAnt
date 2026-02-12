#!/usr/bin/env python3
"""
Batch Genome Download and Indexing Script (Manual Download)
===========================================================
Iterates through all verified `amalgkit_*.yaml` configs in `config/amalgkit/`,
excludes templates/tests, and performs:
1. Downloads genome assets (FASTA/GTF/GFF) using `curl` from the FTP URL in config.
2. Runs `kallisto index` to build the transcriptome index.

This replaces the failed `amalgkit genome` approach.
"""

import os
import glob
import subprocess
import yaml
import sys
from concurrent.futures import ThreadPoolExecutor

CONFIG_DIR = "config/amalgkit"
EXCLUDE = [
    "amalgkit_template.yaml",
    "amalgkit_test.yaml",
    "amalgkit_cross_species.yaml",
    "amalgkit_apis_mellifera_all.yaml",  # Already processing/processed separately
    "tissue_mapping.yaml",
    "tissue_patches.yaml"
]

def get_species_configs():
    configs = []
    for fpath in sorted(glob.glob(os.path.join(CONFIG_DIR, "amalgkit_*.yaml"))):
        bn = os.path.basename(fpath)
        if bn in EXCLUDE:
            continue
        configs.append(fpath)
    return configs

def process_species(config_path):
    bn = os.path.basename(config_path)
    print(f"🐜 Processing {bn}...")
    
    try:
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
            
        genome = cfg.get('genome', {})
        ftp_url = genome.get('ftp_url', '')
        dest_dir = genome.get('dest_dir', '')
        files = genome.get('files', {})
        
        if not ftp_url or not dest_dir:
            print(f"   ❌ [Config] Missing ftp_url or dest_dir in {bn}")
            return False
            
        # Ensure dest_dir exists
        os.makedirs(dest_dir, exist_ok=True)
        
        # 1. Download Files
        print(f"   [Download] Downloading assets for {bn}...")
        download_success = True
        
        # Files to download: usually mapped in 'files' dict
        # keys: genomic_fasta, transcriptome_fasta, cds_fasta, protein_fasta, annotation_gff
        
        for file_type, filename in files.items():
            if not filename:
                continue
                
            local_path = os.path.join(dest_dir, filename)
            remote_url = os.path.join(ftp_url, filename)
            
            if os.path.exists(local_path) and os.path.getsize(local_path) > 1000:
                # print(f"      - {filename} exists, skipping.")
                pass
            else:
                print(f"      ⬇️ Downloading {filename}...")
                cmd_curl = ["curl", "-s", "-f", "-o", local_path, remote_url] # -f fails silently on error
                res = subprocess.run(cmd_curl)
                if res.returncode != 0:
                    print(f"      ❌ Failed to download {remote_url}")
                    download_success = False
        
        if not download_success:
            print(f"   ❌ [Download] Failed incomplete download for {bn}")
            return False
            
        print(f"   ✅ [Download] Complete for {bn}")

        # 2. Build Index
        rna_fasta = files.get('transcriptome_fasta')
        if not rna_fasta:
            print(f"   ⚠️  [Index] No transcriptome_fasta for {bn} (GenBank only?)")
            return True
            
        fasta_path = os.path.join(dest_dir, rna_fasta)
        
        # Handle compressed FASTA
        # Kallisto usually handles .gz input since v0.44.0? Yes.
        
        index_dir = os.path.join(dest_dir, "index")
        os.makedirs(index_dir, exist_ok=True)
        
        index_name = rna_fasta + ".idx"
        index_path = os.path.join(index_dir, index_name)
        
        if os.path.exists(index_path):
            print(f"   ✅ [Index] Exists: {index_name}")
            return True # Skip if exists to save time
        
        print(f"   [Index] Building kallisto index for {bn}...")
        # Start build
        cmd_index = ["kallisto", "index", "-i", index_path, fasta_path]
        
        subprocess.run(cmd_index, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        print(f"   ✅ [Index] Built {index_path}")
        return True
        
    except Exception as e:
        print(f"   ❌ [Process] Failed for {bn}: {e}")
        return False

def main():
    configs = get_species_configs()
    print(f"Found {len(configs)} verified species configs.")
    
    # 4 workers for download/index
    results = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_species, c): c for c in configs}
        for future in futures:
            c = futures[future]
            try:
                success = future.result()
                results[c] = success
            except Exception as e:
                print(f"   ❌ Unhandled exception for {c}: {e}")
                results[c] = False

    print("\nSummary:")
    success_count = sum(1 for s in results.values() if s)
    print(f"Processed {len(results)} species. Success: {success_count}/{len(results)}")
    
    if success_count < len(results):
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()
