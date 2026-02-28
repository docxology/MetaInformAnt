#!/usr/bin/env python3
"""
Fix Harpegnathos Index Complexity.

Filters the Harpegnathos saltator transcriptome to remove:
1. Low complexity/repetitive non-coding RNAs (XR_, NR_ prefixes)
2. Short fragments (< 200bp)
3. Ribosomal RNAs (specifically identified in headers)

Then rebuilds the kallisto index.
"""

import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path

# Paths
GENOME_DIR = Path("output/amalgkit/shared/genome/Harpegnathos_saltator")
INPUT_FASTA = GENOME_DIR / "GCF_003227715.2_Hsal_v8.6_rna_from_genomic.fna.gz"
OUTPUT_FASTA = GENOME_DIR / "Harpegnathos_saltator_filtered.fna"
OUTPUT_INDEX = GENOME_DIR / "index/Harpegnathos_saltator_filtered.idx"

# Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

def filter_fasta(input_path: Path, output_path: Path):
    """Filter FASTA file."""
    logger.info(f"Filtering {input_path} -> {output_path}")
    
    kept = 0
    removed_xr = 0
    removed_short = 0
    removed_rrna = 0
    
    try:
        with gzip.open(input_path, "rt") as f_in, open(output_path, "w") as f_out:
            header = None
            sequence = []
            
            def process_entry(h, s):
                nonlocal kept, removed_xr, removed_short, removed_rrna
                if not h:
                    return
                
                s_str = "".join(s)
                
                # Check 1: Exclude XR_/NR_ (ncRNA)
                if "XR_" in h or "NR_" in h:
                    removed_xr += 1
                    return
                
                # Check 2: Exclude rRNA explicit
                if "ribosomal RNA" in h or "rRNA" in h:
                    removed_rrna += 1
                    return
                
                # Check 3: Length < 200bp
                if len(s_str) < 200:
                    removed_short += 1
                    return
                
                # Keep
                f_out.write(f"{h}\n{s_str}\n")
                kept += 1

            for line in f_in:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        process_entry(header, sequence)
                    header = line
                    sequence = []
                else:
                    sequence.append(line)
            
            # Process last entry
            if header:
                process_entry(header, sequence)
                
    except Exception as e:
        logger.error(f"Failed to filter FASTA: {e}")
        sys.exit(1)
        
    logger.info("Filtering complete:")
    logger.info(f"  Kept: {kept}")
    logger.info(f"  Removed (XR/NR): {removed_xr}")
    logger.info(f"  Removed (rRNA explicit): {removed_rrna}")
    logger.info(f"  Removed (Short <200bp): {removed_short}")
    logger.info(f"  Total Removed: {removed_xr + removed_short + removed_rrna}")

def build_index(fasta_path: Path, index_path: Path):
    """Build Kallisto index."""
    logger.info(f"Building Kallisto index: {index_path}")
    
    cmd = [
        "kallisto", "index",
        "-i", str(index_path),
        str(fasta_path)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        logger.info("Index built successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Kallisto index failed: {e}")
        sys.exit(1)

def main():
    if not INPUT_FASTA.exists():
        logger.error(f"Input FASTA not found: {INPUT_FASTA}")
        sys.exit(1)
        
    filter_fasta(INPUT_FASTA, OUTPUT_FASTA)
    
    # Ensure index dir exists
    OUTPUT_INDEX.parent.mkdir(parents=True, exist_ok=True)
    
    build_index(OUTPUT_FASTA, OUTPUT_INDEX)
    
    logger.info("Done.")
    logger.info(f"New index: {OUTPUT_INDEX}")
    logger.info("Update your amalgkit config to point to this index_dir/index.")

if __name__ == "__main__":
    main()
