#!/usr/bin/env python3
"""
Count samples in Amalgkit stages for each species.
"""

import os
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path("/Volumes/blue/data/amalgkit")

def count_files(directory, pattern="*"):
    if not directory.exists():
        return 0
    return len(list(directory.glob(pattern)))

def main():
    if not DATA_DIR.exists():
        print(f"Data directory {DATA_DIR} not found.")
        return

    species_dirs = [d for d in DATA_DIR.iterdir() if d.is_dir() and d.name != "shared"]
    
    print(f"{'Species':<30} | {'FastQ':<10} | {'GetFastQ':<10} | {'Quant':<10} | {'Work/Quant':<10}")
    print("-" * 80)

    for species_dir in sorted(species_dirs, key=lambda x: x.name):
        name = species_dir.name
        
        # Paths (adjusting based on observed structure)
        # Structure seems to be: species_dir/fastq, species_dir/quant OR species_dir/work/quant
        
        fastq_dir = species_dir / "fastq"
        # getfastq might be inside fastq/getfastq or work/getfastq
        getfastq_dir_1 = species_dir / "fastq" / "getfastq"
        getfastq_dir_2 = species_dir / "work" / "getfastq"
        
        quant_dir = species_dir / "quant"
        work_quant_dir = species_dir / "work" / "quant"
        
        # Counting (assuming one file/dir per sample or distinct file types)
        # FastQ usually has *.fastq.gz
        n_fastq = count_files(fastq_dir, "*.fastq.gz")
        
        # GetFastQ usually has subdirectories per sample or files
        n_getfastq = 0
        if getfastq_dir_1.exists():
            n_getfastq += len([d for d in getfastq_dir_1.iterdir() if d.is_dir()])
        if getfastq_dir_2.exists():
            n_getfastq += len([d for d in getfastq_dir_2.iterdir() if d.is_dir()])

        # Quant usually has subdirectories per sample or abundant tsvs
        n_quant = 0
        if quant_dir.exists():
             n_quant += len([d for d in quant_dir.iterdir() if d.is_dir()])
        
        n_work_quant = 0
        if work_quant_dir.exists():
            n_work_quant += len([d for d in work_quant_dir.iterdir() if d.is_dir()])

        print(f"{name:<30} | {n_fastq:<10} | {n_getfastq:<10} | {n_quant:<10} | {n_work_quant:<10}")

if __name__ == "__main__":
    main()
