#!/usr/bin/env python3
"""
Amalgkit Pipeline Status Monitor
================================

Scans the `blue/amalgkit` directory to report:
- Count of quantified samples per species
- Completion status of downstream steps (Merge, Curate, Sanity)
- Active failures or running processes

Usage:
  python3 scripts/rna/monitor_status.py
"""

import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import time

AMALGKIT_DIR = Path("blue/amalgkit")
SPECIES_ORDER = [
    # Small (<50)
    "anoplolepis_gracilipes", "acromyrmex_echinatior", "formica_exsecta", 
    "pbarbatus", "wasmannia_auropunctata", "temnothorax_nylanderi", 
    "nylanderia_fulva", "temnothorax_americanus", "temnothorax_curvispinosus",
    # Medium (<300)
    "cardiocondyla_obscurior", "linepithema_humile", "atta_cephalotes", 
    "ooceraea_biroi", "camponotus_floridanus",
    # Large (>300)
    "solenopsis_invicta", "apis_mellifera", "odontomachus_brunneus"
]

def get_species_status(species_name):
    sp_dir = AMALGKIT_DIR / species_name
    if not sp_dir.exists():
        return None
    
    # 1. Quant Count
    # Look in work/quant/SRR.../*_abundance.tsv OR quant/SRR.../*_abundance.tsv
    # Fast counting using glob
    quant_count = 0
    # Try finding all abundance files
    # This might be slow if naive, but with specific glob it's okay
    # Actually, we can just count directories in work/quant that have the file?
    # Or use the 'find' strategy via shell if python is slow? 
    # Python's rglob is okay for 3000 files if cached, but shell is faster.
    # Let's stick to Python for portability, but optimize.
    
    quant_dirs = []
    if (sp_dir / "work" / "quant").exists():
        quant_dirs.append(sp_dir / "work" / "quant")
    if (sp_dir / "quant").exists():
        quant_dirs.append(sp_dir / "quant")
        
    quant_files = []
    for qd in quant_dirs:
        quant_files.extend(list(qd.rglob("*_abundance.tsv")))
        quant_files.extend(list(qd.rglob("abundance.tsv")))
    
    # Unique SRRs
    srrs = set(f.parent.name for f in quant_files)
    quant_count = len(srrs)
    
    # 2. Step Completion Flags
    # merge, curate, sanity typically leave markers or specific output files
    # Merge: merged/merged_abundance.tsv
    # Curate: work/curate/Species_curated.tsv (or just curate_completion_flag?)
    # Sanity: work/sanity/Species_sanity.tsv (no, it produces plots/tables)
    
    has_merge = (sp_dir / "merged" / "merged_abundance.tsv").exists()
    
    # Curate check: look for completion flag OR output file
    has_curate = False
    curate_dir = sp_dir / "work" / "curate"
    # Try finding the flag in work/curate/Species/
    # Capitalize species name for directory if needed, or glob
    if list(curate_dir.glob("*/curate_completion_flag.txt")):
         has_curate = True
    elif (curate_dir / "curate_completion_flag.txt").exists():
        has_curate = True
    elif list(curate_dir.glob("*_curated.tsv")):
         has_curate = True
    elif list(curate_dir.glob("*/*_curated.tsv")): # work/curate/Species/Species_curated.tsv
         has_curate = True
         
    # Sanity check
    has_sanity = (sp_dir / "work" / "sanity_completion_flag.txt").exists()
    
    return {
        "name": species_name,
        "quant": quant_count,
        "merge": "✅" if has_merge else "no",
        "curate": "✅" if has_curate else "no",
        "sanity": "✅" if has_sanity else "no"
    }

def main():
    print(f"{'Species':<30} {'Quant':<8} {'Merge':<6} {'Curate':<6} {'Sanity':<6}")
    print("-" * 65)
    
    # Parallelize for speed
    with ThreadPoolExecutor(max_workers=8) as executor:
        results = executor.map(get_species_status, SPECIES_ORDER)
        
    for res in results:
        if res:
            # Add heuristics for status
            status_icon = ""
            if res['sanity'] == "✅":
                status_icon = "✅ DONE"
            elif res['curate'] == "✅":
                status_icon = "✅ Curated"
            elif res['merge'] == "✅":
                status_icon = "⚠️ Merged"
            elif res['quant'] > 0:
                status_icon = "🔄 In Progress"
                
            print(f"{res['name']:<30} {res['quant']:<8} {res['merge']:<6} {res['curate']:<6} {res['sanity']:<6}")
        else:
            # print(f"{species:<30} Not Found")
            pass

if __name__ == "__main__":
    main()
