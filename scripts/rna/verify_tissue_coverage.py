#!/usr/bin/env python3
"""
Verify Tissue Normalization Coverage
====================================
"""
import sys
import yaml
import pandas as pd
from pathlib import Path

CONFIG_DIR = Path("config/amalgkit")
AMALGKIT_DIR = Path("blue/amalgkit")

def main():
    print(f"{'Species':<30} {'Total':<10} {'Unmapped':<10} {'Status'}")
    print("-" * 70)
    
    # Load rules
    with open(CONFIG_DIR / "tissue_mapping.yaml") as f:
        mapping = yaml.safe_load(f)
    with open(CONFIG_DIR / "tissue_patches.yaml") as f:
        patches = yaml.safe_load(f)
        
    # Build inverted map
    inv_map = {}
    for canon, syns in mapping.items():
        if syns:
            for s in syns:
                inv_map[str(s).lower()] = canon

    # Scan species
    configs = sorted(list(CONFIG_DIR.glob("amalgkit_*.yaml")))
    for cfg in configs:
        name = cfg.stem
        if name in ["amalgkit_template", "amalgkit_test", "amalgkit_cross_species"]:
            continue
            
        species = name.replace("amalgkit_", "")
        
        # Check metadata
        meta_path = AMALGKIT_DIR / species / "work/metadata/metadata.tsv"
        if not meta_path.exists():
            meta_path = AMALGKIT_DIR / species / "work/metadata/metadata_selected.tsv"
            
        if not meta_path.exists():
            print(f"{species:<30} {'-':<10} {'-':<10} No Metadata (Pending Pipeline)")
            continue
            
        try:
            df = pd.read_csv(meta_path, sep="\t", low_memory=False)
        except:
            print(f"{species:<30} {'Error':<10} {'-':<10} Read Failed")
            continue
            
        col = "tissue"
        if col not in df.columns:
            print(f"{species:<30} {len(df):<10} {'All':<10} Missing 'tissue' column")
            continue
            
        unmapped = 0
        examples = []
        
        for _, row in df.iterrows():
            val = str(row.get(col, "")).strip()
            if not val or val.lower() == "nan":
                val = "missing"
                
            # Check patch
            srr = row.get("run") or row.get("run_accession")
            bp = row.get("bioproject")
            
            patched = False
            samples_patch = patches.get("samples") or {}
            bioprojects_patch = patches.get("bioprojects") or {}
            
            if srr in samples_patch: patched = True
            elif bp in bioprojects_patch: patched = True
            
            if patched: continue
            
            # Check mapping
            if val.lower() in inv_map: continue
            
            unmapped += 1
            if len(examples) < 3: examples.append(val)
            
        status = "OK"
        if unmapped > 0:
            status = f"MISSING: {', '.join(examples)}"
            
        print(f"{species:<30} {len(df):<10} {unmapped:<10} {status}")

if __name__ == "__main__":
    main()
