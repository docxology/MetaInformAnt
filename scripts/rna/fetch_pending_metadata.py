import argparse
import glob
import os
import subprocess
import sys

def main():
    configs = sorted(glob.glob("config/amalgkit/amalgkit_*.yaml"))
    to_run = []
    
    for cfg in configs:
        species = os.path.basename(cfg).replace("amalgkit_", "").replace(".yaml", "")
        if any(skip in species for skip in ["template", "test", "cross_species", "apis_mellifera_all"]):
            continue
            
        meta_selected = f"output/amalgkit/{species}/work/metadata/metadata_selected.tsv"
        meta_unselected = f"output/amalgkit/{species}/work/metadata/metadata.tsv"
        
        if not os.path.exists(meta_selected) and not os.path.exists(meta_unselected):
            to_run.append((species, cfg))
            
    if not to_run:
        print("All species already have metadata files!")
        return
        
    print(f"Found {len(to_run)} species pending metadata.")
    
    for species, cfg in to_run:
        print(f"--- Fetching metadata for {species} ---")
        cmd = [
            sys.executable,
            "scripts/rna/run_workflow.py",
            "--config", cfg,
            "--steps", "metadata", "select"
        ]
        
        try:
            print(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            print(f"✓ Metadata generated for {species}")
        except subprocess.CalledProcessError as e:
            print(f"✗ Error fetching metadata for {species}: {e}")

if __name__ == "__main__":
    main()
