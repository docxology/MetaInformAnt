import sys
import yaml
import pandas as pd
from pathlib import Path

CONFIG_DIR = Path("config/amalgkit")
AMALGKIT_DIR = Path("output/amalgkit")

def main():
    print(f"{'Species':<30} {'BioProject':<20} {'Samples':<10} {'SRA Examples'}")
    print("-" * 80)
    
    with open(CONFIG_DIR / "tissue_patches.yaml") as f:
        patches = yaml.safe_load(f)
        
    configs = sorted(list(CONFIG_DIR.glob("amalgkit_*.yaml")))
    for cfg in configs:
        name = cfg.stem
        if name in ["amalgkit_template", "amalgkit_test", "amalgkit_cross_species"]:
            continue
            
        species = name.replace("amalgkit_", "")
        
        meta_path = AMALGKIT_DIR / species / "work/metadata/metadata.tsv"
        if not meta_path.exists():
            meta_path = AMALGKIT_DIR / species / "work/metadata/metadata_selected.tsv"
            
        if not meta_path.exists():
            continue
            
        try:
            df = pd.read_csv(meta_path, sep="\t", low_memory=False)
        except:
            continue
            
        col = "tissue"
        missing_bps = {}
        
        for _, row in df.iterrows():
            val = str(row.get(col, "")).strip()
            # If the raw metadata is empty/nan/missing
            is_missing = not val or val.lower() in ("nan", "missing", "not applicable", "n/a", "na", "not_applicable", "unknown", "unknow")
            
            srr = str(row.get("run") or row.get("run_accession"))
            bp = str(row.get("bioproject"))
            
            patched = False
            if srr in (patches.get("samples") or {}): patched = True
            if bp in (patches.get("bioprojects") or {}): patched = True
            
            if is_missing and not patched:
                if bp not in missing_bps:
                    missing_bps[bp] = []
                missing_bps[bp].append(srr)
                
        for bp, srrs in missing_bps.items():
            examples = ", ".join(srrs[:3]) + ("..." if len(srrs)>3 else "")
            print(f"{species:<30} {bp:<20} {len(srrs):<10} {examples}")

if __name__ == "__main__":
    main()
