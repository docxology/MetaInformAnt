import yaml
import pandas as pd
from pathlib import Path

CONFIG_DIR = Path("config/amalgkit")
AMALGKIT_DIR = Path("output/amalgkit")

missing_bps = {
    "harpegnathos_saltator": ["PRJNA747797", "PRJNA722884", "PRJNA327093", "PRJNA327090", "PRJNA327094"],
    "monomorium_pharaonis": ["PRJNA327091", "PRJDB3164"],
    "ooceraea_biroi": ["PRJNA1010363", "PRJNA230648"],
    "pbarbatus": ["PRJDB4312", "PRJDB3493"],
    "solenopsis_invicta": ["PRJNA847170", "PRJNA522069", "PRJNA49629"],
    "temnothorax_americanus": ["PRJEB76961"],
    "temnothorax_longispinosus": ["PRJEB76961", "PRJEB68328", "PRJEB4368"],
    "vollenhovia_emeryi": ["PRJDB3517"],
    "wasmannia_auropunctata": ["PRJDB3442"],
}

for species, bps in missing_bps.items():
    meta_path = AMALGKIT_DIR / species / "work/metadata/metadata.tsv"
    if not meta_path.exists():
        continue
        
    try:
        df = pd.read_csv(meta_path, sep="\t", low_memory=False)
    except Exception as e:
        print(f"Error reading {species}: {e}")
        continue
        
    for bp in bps:
        proj_df = df[df['bioproject'] == bp]
        if proj_df.empty:
            continue
            
        print(f"\n{'='*80}")
        print(f"Species: {species} | BioProject: {bp} | Samples: {len(proj_df)}")
        
        # Look at the first valid non-empty value for useful columns
        for col in ["study_title", "study_abstract", "title", "sample_title", "description"]:
            if col in proj_df.columns:
                valid_vals = proj_df[col].dropna().unique()
                if len(valid_vals) > 0:
                    print(f"--- {col.upper()} ---")
                    # Print first 2 unique values to get an idea
                    for v in valid_vals[:2]:
                        print(str(v)[:300])
