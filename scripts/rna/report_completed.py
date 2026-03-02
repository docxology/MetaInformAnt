import glob, os
print("| Species | Completed Samples | Total SRA Samples |")
print("|---|---|---|")
for cfg in sorted(glob.glob("config/amalgkit/amalgkit_*.yaml")):
    species = os.path.basename(cfg).replace("amalgkit_", "").replace(".yaml", "")
    if any(skip in species for skip in ["template", "test", "cross_species", "apis_mellifera_all"]):
        continue
        
    meta = f"output/amalgkit/{species}/work/metadata/metadata.tsv"
    total = "Pending Metadata"
    if os.path.isfile(meta):
        with open(meta) as f:
            lines = sum(1 for line in f)
            total = str(max(0, lines - 1))
            
    comp = 0
    quant_dir = f"output/amalgkit/{species}/quant/quant"
    if os.path.isdir(quant_dir):
        comp = len(glob.glob(f"{quant_dir}/*/*_abundance.tsv"))
    # Also check the old/work directory in case
    else:
        comp = len(glob.glob(f"output/amalgkit/{species}/**/*_abundance.tsv", recursive=True))
        
    print(f"| `{species}` | **{comp}** | {total} |")
