import subprocess
from pathlib import Path

import pandas as pd

species_list = ["pogonomyrmex_barbatus", "temnothorax_americanus"]
data_dir = Path("projects/hymenoptera_amalgkit/data")

for species in species_list:
    work_dir = data_dir / species / "work"
    meta_file = work_dir / "metadata" / "metadata_selected.tsv"

    if not meta_file.exists():
        print(f"Metadata not found for {species}")
        continue

    df = pd.read_csv(meta_file, sep="\t", dtype=str)

    # Drop rows without abundance.tsv
    quant_dir = work_dir / "quant"

    def has_quant(row):
        run_id = row["run"]
        abundance_file = quant_dir / run_id / f"{run_id}_abundance.tsv"
        return abundance_file.exists()

    df_filtered = df[df.apply(has_quant, axis=1)].copy()

    print(f"{species}: Filtered from {len(df)} to {len(df_filtered)} samples.")

    # Backup and save
    meta_file.rename(str(meta_file) + ".bak")
    df_filtered.to_csv(meta_file, sep="\t", index=False)

    # Run tissue patching
    print(f"Running tissue patching for {species}...")
    subprocess.run(
        [
            ".venv/bin/python",
            "scripts/rna/normalize_tissue_metadata.py",
            "--input",
            str(meta_file),
            "--output",
            str(meta_file),
        ]
    )

print("Done.")
