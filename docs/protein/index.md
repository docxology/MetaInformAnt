### Protein: Overview

Modular tools for protein sequences and structures. Examples write to `output/` by default and accept caller-provided paths to override.

- Sequences (`metainformant.protein.sequence.sequences`)
  - Parse FASTA, validate sequences, amino-acid composition, k-mer frequencies
  - Example:
    ```python-snippet
    from pathlib import Path
    from metainformant.protein.sequence.sequences import calculate_aa_composition, read_fasta

    recs = read_fasta(Path("data/protein/example.faa"))
    for rid, seq in recs.items():
        comp = calculate_aa_composition(seq)
        # write JSON to output/
        from metainformant.core.io import dump_json, ensure_directory
        ensure_directory("output/protein")
        dump_json(comp, Path("output/protein") / f"{rid}.comp.json", indent=2)
    ```

- UniProt (`metainformant.protein.database.uniprot`)
  - ID mapping API polling, and FASTA retrieval
  - Example:
    ```python-snippet
    from pathlib import Path
    from metainformant.protein.database.uniprot import fetch_uniprot_fasta, map_ids_uniprot
    mapping = map_ids_uniprot(["P69905"])  # {"P69905": "P69905"}
    fasta_text = fetch_uniprot_fasta("P69905")
    Path("output/protein").mkdir(parents=True, exist_ok=True)
    (Path("output/protein")/"P69905.faa").write_text(fasta_text)
    ```

- PDB (`metainformant.protein.structure.pdb`)
  - Download structures in PDB/CIF
  - Example:
    ```python-snippet
    from pathlib import Path
    from metainformant.protein.structure.pdb import fetch_pdb_structure
    out = fetch_pdb_structure("1CRN", Path("output/protein/structures"), fmt="pdb")
    ```

- AlphaFold (`metainformant.protein.structure.alphafold`)
  - Build model URLs and fetch models
  - Example:
    ```python-snippet
    from pathlib import Path
    from metainformant.protein.structure.alphafold import fetch_alphafold_model
    fetch_alphafold_model("P69905", Path("output/protein/alphafold"), version=4, fmt="pdb")
    ```

- InterPro (`metainformant.protein.database.interpro`)
  - Fetch InterPro entries for a UniProt accession; hierarchy/statistics/similar-entry helpers raise `NotImplementedError`
  - Example:
    ```python-snippet
    from metainformant.protein.database.interpro import fetch_interpro_domains
    entries = fetch_interpro_domains("P69905")
    ```

- Structure alignment (`metainformant.protein.structure.general`, `structure.io`)
  - Read CA coordinates from PDB; compute Kabsch RMSD
  - Example:
    ```python-snippet
    import numpy as np
    from pathlib import Path
    from metainformant.protein.structure.general import compute_rmsd_kabsch
    from metainformant.protein.structure.io import read_pdb_ca_coordinates

    ca = read_pdb_ca_coordinates(Path("data/protein/example.pdb"))
    A = np.array(ca)
    # ... build B ... then
    rmsd = compute_rmsd_kabsch(A, B)
    ```

- CLI
  - Taxon IDs: `uv run python -m metainformant protein taxon-ids --file tests/data/protein/taxon_id_list.txt`
  - Composition: `uv run python -m metainformant protein comp --fasta data/protein/example.faa`

Notes
- Networked APIs (UniProt, PDB, InterPro, AlphaFold) are used via simple, modular functions so you can control I/O. `PROT_TIMEOUT` sets the HTTP timeout in seconds for these protein clients. Tests must use the real network; if offline, skip gracefully.
- For reproducibility, prefer writing artifacts under `output/` and keep deterministic seeds where applicable.
