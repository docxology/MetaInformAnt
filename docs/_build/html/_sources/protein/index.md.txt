### Protein: Overview

Modular tools for protein sequences and structures. Examples write to `output/` by default and accept caller-provided paths to override.

- Sequences (`metainformant.protein.sequences`)
  - Parse FASTA, validate sequences, amino-acid composition, k-mer frequencies
  - Example:
    ```python
    from pathlib import Path
    from metainformant.protein.sequences import parse_fasta, calculate_aa_composition

    recs = parse_fasta(Path("data/protein/example.faa"))
    for rid, seq in recs.items():
        comp = calculate_aa_composition(seq)
        # write JSON to output/
        from metainformant.core.io import dump_json, ensure_directory
        ensure_directory("output/protein")
        dump_json(comp, Path("output/protein") / f"{rid}.comp.json", indent=2)
    ```

- UniProt (`metainformant.protein.uniprot`)
  - ID mapping API polling, and FASTA retrieval
  - Example:
    ```python
    from metainformant.protein.uniprot import map_ids_uniprot, fetch_uniprot_fasta
    mapping = map_ids_uniprot(["P69905"])  # {"P69905": "P69905"}
    fasta_text = fetch_uniprot_fasta("P69905")
    Path("output/protein").mkdir(parents=True, exist_ok=True)
    (Path("output/protein")/"P69905.faa").write_text(fasta_text)
    ```

- PDB (`metainformant.protein.pdb`)
  - Download structures in PDB/CIF
  - Example:
    ```python
    from pathlib import Path
    from metainformant.protein.pdb import fetch_pdb_structure
    out = fetch_pdb_structure("1CRN", Path("output/protein/structures"), fmt="pdb")
    ```

- AlphaFold (`metainformant.protein.alphafold`)
  - Build model URLs and fetch models
  - Example:
    ```python
    from pathlib import Path
    from metainformant.protein.alphafold import fetch_alphafold_model
    fetch_alphafold_model("P69905", Path("output/protein/alphafold"), version=4, fmt="pdb")
    ```

- InterPro (`metainformant.protein.interpro`)
  - Fetch InterPro entries for a UniProt accession
  - Example:
    ```python
    from metainformant.protein.interpro import fetch_interpro_domains
    entries = fetch_interpro_domains("P69905")
    ```

- Structure alignment (`metainformant.protein.structure`, `structure_io`)
  - Read CA coordinates from PDB; compute Kabsch RMSD
  - Example:
    ```python
    import numpy as np
    from pathlib import Path
    from metainformant.protein.structure import compute_rmsd_kabsch
    from metainformant.protein.structure_io import read_pdb_ca_coordinates

    ca = read_pdb_ca_coordinates(Path("data/protein/example.pdb"))
    A = np.array(ca)
    # ... build B ... then
    rmsd = compute_rmsd_kabsch(A, B)
    ```

- CLI
  - Taxon IDs: `uv run python -m metainformant protein taxon-ids --file tests/data/protein/taxon_id_list.txt`
  - Composition: `uv run python -m metainformant protein comp --fasta data/protein/example.faa`

Notes
- Networked APIs (UniProt, PDB, InterPro, AlphaFold) are used via simple, modular functions so you can control I/O. Tests must use the real network; if offline, skip gracefully.
- For reproducibility, prefer writing artifacts under `output/` and keep deterministic seeds where applicable.
