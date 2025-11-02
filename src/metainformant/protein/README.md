# Protein Analysis Module

The `protein` module provides comprehensive tools for proteomic analysis, including sequence manipulation, structure integration, and functional annotation.

## Overview

This module handles protein sequence analysis and integrates with various protein databases and structural data sources.

## Submodules

### Proteome Retrieval (`proteomes.py`)
Tools for retrieving and working with complete proteome datasets.

**Key Features:**
- UniProt proteome retrieval and parsing
- Proteome annotation and metadata extraction
- Cross-species proteome comparison
- Proteome quality assessment

**Usage:**
```python
from metainformant.protein import read_taxon_ids
from pathlib import Path

# Read taxon IDs from file
taxon_file = Path("taxon_ids.txt")
taxon_ids = read_taxon_ids(taxon_file)

# For full proteome analysis, use protein sequences module
from metainformant.protein import parse_fasta
proteome_dict = parse_fasta(Path("proteome.fasta"))
```

### Protein Sequences (`sequences.py`)
Protein sequence manipulation and analysis.

**Key Features:**
- Amino acid composition analysis
- Physicochemical property calculation
- Sequence motif detection
- Post-translational modification prediction

**Usage:**
```python
from metainformant.protein import (
    is_valid_protein_sequence,
    calculate_aa_composition,
    kmer_frequencies,
    parse_fasta
)
from pathlib import Path

# Validate sequence
seq = "MKVLWAALLVTFLAGCQAKVE"
is_valid = is_valid_protein_sequence(seq)

# Calculate amino acid composition
composition = calculate_aa_composition(seq)

# Calculate k-mer frequencies
kmers = kmer_frequencies(seq, k=3)

# Parse FASTA file
sequences = parse_fasta(Path("proteins.fasta"))
```

### Structure Integration (`structure.py`, `pdb.py`)
Protein structure data integration and analysis.

**Key Features:**
- PDB file parsing and manipulation
- Structure validation and quality assessment
- Structure-based function prediction
- Structure comparison and alignment

**Usage:**
```python
from metainformant.protein import fetch_pdb_structure
from pathlib import Path

# Download PDB structure
pdb_id = "1A2B"
output_dir = Path("output/pdb")
pdb_path = fetch_pdb_structure(pdb_id, output_dir, fmt="pdb")

# For structure analysis, use compute_rmsd_kabsch
from metainformant.protein import compute_rmsd_kabsch
import numpy as np

# Compare two structures (requires coordinate arrays)
# coords_ref = np.array([...])  # Reference coordinates
# coords_mobile = np.array([...])  # Mobile coordinates
# rmsd = compute_rmsd_kabsch(coords_ref, coords_mobile)
```

### Functional Annotation (`interpro.py`, `uniprot.py`)
Integration with InterPro and UniProt for functional annotation.

**Key Features:**
- InterPro domain and family annotation
- UniProt keyword and GO term mapping
- Enzyme classification and pathway mapping
- Literature reference integration

**Usage:**
```python
from metainformant.protein import map_ids_uniprot
from metainformant.protein.interpro import fetch_interpro_domains

# Map IDs to UniProt accessions
protein_ids = ["P12345", "Q67890"]
id_mapping = map_ids_uniprot(protein_ids)

# Fetch InterPro domains for a UniProt accession
uniprot_acc = "P12345"
domains = fetch_interpro_domains(uniprot_acc)
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import translation
from metainformant.protein import calculate_aa_composition, is_valid_protein_sequence

# Translate DNA and analyze protein
dna_sequence = "ATGGCC..."
protein_sequence = translation.translate_dna(dna_sequence)
is_valid = is_valid_protein_sequence(protein_sequence)
composition = calculate_aa_composition(protein_sequence)
```

### With Visualization Module
```python
from metainformant.protein import simple_helix_coil_propensity
from metainformant.visualization import lineplot

# Visualize protein properties
seq = "MKVLWAALLVTFLAGCQAKVE"
propensity = simple_helix_coil_propensity(seq)
ax = lineplot(None, propensity)
ax.set_title("Helix-Coil Propensity Profile")
```

## Performance Features

- Memory-efficient processing of large proteomes
- Parallel sequence analysis
- Caching of expensive computations
- Streaming processing for large datasets

## Testing

Comprehensive tests cover:
- Sequence format validation
- Annotation accuracy verification
- Structure parsing correctness
- Integration with external databases

## Dependencies

- Biopython for sequence objects
- Optional: UniProt API access, PDB parsing libraries

This module provides essential tools for protein sequence analysis and functional annotation.
