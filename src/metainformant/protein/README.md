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
from metainformant.protein import proteomes

# Retrieve proteome data
proteome = proteomes.get_proteome("UP000005640")  # Human proteome
sequences = proteomes.extract_sequences(proteome)

# Analyze proteome statistics
stats = proteomes.proteome_statistics(proteome)
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
from metainformant.protein import sequences

# Sequence analysis
seq = "MKVLWAALLVTFLAGCQAKVE"
composition = sequences.aa_composition(seq)
molecular_weight = sequences.molecular_weight(seq)
hydrophobicity = sequences.hydrophobicity_profile(seq)
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
from metainformant.protein import pdb

# Load and analyze structure
structure = pdb.load_pdb("protein.pdb")
chains = pdb.extract_chains(structure)
resolution = pdb.get_resolution(structure)
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
from metainformant.protein import uniprot, interpro

# Retrieve annotations
protein_id = "P12345"
uniprot_data = uniprot.get_protein_data(protein_id)
domains = interpro.get_domains(uniprot_data["sequence"])
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import translation
from metainformant.protein import proteomes

# Translate DNA and analyze protein
dna_sequence = "ATGGCC..."
protein_sequence = translation.translate(dna_sequence)
proteome_analysis = proteomes.analyze_protein(protein_sequence)
```

### With Visualization Module
```python
from metainformant.protein import sequences
from metainformant.visualization import lineplot

# Visualize protein properties
seq = "MKVLWAALLVTFLAGCQAKVE"
hydrophobicity = sequences.hydrophobicity_profile(seq)
lineplot(hydrophobicity, title="Hydrophobicity Profile")
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
