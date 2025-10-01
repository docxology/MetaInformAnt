# Ontology Module

The `ontology` module provides tools for functional annotation and semantic analysis using biological ontologies like Gene Ontology (GO).

## Overview

This module handles functional annotation, semantic similarity, and ontology-based analysis of biological data.

## Submodules

### Gene Ontology (`go.py`)
Gene Ontology annotation and enrichment analysis.

**Key Features:**
- GO term annotation retrieval
- Enrichment analysis algorithms
- Semantic similarity calculation
- GO hierarchy navigation

**Usage:**
```python
from metainformant.ontology import go

# GO enrichment analysis
gene_list = ["GENE1", "GENE2", "GENE3"]
enriched = go.enrich_genes(gene_list, background_genes)
similarity = go.semantic_similarity("GO:0008150", "GO:0009987")
```

### OBO Format (`obo.py`)
OBO (Open Biological and Biomedical Ontologies) format parsing.

**Key Features:**
- OBO file parsing and validation
- Ontology structure extraction
- Term relationship analysis
- Custom ontology support

**Usage:**
```python
from metainformant.ontology import obo

# Load ontology
ontology = obo.load_obo("go.obo")
terms = obo.extract_terms(ontology)
relationships = obo.get_relationships(ontology)
```

### Query Interface (`query.py`)
Unified interface for ontology queries and operations.

**Key Features:**
- Cross-ontology queries
- Annotation retrieval
- Term mapping and conversion
- Batch processing

**Usage:**
```python
from metainformant.ontology import query

# Multi-ontology queries
results = query.query_ontologies(["GO:0008150", "MP:0001262"])
annotations = query.get_annotations("GENE1")
```

### Type System (`types.py`)
Ontology type definitions and validation.

**Key Features:**
- Type definitions for ontology objects
- Validation and type checking
- Serialization support
- Type conversion utilities

**Usage:**
```python
from metainformant.ontology import types

# Type validation
term = types.GOTerm(id="GO:0008150", name="biological_process")
is_valid = types.validate_term(term)
```

## Integration with Other Modules

### With Networks Module
```python
from metainformant.networks import ppi
from metainformant.ontology import go

# Functional analysis of network modules
modules = ppi.find_modules(protein_network)
functional_annotation = go.analyze_modules(modules)
```

### With Protein Module
```python
from metainformant.protein import proteomes
from metainformant.ontology import go

# Functional annotation of proteomes
proteins = proteomes.get_proteome("UP000005640")
annotations = go.annotate_proteins(proteins)
```

## Performance Features

- Efficient ontology traversal algorithms
- Caching of expensive computations
- Batch processing for large gene lists
- Memory-optimized data structures

## Testing

Comprehensive tests cover:
- Ontology parsing accuracy
- Enrichment analysis correctness
- Semantic similarity validation
- Integration with external databases

## Dependencies

- Optional: GO database access, OBO parsing libraries

This module provides essential tools for functional annotation and semantic analysis in biological research.
