# Ontology Module

The `ontology` module provides tools for functional annotation and semantic analysis using biological ontologies like Gene Ontology (GO).

## Overview

This module handles ontology parsing, hierarchy traversal, and term queries. Provides lightweight, efficient tools for working with OBO-format ontologies without requiring external database connections.

## Key Components

### Gene Ontology Loading (`go.py`)
Load and work with Gene Ontology from OBO files.

**Usage:**
```python
from metainformant.ontology import load_go_obo, write_go_summary
from pathlib import Path

# Load GO ontology from OBO file
onto = load_go_obo("data/go-basic.obo")

# Get basic statistics
print(f"Number of terms: {onto.num_terms()}")

# Write summary
summary_path = write_go_summary(onto)
print(f"Summary written to: {summary_path}")
```

### Ontology Types (`types.py`)
Core data structures for ontology representation.

**Usage:**
```python
from metainformant.ontology.types import Term, Ontology

# Create a term
term = Term(
    term_id="GO:0008150",
    name="biological_process",
    namespace="biological_process",
    definition="Any process accomplished by biological systems",
    is_a_parents=["GO:0003674"]
)

# Create ontology and add term
onto = Ontology()
onto.add_term(term)

# Check term existence
if onto.has_term("GO:0008150"):
    print("Term exists in ontology")
```

### Ontology Queries (`query.py`)
Traverse ontology hierarchies and extract subgraphs.

**Usage:**
```python
from metainformant.ontology.query import ancestors, descendants, subgraph

# Get all ancestor terms (broader terms)
onto = load_go_obo("go.obo")
ancestors_set = ancestors(onto, "GO:0008150")
print(f"Ancestors of biological_process: {len(ancestors_set)} terms")

# Get all descendant terms (more specific terms)
descendants_set = descendants(onto, "GO:0008150")
print(f"Descendants: {len(descendants_set)} terms")

# Extract subgraph rooted at specific terms
roots = ["GO:0008150", "GO:0003674"]
sub_onto = subgraph(onto, roots)
print(f"Subgraph size: {sub_onto.num_terms()} terms")
```

### OBO Parsing (`obo.py`)
Parse OBO format files into Ontology objects.

**Usage:**
```python
from metainformant.ontology.obo import parse_obo

# Parse OBO file
onto = parse_obo("go-basic.obo")

# Access terms
for term_id, term in onto.terms.items():
    print(f"{term_id}: {term.name}")
    if term.is_a_parents:
        print(f"  Parents: {term.is_a_parents}")
```

**Supported OBO Fields:**
- `id`: Term identifier
- `name`: Term name
- `namespace`: Ontology namespace
- `def`: Term definition
- `alt_id`: Alternative identifiers
- `is_a`: Parent relationships

## Integration with Other Modules

### With Networks Module
```python
from metainformant.networks import detect_communities
from metainformant.ontology import load_go_obo

# Functional analysis of network modules
communities = detect_communities(protein_network)

# Load GO for enrichment analysis
go_onto = load_go_obo("go-basic.obo")
# Use GO for functional annotation
```

### With Protein Module
```python
from metainformant.protein import parse_fasta
from metainformant.ontology import load_go_obo

# Load proteome
proteins = parse_fasta(Path("proteome.fasta"))

# Load GO for functional annotation
go_onto = load_go_obo("go-basic.obo")
# Use GO for protein functional annotation
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
