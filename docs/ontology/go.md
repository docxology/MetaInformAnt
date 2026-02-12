# Gene Ontology (GO) Analysis

The ontology module provides comprehensive tools for working with Gene Ontology (GO) ontologies.

## Loading GO Ontologies

```python
from metainformant.ontology.core.go import load_go_obo, validate_go_ontology, write_go_summary

# Load GO ontology from OBO file
onto = load_go_obo("data/go-basic.obo")

# Validate ontology structure
is_valid, errors = validate_go_ontology(onto)
if not is_valid:
    print(f"Validation issues: {errors}")

# Write summary statistics
summary_path = write_go_summary(onto)
```

## GO Term Queries

```python
from metainformant.ontology.core.go import (
    ancestors, descendants, common_ancestors,
    path_to_root, distance, find_term_by_name,
    filter_by_namespace, get_roots, get_leaves
)

onto = load_go_obo("go.obo")

# Find terms by name
matches = find_term_by_name(onto, "biological process")
print(f"Found {len(matches)} matching terms")

# Get namespace-specific terms
bp_onto = filter_by_namespace(onto, "biological_process")
print(f"Biological process terms: {bp_onto.num_terms()}")

# Get root and leaf terms
roots = get_roots(onto, namespace="biological_process")
leaves = get_leaves(onto, namespace="biological_process")
```

## GO Annotation File Processing

The module includes utilities for counting GO annotation scripts:

```python
from pathlib import Path
from metainformant.ontology.core.go import count_go_scripts

n = count_go_scripts(Path("tests/data/ontology/GO_v3"))
print(f"Found {n} GO annotation scripts")
```

## Planned Features

The following features are planned but not yet implemented:
- GO annotation file parsing (GAF/GPAD format)
- Gene enrichment analysis
- Semantic similarity calculations

See the main module README for current functionality and limitations.
