# Ontology Core

Core ontology data structures, OBO format parsing, and Gene Ontology analysis including enrichment, semantic similarity, and hierarchy traversal.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `go`, `obo`, `types` |
| `types.py` | Core dataclasses: `Term`, `Relationship`, `Ontology` |
| `obo.py` | OBO format file parsing into `Ontology` objects |
| `go.py` | GO enrichment, semantic similarity, hierarchy utilities |

## Key Classes

| Class | Description |
|-------|-------------|
| `Term` | Ontology term with metadata, synonyms, and relationships |
| `Relationship` | Typed relationship between two ontology terms |
| `Ontology` | Container for terms and relationships with lookup methods |

## Key Functions

| Function | Description |
|----------|-------------|
| `parse_obo()` | Parse an OBO format file into an `Ontology` object |
| `create_term()` | Factory for creating `Term` instances |
| `create_relationship()` | Factory for creating `Relationship` instances |
| `create_ontology()` | Factory for creating `Ontology` instances |
| `load_go_obo()` | Load a Gene Ontology OBO file |
| `enrich_genes()` | GO enrichment analysis for a gene list |
| `semantic_similarity()` | Semantic similarity between GO terms |
| `calculate_term_ic()` | Information content of ontology terms |
| `build_hierarchy_dict()` | Build parent-child hierarchy dictionary |

## Usage

```python
from metainformant.ontology.core import go, obo, types

ontology = obo.parse_obo("data/go.obo")
term = types.create_term(term_id="GO:0008150", name="biological_process")
enrichment = go.enrich_genes(gene_list, ontology, background=bg_genes)
sim = go.semantic_similarity(ontology, "GO:0008150", "GO:0003674")
```
