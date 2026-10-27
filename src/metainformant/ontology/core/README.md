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
| `parse_obo_definition()` | Split an OBO `def:` value into text and dbxrefs |
| `parse_obo_synonym()` | Split an OBO `synonym:` value into text, scope (`EXACT`, `BROAD`, `NARROW`, `RELATED`) and dbxrefs |
| `create_term()` | Factory for creating `Term` instances |
| `create_relationship()` | Factory for creating `Relationship` instances |
| `create_ontology()` | Factory for creating `Ontology` instances |
| `load_go_obo()` | Load a Gene Ontology OBO file |
| `enrich_genes()` | GO over-representation analysis for a gene list: one row per tested term with a Bonferroni-adjusted p-value and a boolean `significant` column (`alpha`); only genes in the background count, for both the Fisher and the hypergeometric method |
| `semantic_similarity()` | Semantic similarity between GO terms: Resnik is an information content in bits (not bounded by 1), Lin and Jiang-Conrath lie in [0, 1] |
| `calculate_term_ic()` | Information content of ontology terms |
| `build_hierarchy_dict()` | Build parent-child hierarchy dictionary |

## Usage

```python
from metainformant.ontology.core import go, obo, types

ontology = obo.parse_obo("data/go.obo")
term = types.create_term(id="GO:0008150", name="biological_process")

# go_annotations: GO term ID -> set of gene IDs; bg_genes is the gene universe (None uses every annotated gene)
enrichment = go.enrich_genes(gene_list, bg_genes, go_annotations, alpha=0.05, method="fisher")
significant = enrichment[enrichment["significant"]]  # bonferroni_p <= alpha; every tested term stays in the frame

# term_counts: GO term ID -> annotation count; the hierarchy maps each term to its ancestors
term_ic = go.calculate_term_ic(ontology, term_counts, total_annotations)
hierarchy = go.build_hierarchy_dict(ontology)
sim = go.semantic_similarity("GO:0008150", "GO:0003674", term_ic, hierarchy, method="lin")
```

`enrich_genes` returns the columns `go_term`, `p_value`, `enrichment`, `observed`, `expected`, `term_size` (the term's
genes that are in the background), `query_size`, `bonferroni_p` and `significant`; without scipy the p-values are the
exact upper tail. A `synonym:` line is read with `parse_obo_synonym`, so quotes, backslashes and a `!` inside the
synonym text survive; ontologies parsed from OBO round-trip through `query.serialize.save_ontology` and
`load_ontology` in JSON as well as in OBO format.
