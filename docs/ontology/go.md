# Gene Ontology (GO) Analysis

The ontology module provides tools for working with Gene Ontology (GO) ontologies: loading and validating OBO files, querying the term graph, over-representation analysis, and information-content semantic similarity. GO loading, enrichment and similarity live in `metainformant.ontology.core.go`; graph queries live in `metainformant.ontology.query`.

## Loading GO Ontologies

```python
from pathlib import Path
from metainformant.ontology.core.go import load_go_obo, validate_go_ontology, write_go_summary

# A tiny GO-style OBO file; in practice point load_go_obo at a real file such as go-basic.obo.
# The examples below reuse mini_go.obo.
Path("mini_go.obo").write_text("""format-version: 1.2
ontology: go
remark: Gene Ontology (tiny example)

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0006810
name: transport
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0055085
name: transmembrane transport
namespace: biological_process
is_a: GO:0006810 ! transport
is_a: GO:0009987 ! cellular process
""")

# Load GO ontology from OBO file
onto = load_go_obo("mini_go.obo")

# Validate ontology structure: returns (is_valid, list_of_issues)
is_valid, errors = validate_go_ontology(onto)
if not is_valid:
    print(f"Validation issues: {errors}")

# Write summary statistics (defaults to output/go_summary.txt; a .json path writes JSON)
summary_path = write_go_summary(onto)
```

## GO Term Queries

```python
from metainformant.ontology.core.go import load_go_obo
from metainformant.ontology.query import (
    ancestors, descendants, common_ancestors,
    path_to_root, distance, find_term_by_name,
    filter_by_namespace, get_roots, get_leaves
)

onto = load_go_obo("mini_go.obo")

# Find terms by name: a case-insensitive substring match that returns term IDs
matches = find_term_by_name(onto, "transport")
print(f"Found {len(matches)} matching terms")  # GO:0006810 and GO:0055085

# Get namespace-specific terms: the filtered ontology keeps the parent/child links among its terms
bp_onto = filter_by_namespace(onto, "biological_process")
print(f"Biological process terms: {len(bp_onto)}")  # 4

# Root and leaf terms of one namespace
roots = get_roots(bp_onto)      # {'GO:0008150'}
leaves = get_leaves(bp_onto)    # {'GO:0055085'}

# Traverse the is_a graph; a term is not counted among its own ancestors or descendants
ancestors(onto, "GO:0055085")        # {'GO:0006810', 'GO:0009987', 'GO:0008150'}
descendants(onto, "GO:0008150")      # {'GO:0006810', 'GO:0009987', 'GO:0055085'}
common_ancestors(onto, "GO:0055085", "GO:0009987")  # {'GO:0008150'}
path_to_root(onto, "GO:0055085")     # a shortest path to a root, e.g. ['GO:0055085', 'GO:0006810', 'GO:0008150']
distance(onto, "GO:0055085", "GO:0008150")  # 2 edges; -1 when the terms are not connected
```

`get_roots` and `get_leaves` take only the ontology (and a `relation_type`, `"is_a"` by default), not a namespace; to look at one namespace, filter first. `filter_by_namespace` returns a new `Ontology` holding that namespace's terms, the relationships among them and the parent/child links among them, so hierarchy queries such as `get_roots`, `get_leaves` and `ancestors` work on the filtered ontology. A term whose parents all lie in other namespaces becomes a root of it. `subgraph` and `get_subontology` carry the parent/child links among their terms in the same way.

## GO Enrichment

`enrich_genes(genes, background, annotations, alpha=0.05, method="fisher")` runs over-representation analysis. `annotations` maps GO term IDs to sets of gene IDs, and `background` is the gene universe (`None` uses every gene in `annotations`). Only genes in the background count: query genes outside it are dropped with a warning, and so are a term's members outside it. The method is `"fisher"` (one-sided Fisher's exact test) or `"hypergeometric"`; both test the same upper tail on the same background-restricted counts, so they give the same p-value, and without scipy the p-value is the exact upper tail.

```python
from metainformant.ontology.core.go import enrich_genes

annotations = {
    "GO:0006810": {"g1", "g2", "g3", "g4"},
    "GO:0055085": {"g1", "g2"},
    "GO:0009987": {"g1", "g2", "g5", "g6", "g7", "g8"},
}
background = [f"g{i}" for i in range(1, 21)]

results = enrich_genes(["g1", "g2", "g3"], background, annotations, alpha=0.05, method="fisher")
significant = results[results["significant"]]
print(results[["go_term", "p_value", "enrichment", "observed", "bonferroni_p", "significant"]])
```

The result is a pandas `DataFrame` sorted by `p_value` with one row per tested term and the columns `go_term`, `p_value`, `enrichment`, `observed`, `expected`, `term_size` (the term's genes that are in the background), `query_size`, `bonferroni_p` and `significant`. Only terms that annotate at least one query gene are tested. `bonferroni_p` is `p_value` times the number of tested terms, capped at 1, and `alpha` sets the boolean `significant` column (`bonferroni_p <= alpha`); every tested term is still returned, so near misses stay visible. The frame is empty when nothing can be tested, for example when no query gene is in the background.

## Semantic Similarity

Similarity between two terms uses their information content (IC). `calculate_term_ic` turns annotation counts into IC values (-log2 of the term's share of all annotations), `build_hierarchy_dict` maps each term to its ancestors, and `semantic_similarity` combines them through the most informative common ancestor.

```python
from metainformant.ontology.core.go import (
    load_go_obo, calculate_term_ic, build_hierarchy_dict, semantic_similarity,
)

onto = load_go_obo("mini_go.obo")

term_counts = {"GO:0008150": 20, "GO:0009987": 6, "GO:0006810": 4, "GO:0055085": 2}
term_ic = calculate_term_ic(onto, term_counts, total_annotations=20)
hierarchy = build_hierarchy_dict(onto)

resnik = semantic_similarity("GO:0055085", "GO:0006810", term_ic, hierarchy)  # 2.32 bits
lin = semantic_similarity("GO:0055085", "GO:0006810", term_ic, hierarchy, method="lin")  # 0.82
```

The methods are `"resnik"` (the default), `"lin"` and `"jiang-conrath"`. Resnik is the IC of the most informative common ancestor in the units of `term_ic`, which are bits when `term_ic` comes from `calculate_term_ic`, so it is not bounded by 1; divide by `max(term_ic.values())` for a 0-1 scale. Lin (`2 * IC(MICA) / (IC(term1) + IC(term2))`) and Jiang-Conrath (the distance `IC(term1) + IC(term2) - 2 * IC(MICA)` rescaled by twice the largest IC in `term_ic`) lie between 0 and 1. A score is 0.0 when a term is missing from `term_ic` or the two terms share no ancestor, and an unknown method raises `ValueError`.

## GO Annotation Files

`genes_to_go_annotations` in `metainformant.ontology.annotation` builds the `{go_id: {genes}}` mapping that `enrich_genes` expects. With `source="gaf"` it reads a local GAF 2.x file (plain or `.gz`) and keeps only the rows for the genes you pass; the default `source="api"` queries QuickGO instead and needs network access.

```python
from pathlib import Path
from metainformant.ontology.annotation.annotate import genes_to_go_annotations

# Two GAF rows: gene symbol in column 3, GO ID in column 5, aspect (P/F/C) in column 9
rows = [
    ["UniProtKB", "P1", "g1", "", "GO:0006810", "PMID:1", "IEA", "", "P", "", "", "protein", "taxon:7460", "20240101", "UniProt"],
    ["UniProtKB", "P2", "g2", "", "GO:0005215", "PMID:1", "IEA", "", "F", "", "", "protein", "taxon:7460", "20240101", "UniProt"],
]
Path("mini.gaf").write_text("!gaf-version: 2.2\n" + "\n".join("\t".join(row) for row in rows) + "\n")

annotations = genes_to_go_annotations(["g1", "g2"], source="gaf", annotation_file="mini.gaf")
# {'GO:0006810': {'g1'}, 'GO:0005215': {'g2'}}

# Keep one aspect only
function_only = genes_to_go_annotations(
    ["g1", "g2"], source="gaf", annotation_file="mini.gaf", go_aspects=["molecular_function"]
)
```

## GO Analysis Scripts

`count_go_scripts` counts the GO-related analysis scripts (`.py`, `.r`, `.sh` and `.pl` files whose names contain "go" or "ontology") in a directory, which must exist:

```python
from pathlib import Path
from metainformant.ontology.core.go import count_go_scripts

n = count_go_scripts(Path("tests/data/ontology/GO_v3"))
print(f"Found {n} GO analysis scripts")  # 3 when run from the repository root
```

## Not Yet Implemented

- GPAD-format annotation parsing, and a standalone GAF/GPAD file parser (GAF files are read only through `genes_to_go_annotations`, filtered to the genes you pass)

See the [ontology README](README.md) for the rest of the module.
