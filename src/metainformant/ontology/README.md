# Ontology Module

The `ontology` module provides tools for functional annotation and ontology analysis using biological ontologies like Gene Ontology (GO).

## Overview

This module handles ontology parsing, hierarchy traversal, term queries, and serialization. Provides lightweight, efficient tools for working with OBO-format ontologies without requiring external database connections. Includes error handling, validation, caching, and support for multiple relationship types.

### Module Architecture

```mermaid
graph TB
    subgraph "Ontology Module"
        GO[go<br/>Gene Ontology]
        OBO[obo<br/>OBO Parser]
        Query[query<br/>Term Queries]
        Serialize[serialize<br/>Serialization]
        Types[types<br/>Type System]
    end
    
    subgraph "Input"
        OBOFile[OBO File]
        Terms[GO Terms]
    end
    
    subgraph "Other Modules"
        Protein_Mod[protein]
        Networks_Mod[networks]
        Info_Mod[information]
    end
    
    OBOFile --> OBO
    OBO --> GO
    Terms --> Query
    GO --> Query
    Query --> Serialize
    GO --> Protein_Mod
    GO --> Networks_Mod
    GO --> Info_Mod
```

### Ontology Parsing and Loading

```mermaid
graph TD
    A[OBO File] --> B[Parse Header]
    B --> C[Extract Ontology Metadata]

    A --> D[Parse Terms]
    D --> E[Term Definitions]
    E --> F[Term Relationships]
    F --> G[Synonyms & Cross-references]

    A --> H[Parse Typedefs]
    H --> I[Relationship Types]
    I --> J[Relationship Properties]

    C --> K[Ontology Graph]
    F --> K
    I --> K

    K --> L[Validation]
    L --> M{Valid Structure?}

    M -->|Yes| N[Build NetworkX Graph]
    M -->|No| O[Error Reporting]

    N --> P[Index Terms]
    P --> Q[Term Lookup Tables]
    Q --> R[Ready for Queries]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style K fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style R fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "OBO Components"
        S[Header Section] -.-> B
        T[Term Stanzas] -.-> D
        U[Typedef Stanzas] -.-> H
        V[Instance Stanzas] -.-> D
    end

    subgraph "Relationship Types"
        W[is_a] -.-> F
        X[part_of] -.-> F
        Y[regulates] -.-> F
        Z[has_part] -.-> F
    end
```

### Semantic Similarity Analysis

```mermaid
graph TD
    A[Two GO Terms] --> B[Information Content]
    A --> C[Common Ancestors]

    B --> D[Term IC Calculation]
    C --> E[Lowest Common Ancestor]

    D --> F[Similarity Method]
    E --> F

    F --> G{Similarity Measure}
    G -->|Resnik| H[Shared IC]
    G -->|Lin| I[Shared IC / Average IC]
    G -->|Jiang-Conrath| J[1 - Shared IC Distance]
    G -->|Wang| K[Semantic Contribution]

    H --> L[Similarity Score]
    I --> L
    J --> L
    K --> L

    L --> M[Threshold Comparison]
    M --> N{Similar?}

    N -->|Yes| O[Functional Relationship]
    N -->|No| P[Different Functions]

    O --> Q[Annotation Transfer]
    P --> R[Distinct Functions]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style F fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style Q fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "IC Calculation"
        S[Term Frequency] -.-> D
        T[Corpus Size] -.-> D
        U[Annotation Corpus] -.-> D
    end

    subgraph "Applications"
        V[Gene Clustering] -.-> O
        W[Function Prediction] -.-> O
        X[Ortholog Detection] -.-> O
    end
```

### Functional Enrichment Analysis

```mermaid
graph TD
    A[Gene List] --> B[Background Genes]
    B --> C[GO Annotations]

    A --> D[Annotated Genes]
    D --> E[Term Counts]

    C --> F[Background Counts]
    F --> G[Term Statistics]

    E --> H{Enrichment Test}
    H -->|Fisher's Exact| I[Hypergeometric Test]
    H -->|Chi-square| J[Chi-square Test]
    H -->|Binomial| K[Binomial Test]

    I --> L[P-value Calculation]
    J --> L
    K --> L

    L --> M[Multiple Testing Correction]
    M --> N{Method}
    N -->|Bonferroni| O[Family-wise Error]
    N -->|BH FDR| P[False Discovery Rate]
    N -->|Holm| Q[Step-down Bonferroni]

    O --> R[Adjusted P-values]
    P --> R
    Q --> R

    R --> S[Significant Terms]
    S --> T[Enrichment Results]

    T --> U[Pathway Visualization]
    U --> V[Biological Insights]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style H fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style V fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Contingency Table"
        W[Annotated in List] -.-> E
        X[Not in List] -.-> E
        Y[Annotated in Background] -.-> F
        Z[Not Annotated] -.-> F
    end

    subgraph "GO Aspects"
        AA[Biological Process] -.-> C
        BB[Molecular Function] -.-> C
        CC[Cellular Component] -.-> C
    end
```

### Ontology-Based Clustering

```mermaid
graph TD
    A[Gene Set] --> B[GO Annotations]
    B --> C[Semantic Similarity Matrix]

    C --> D{Clustering Method}
    D -->|Hierarchical| E[Agglomerative Clustering]
    D -->|K-means| F[K-means on Similarity]
    D -->|Community| G[Community Detection]

    E --> H[Distance Matrix]
    F --> I[Centroid-based]
    G --> J[Modularity Optimization]

    H --> K[Cluster Assignments]
    I --> K
    J --> K

    K --> L[Cluster Validation]
    L --> M[Silhouette Score]
    L --> N[Adjusted Rand Index]

    M --> O[Optimal Clustering]
    N --> O

    O --> P[Functional Interpretation]
    P --> Q[Cluster Enrichment]
    Q --> R[GO Term Summaries]

    R --> S[Cluster Characterization]
    S --> T[Gene Function Groups]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style D fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style T fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Similarity Measures"
        U[Resnik] -.-> C
        V[Lin] -.-> C
        W[Jiang-Conrath] -.-> C
        X[Wang] -.-> C
    end

    subgraph "Distance Metrics"
        Y[1 - Similarity] -.-> H
        Z[Euclidean] -.-> H
        AA[Manhattan] -.-> H
    end
```

### Ontology Integration Framework

```mermaid
graph TD
    A[Multiple Ontologies] --> B[Ontology Mapping]
    B --> C[Term Alignment]
    C --> D[Cross-ontology Relationships]

    A --> E[Integrated Graph]
    D --> E

    E --> F[Unified Queries]
    F --> G[Term Resolution]
    G --> H[Annotation Integration]

    H --> I[Unified Enrichment]
    I --> J[Cross-ontology Analysis]

    J --> K[Comprehensive Annotations]
    K --> L[Enhanced Biological Insights]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style L fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Ontology Sources"
        M[Gene Ontology] -.-> A
        N[Human Disease] -.-> A
        O[Plant Ontology] -.-> A
        P[Cell Ontology] -.-> A
    end

    subgraph "Mapping Methods"
        Q[Exact Match] -.-> C
        R[Semantic Similarity] -.-> C
        S[Manual Curation] -.-> C
        T[Machine Learning] -.-> C
    end

    subgraph "Integration Benefits"
        U[Broader Coverage] -.-> K
        V[Contextual Information] -.-> K
        W[Cross-domain Insights] -.-> K
    end
```

## Key Components

### Gene Ontology Loading (`go.py`)
Load and work with Gene Ontology from OBO files.

**Usage:**
```python
from metainformant.ontology import load_go_obo, validate_go_ontology, write_go_summary
from pathlib import Path

# Load GO ontology from OBO file
onto = load_go_obo("data/go-basic.obo")

# Get basic statistics
print(f"Number of terms: {onto.num_terms()}")

# Validate ontology structure
is_valid, errors = validate_go_ontology(onto)
if not is_valid:
    print(f"Validation errors: {errors}")

# Write summary with statistics
summary_path = write_go_summary(onto)
print(f"Summary written to: {summary_path}")
```

### Ontology Types (`types.py`)
Core data structures for ontology representation.

**Usage:**
```python
from metainformant.ontology.types import Term, Ontology

# Create a term with relationships
term = Term(
    term_id="GO:0008150",
    name="biological_process",
    namespace="biological_process",
    definition="Any process accomplished by biological systems",
    is_a_parents=["GO:0003674"],
    relationships={"part_of": ["GO:001"]},  # Additional relationship types
    synonyms=["biological process", "BP"],
    xrefs=["Wikipedia:biological_process"],
    subsets=["goslim_generic"]
)

# Create ontology and add term
onto = Ontology()
onto.add_term(term)

# Check term existence and get term
if onto.has_term("GO:0008150"):
    term = onto.get_term("GO:0008150")
    print(f"Term namespace: {onto.get_namespace('GO:0008150')}")

# Validate ontology integrity
is_valid, errors = onto.validate()
if not is_valid:
    print(f"Validation errors: {errors}")

# Get relationships for a term
all_rels = onto.get_relationships("GO:0008150")
part_of_rels = onto.get_relationships("GO:0008150", rel_type="part_of")
```

### Ontology Queries (`query.py`)
Traverse ontology hierarchies and extract subgraphs.

**Usage:**
```python
from metainformant.ontology.query import (
    ancestors, descendants, subgraph, common_ancestors,
    path_to_root, distance, find_term_by_name,
    filter_by_namespace, get_roots, get_leaves
)

onto = load_go_obo("go.obo")

# Get all ancestor terms (broader terms)
ancestors_set = ancestors(onto, "GO:0008150")
print(f"Ancestors: {len(ancestors_set)} terms")

# Get all descendant terms (more specific terms)
descendants_set = descendants(onto, "GO:0008150")
print(f"Descendants: {len(descendants_set)} terms")

# Find common ancestors of two terms
common = common_ancestors(onto, "GO:0009987", "GO:0008150")
print(f"Common ancestors: {len(common)} terms")

# Get path from term to root
path = path_to_root(onto, "GO:0009987")
print(f"Path to root: {path}")

# Calculate distance between two terms
dist = distance(onto, "GO:0009987", "GO:0008150")
print(f"Distance: {dist}")

# Find terms by name
matches = find_term_by_name(onto, "biological process")
print(f"Found {len(matches)} matching terms")

# Filter by namespace
bp_onto = filter_by_namespace(onto, "biological_process")
print(f"Biological process terms: {bp_onto.num_terms()}")

# Get root and leaf terms
roots = get_roots(onto)
leaves = get_leaves(onto)
print(f"Root terms: {len(roots)}, Leaf terms: {len(leaves)}")

# Extract subgraph rooted at specific terms
sub_onto = subgraph(onto, ["GO:0008150"])
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
- `is_a`: Parent relationships (is_a)
- `synonym`: Alternative names
- `xref`: Cross-references to other databases
- `subset`: GO subsets/slims
- Relationship types: `part_of`, `regulates`, `has_part`, and other custom relationships

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

### With Phenotype Module
```python
from metainformant.phenotype import load_antwiki_json
from metainformant.ontology import load_go_obo, ancestors

# Functional annotation of phenotypic traits
phenotype_data = load_antwiki_json(Path("antwiki_species.json"))

# Load GO for trait functional annotation
go_onto = load_go_obo("go-basic.obo")

# Map traits to GO terms and find broader categories
trait_term = "GO:0008150"  # biological_process
broader_terms = ancestors(go_onto, trait_term)
# Use GO hierarchy for trait categorization
```

## Serialization

Save and load ontologies to/from JSON format:

```python
from metainformant.ontology import save_ontology, load_ontology

# Save ontology
onto = load_go_obo("go.obo")
save_ontology(onto, "output/go_saved.json")

# Load ontology
loaded_onto = load_ontology("output/go_saved.json")
```

## Performance Features

- Efficient ontology traversal algorithms with BFS
- In-memory caching for expensive operations (ancestors/descendants)
- Memory-optimized data structures
- Optional caching controls

**Caching:**
```python
from metainformant.ontology import set_cache_enabled, set_cache_ttl, clear_cache

# Enable/disable caching
set_cache_enabled(True)

# Set cache TTL (default: 3600 seconds)
set_cache_ttl(7200)  # 2 hours

# Clear cache
clear_cache()

# Disable caching for individual calls
ancestors_set = ancestors(onto, "GO:0008150", use_cache=False)
```

## Error Handling

All functions include error handling:

```python
from metainformant.ontology import load_go_obo, ancestors
from metainformant.core.errors import ValidationError, IOError

try:
    onto = load_go_obo("nonexistent.obo")
except IOError as e:
    print(f"File error: {e}")

try:
    ancestors_set = ancestors(onto, "INVALID_TERM")
except ValueError as e:
    print(f"Invalid term: {e}")
```

## Testing

Comprehensive tests cover:
- Ontology parsing accuracy with various OBO formats
- Hierarchy traversal correctness
- Error handling and validation
- Serialization round-trip integrity
- Edge cases (empty ontologies, cycles, orphaned terms)

Run tests:
```bash
python -m pytest tests/test_ontology_*.py -v
```

## Advanced Features (Require External Dependencies)

Enrichment analysis and semantic similarity functions are available but require additional dependencies:

```python
from metainformant.ontology import enrich_genes, semantic_similarity

# Requires scipy: pip install scipy
# Requires gene-to-term mappings (GAF/GPAD format)

gene_to_terms = {"GENE1": {"GO:0008150"}, "GENE2": {"GO:0008150"}}
genes = ["GENE1", "GENE2"]

# Enrichment analysis (placeholder - requires scipy)
# result = enrich_genes(genes, None, onto, gene_to_terms)

# Semantic similarity (placeholder - requires scipy and annotations)
# similarity = semantic_similarity(onto, "GO:0008150", "GO:0009987")
```

**Note**: These functions are placeholders that raise `ImportError` if `scipy` is not installed. Full implementation requires:
- `scipy` package for statistical tests
- GAF/GPAD file parsers for gene annotations
- Information content calculations

## Limitations

- Enrichment analysis and semantic similarity are placeholder implementations requiring scipy
- Complex OBO qualifiers and advanced features may require specialized libraries
- Caching is in-memory only (not persistent across sessions)

## Dependencies

- Core utilities: `metainformant.core` (io, logging, paths, errors, validation)
- No external dependencies required for basic functionality

## See Also

- **[AGENTS.md](AGENTS.md)**: AI agent contributions and development details for the ontology module

## Related Modules

The Ontology module integrates with several other METAINFORMANT modules:

- **Protein Module**: Gene Ontology (GO) term assignment to protein sequences; functional annotation and enrichment analysis
- **Networks Module**: Ontology-based network construction, functional enrichment of network modules, and semantic similarity networks
- **Information Module**: Information content calculations, semantic similarity measures, and ontology complexity analysis
- **RNA Module**: Functional annotation of genes from RNA-seq data; pathway enrichment analysis
- **GWAS Module**: Functional enrichment of GWAS results; gene set analysis using ontology hierarchies
- **Multi-omics Module**: Cross-omics functional integration using ontology-based annotations
- **Quality Module**: Validation of ontology annotations and enrichment analysis quality control
- **Visualization Module**: Ontology hierarchy visualization, enrichment plot generation, and semantic similarity heatmaps
- **ML Module**: Ontology-based feature engineering and functional prediction models
- **Epigenome Module**: Functional annotation of epigenetic modifications and regulatory elements

This module provides essential tools for ontology analysis and functional annotation in biological research.
