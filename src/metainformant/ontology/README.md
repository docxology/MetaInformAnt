# Ontology Module

The `ontology` module provides tools for functional annotation and ontology analysis using biological ontologies like Gene Ontology (GO).

## Overview

This module handles ontology parsing, hierarchy traversal, term queries, and serialization. Provides lightweight, efficient tools for working with OBO-format ontologies without requiring external database connections. Includes error handling, validation, caching, and support for multiple relationship types.

### Module Architecture

```mermaid
graph TB
    subgraph "Ontology Module"
        GOgoGeneOntology[go_Gene Ontology]
        OBOoboOboParser[obo_OBO Parser]
        QueryqueryTermQueries[query_Term Queries]
        SerializeserializeSerialization[serialize_Serialization]
        TypestypesTypeSystem[types_Type System]
    end
    
    subgraph "Input"
        OBOFileoboFile[OBO File]
        TermsgoTerms[GO Terms]
    end
    
    subgraph "Other Modules"
        protein[protein]
        networks[networks]
        information[information]
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
    AoboFile[OBO File] --> BparseHeader[Parse Header]
    B --> CextractOntologyMetadata[Extract Ontology Metadata]

    A --> DparseTerms[Parse Terms]
    D --> EtermDefinitions[Term Definitions]
    E --> FtermRelationships[Term Relationships]
    F --> Gsynonyms&Cross-references[Synonyms & Cross-references]

    A --> HparseTypedefs[Parse Typedefs]
    H --> IrelationshipTypes[Relationship Types]
    I --> JrelationshipProperties[Relationship Properties]

    C --> KontologyGraph[Ontology Graph]
    F --> K
    I --> K

    K --> L[Validation]
    L --> M{Valid Structure?}

    M -->|Yes| NbuildNetworkxGraph[Build NetworkX Graph]
    M -->|No| OerrorReporting[Error Reporting]

    N --> PindexTerms[Index Terms]
    P --> QtermLookupTables[Term Lookup Tables]
    Q --> RreadyForQueries[Ready for Queries]


    subgraph "OBO Components"
        SheaderSection[Header Section] -.-> B
        TtermStanzas[Term Stanzas] -.-> D
        UtypedefStanzas[Typedef Stanzas] -.-> H
        VinstanceStanzas[Instance Stanzas] -.-> D
    end

    subgraph "Relationship Types"
        WisA[is_a] -.-> F
        XpartOf[part_of] -.-> F
        Y[regulates] -.-> F
        ZhasPart[has_part] -.-> F
    end
```

### Semantic Similarity Analysis

```mermaid
graph TD
    AtwoGoTerms[Two GO Terms] --> BinformationContent[Information Content]
    A --> CcommonAncestors[Common Ancestors]

    B --> DtermIcCalculation[Term IC Calculation]
    C --> ElowestCommonAncestor[Lowest Common Ancestor]

    D --> FsimilarityMethod[Similarity Method]
    E --> F

    F --> G{Similarity Measure}
    G -->|Resnik| HsharedIc[Shared IC]
    G -->|Lin| IsharedIc/AverageIc[Shared IC / Average IC]
    G -->|Jiang-Conrath| J1-SharedIcDistance[1 - Shared IC Distance]
    G -->|Wang| KsemanticContribution[Semantic Contribution]

    H --> LsimilarityScore[Similarity Score]
    I --> L
    J --> L
    K --> L

    L --> MthresholdComparison[Threshold Comparison]
    M --> N{Similar?}

    N -->|Yes| OfunctionalRelationship[Functional Relationship]
    N -->|No| PdifferentFunctions[Different Functions]

    O --> QannotationTransfer[Annotation Transfer]
    P --> RdistinctFunctions[Distinct Functions]


    subgraph "IC Calculation"
        StermFrequency[Term Frequency] -.-> D
        TcorpusSize[Corpus Size] -.-> D
        UannotationCorpus[Annotation Corpus] -.-> D
    end

    subgraph "Applications"
        VgeneClustering[Gene Clustering] -.-> O
        WfunctionPrediction[Function Prediction] -.-> O
        XorthologDetection[Ortholog Detection] -.-> O
    end
```

### Functional Enrichment Analysis

```mermaid
graph TD
    AgeneList[Gene List] --> BbackgroundGenes[Background Genes]
    B --> CgoAnnotations[GO Annotations]

    A --> DannotatedGenes[Annotated Genes]
    D --> EtermCounts[Term Counts]

    C --> FbackgroundCounts[Background Counts]
    F --> GtermStatistics[Term Statistics]

    E --> H{Enrichment Test}
    H -->|Fisher's Exact| IhypergeometricTest[Hypergeometric Test]
    H -->|Chi-square| Jchi-squareTest[Chi-square Test]
    H -->|Binomial| KbinomialTest[Binomial Test]

    I --> Lp-valueCalculation[P-value Calculation]
    J --> L
    K --> L

    L --> MmultipleTestingCorrection[Multiple Testing Correction]
    M --> N{Method}
    N -->|Bonferroni| Ofamily-wiseError[Family-wise Error]
    N -->|BH FDR| PfalseDiscoveryRate[False Discovery Rate]
    N -->|Holm| Qstep-downBonferroni[Step-down Bonferroni]

    O --> RadjustedP-values[Adjusted P-values]
    P --> R
    Q --> R

    R --> SsignificantTerms[Significant Terms]
    S --> TenrichmentResults[Enrichment Results]

    T --> UpathwayVisualization[Pathway Visualization]
    U --> VbiologicalInsights[Biological Insights]


    subgraph "Contingency Table"
        WannotatedInList[Annotated in List] -.-> E
        XnotInList[Not in List] -.-> E
        YannotatedInBackground[Annotated in Background] -.-> F
        ZnotAnnotated[Not Annotated] -.-> F
    end

    subgraph "GO Aspects"
        AAbiologicalProcess[Biological Process] -.-> C
        BBmolecularFunction[Molecular Function] -.-> C
        CCcellularComponent[Cellular Component] -.-> C
    end
```

### Ontology-Based Clustering

```mermaid
graph TD
    AgeneSet[Gene Set] --> BgoAnnotations[GO Annotations]
    B --> CsemanticSimilarityMatrix[Semantic Similarity Matrix]

    C --> D{Clustering Method}
    D -->|Hierarchical| EagglomerativeClustering[Agglomerative Clustering]
    D -->|K-means| Fk-meansOnSimilarity[K-means on Similarity]
    D -->|Community| GcommunityDetection[Community Detection]

    E --> HdistanceMatrix[Distance Matrix]
    F --> I[Centroid-based]
    G --> JmodularityOptimization[Modularity Optimization]

    H --> KclusterAssignments[Cluster Assignments]
    I --> K
    J --> K

    K --> LclusterValidation[Cluster Validation]
    L --> MsilhouetteScore[Silhouette Score]
    L --> NadjustedRandIndex[Adjusted Rand Index]

    M --> OoptimalClustering[Optimal Clustering]
    N --> O

    O --> PfunctionalInterpretation[Functional Interpretation]
    P --> QclusterEnrichment[Cluster Enrichment]
    Q --> RgoTermSummaries[GO Term Summaries]

    R --> SclusterCharacterization[Cluster Characterization]
    S --> TgeneFunctionGroups[Gene Function Groups]


    subgraph "Similarity Measures"
        U[Resnik] -.-> C
        V[Lin] -.-> C
        W[Jiang-Conrath] -.-> C
        X[Wang] -.-> C
    end

    subgraph "Distance Metrics"
        Y1-Similarity[1 - Similarity] -.-> H
        Z[Euclidean] -.-> H
        AA[Manhattan] -.-> H
    end
```

### Ontology Integration Framework

```mermaid
graph TD
    AmultipleOntologies[Multiple Ontologies] --> BontologyMapping[Ontology Mapping]
    B --> CtermAlignment[Term Alignment]
    C --> Dcross-ontologyRelationships[Cross-ontology Relationships]

    A --> EintegratedGraph[Integrated Graph]
    D --> E

    E --> FunifiedQueries[Unified Queries]
    F --> GtermResolution[Term Resolution]
    G --> HannotationIntegration[Annotation Integration]

    H --> IunifiedEnrichment[Unified Enrichment]
    I --> Jcross-ontologyAnalysis[Cross-ontology Analysis]

    J --> KcomprehensiveAnnotations[Comprehensive Annotations]
    K --> LenhancedBiologicalInsights[Enhanced Biological Insights]


    subgraph "Ontology Sources"
        MgeneOntology[Gene Ontology] -.-> A
        NhumanDisease[Human Disease] -.-> A
        OplantOntology[Plant Ontology] -.-> A
        PcellOntology[Cell Ontology] -.-> A
    end

    subgraph "Mapping Methods"
        QexactMatch[Exact Match] -.-> C
        RsemanticSimilarity[Semantic Similarity] -.-> C
        SmanualCuration[Manual Curation] -.-> C
        TmachineLearning[Machine Learning] -.-> C
    end

    subgraph "Integration Benefits"
        UbroaderCoverage[Broader Coverage] -.-> K
        VcontextualInformation[Contextual Information] -.-> K
        Wcross-domainInsights[Cross-domain Insights] -.-> K
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

# Requires scipy: uv pip install scipy
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
