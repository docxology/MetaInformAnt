# Ontology Analysis Scripts

Functional annotation and ontology analysis workflow orchestrators.

## Directory Structure

```
scripts/ontology/
├── run_ontology_analysis.py       # Ontology analysis workflow orchestrator
└── README.md                      # This file
```

## Ontology Analysis Workflow (`run_ontology_analysis.py`)

Comprehensive ontology analysis workflow orchestrator for Gene Ontology (GO) and other biological ontologies, including term queries, enrichment analysis, and semantic similarity calculations.

**Features:**
- Ontology loading and parsing (OBO format)
- Term querying and navigation
- Semantic similarity calculations
- Ontology enrichment analysis
- Subgraph extraction and analysis

**Usage:**
```bash
# Load GO and generate summary
python3 scripts/ontology/run_ontology_analysis.py --go go.obo --write-summary --output output/ontology/go_summary

# Query specific term
python3 scripts/ontology/run_ontology_analysis.py --go go.obo --query-term GO:0008150 --ancestors --descendants

# Extract subgraph
python3 scripts/ontology/run_ontology_analysis.py --go go.obo --subgraph GO:0008150,GO:0003674 --output output/ontology/subgraph
```

**Options:**
- `--go`: Input Gene Ontology OBO file
- `--query-term`: Query specific ontology term
- `--ancestors`: Include ancestor terms in analysis
- `--descendants`: Include descendant terms in analysis
- `--subgraph`: Extract subgraph for specified terms
- `--write-summary`: Generate ontology summary statistics
- `--output`: Output directory (defaults to output/ontology/)
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/ontology/
├── ontology_data/                 # Parsed ontology structures
│   ├── go_terms.json
│   ├── go_relationships.json
│   └── ontology_metadata.json
├── term_queries/                  # Term-specific analysis results
│   ├── term_details.json
│   ├── ancestors.json
│   ├── descendants.json
│   └── related_terms.json
├── subgraphs/                     # Extracted ontology subgraphs
│   ├── subgraph_terms.json
│   ├── subgraph_edges.json
│   └── subgraph_visualization.png
├── semantic_similarity/           # Similarity analysis results
│   ├── term_similarities.json
│   └── similarity_matrix.json
├── enrichment_analysis/           # Enrichment analysis results
│   ├── enrichment_results.json
│   └── enrichment_plots/
├── ontology_plots/                # Generated visualizations
│   ├── ontology_structure.png
│   ├── term_hierarchy.png
│   └── enrichment_barplot.png
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.ontology**: Core ontology analysis functionality
- **GOATOOLS/pronto**: Ontology parsing and analysis
- **Core utilities**: I/O, logging, path management
- **Visualization**: Ontology plotting and graph visualization

## Dependencies

- **metainformant.ontology**: Ontology analysis module
- **pronto**: Ontology parsing
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Ontology Analysis Documentation](../../docs/ontology/README.md)
- [Gene Ontology Guide](../../docs/ontology/go_guide.md)
- [Semantic Similarity](../../docs/ontology/semantic_similarity.md)
- [METAINFORMANT CLI](../../docs/cli.md)
