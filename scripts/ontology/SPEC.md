# SPEC: Ontology Scripts

Orchestration of functional annotation and semantic similarity analysis.

## Workflows

- `enrich_gene_annotations.py`: Performs GO enrichment analysis on gene lists.
- `calculate_semantic_similarity.py`: Quantitative comparison of functional profiles using ontology hierarchies.

## Standards

- **Reference Data**: Uses official OBO files from Gene Ontology or similar repositories.
- **Reporting**: Exports enrichment results with p-values and FDR-corrected statistics.
