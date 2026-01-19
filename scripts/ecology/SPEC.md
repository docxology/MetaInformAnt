# SPEC: Ecology Scripts

Scripts for community diversity and biodiversity analysis.

## Core Workflows

- `calculate_biodiversity.py`: Computes Shannon, Simpson, and other diversity indices for provided species matrices.
- `analyze_community_composition.py`: Performs ordination and clustering on ecological data.

## Standards

- **Data Formats**: Supports CSV and Parquet inputs via `metainformant.core.io`.
- **Validation**: Ensures species matrices are non-negative and correctly formatted.
