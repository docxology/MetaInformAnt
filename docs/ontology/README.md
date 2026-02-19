# ONTOLOGY

## Overview
Gene ontology and functional annotation module for METAINFORMANT.

## 📦 Contents
- `core/` — GO ontology loading, OBO parsing, type definitions (`go.py`, `obo.py`, `types.py`)
- `query/` — Ontology querying and serialization (`query.py`, `serialize.py`)
- `visualization/` — Ontology visualization (`visualization.py`)
- `pathway_enrichment/` — Pathway enrichment analysis

## 📊 Structure

```mermaid
graph TD
    subgraph "Ontology Module"
        C[core/] --> GO[GO Loading & Enrichment]
        C --> OBO[OBO Format Parsing]

        Q[query/] --> QR[Ancestor/Descendant Traversal]
        Q --> S[Save/Load OBO & JSON]

        PE[pathway_enrichment/] --> EN[ORA, GSEA, FDR Correction]

        V[visualization/] --> VZ[DAG Plots, Enrichment Charts]
    end

    OBO --> GO
    GO --> QR
    EN --> QR
    VZ --> GO
```

## Usage
Import module:
```python
from metainformant.ontology import ...
```
