# NETWORKS

## Overview
Network analysis module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)** — Graph algorithms, community detection, pathway analysis
- **[interaction/](interaction/)** — Protein-protein and regulatory interactions
- **[regulatory/](regulatory/)** — Gene regulatory network analysis
- **[config/](config/)** — Network analysis configuration
- **[workflow/](workflow/)** — Network analysis workflows

## 📊 Structure

```mermaid
graph TD
    subgraph "Networks Module"
        AG[analysis/graph.py] --> |BiologicalNetwork| NET[Network Construction]
        AG --> |centrality_measures| ALG[Graph Algorithms]
        AG --> |network_metrics| MET[Topology Metrics]

        AC[analysis/community.py] --> |louvain_communities| LV[Louvain Detection]
        AC --> |label_propagation| LP[Label Propagation]

        AP[analysis/pathway.py] --> |pathway_enrichment_analysis| PE[Pathway Enrichment]

        IP[interaction/ppi.py] --> |load_ppi_network| PPI[PPI Networks]

        IR[interaction/regulatory.py] --> |GeneRegulatoryNetwork| GRN[Regulatory Networks]
        IR --> |infer_grn| INF[GRN Inference]
    end
```

## Usage
Import module:
```python
from metainformant.networks import ...
```
