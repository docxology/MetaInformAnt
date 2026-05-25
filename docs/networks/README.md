# NETWORKS

## Overview
Network analysis module for METAINFORMANT.

## Contents
- **[analysis/](graph.md)** — Graph algorithms, community detection, pathway analysis
- **[interaction/](ppi.md)** — Protein-protein and regulatory interactions
- **[regulatory/](regulatory.md)** — Gene regulatory network analysis
- **config/** — Network analysis configuration
- **workflow/** — Network analysis workflows

## Structure

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
```python-snippet
from metainformant.networks import ...
```
