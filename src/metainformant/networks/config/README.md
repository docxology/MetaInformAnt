# Network Configuration

Typed dataclass configuration for network analysis, community detection, pathway enrichment, PPI analysis, GRN inference, and workflow orchestration.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `config` module |
| `config.py` | Dataclass configurations for all network analysis components |

## Key Classes

| Class | Description |
|-------|-------------|
| `NetworkConfig` | Base configuration (directed, weighted, min edge weight, max nodes) |
| `CommunityDetectionConfig` | Community detection method, resolution, number of communities |
| `PathwayEnrichmentConfig` | Enrichment method, correction, overlap and p-value thresholds |
| `PPIConfig` | Protein-protein interaction analysis settings |
| `GRNConfig` | Gene regulatory network inference parameters |
| `NetworkWorkflowConfig` | High-level workflow orchestration configuration |

## Usage

```python
from metainformant.networks.config.config import (
    NetworkConfig,
    CommunityDetectionConfig,
    NetworkWorkflowConfig,
)

net_cfg = NetworkConfig(directed=True, weighted=True, min_edge_weight=0.5)
comm_cfg = CommunityDetectionConfig(method="greedy", resolution=1.0)
wf_cfg = NetworkWorkflowConfig()
```
