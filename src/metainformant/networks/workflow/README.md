# Network Workflow

High-level workflow orchestration for multi-step network analyses combining graph construction, community detection, pathway enrichment, and PPI analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `workflow` module |
| `workflow.py` | `NetworkWorkflow` class for pipeline orchestration |

## Key Classes

| Class | Description |
|-------|-------------|
| `NetworkWorkflow` | Orchestrates multi-step network analysis pipelines with chainable methods |

## Key Methods

| Method | Description |
|--------|-------------|
| `NetworkWorkflow.build_network()` | Build a network from edge lists or correlation matrices |
| `NetworkWorkflow.detect_communities()` | Run community detection on the built network |
| `NetworkWorkflow.enrich_pathways()` | Perform pathway enrichment analysis |
| `NetworkWorkflow.analyze_ppi()` | Run protein-protein interaction analysis |

## Usage

```python
from metainformant.networks.workflow.workflow import NetworkWorkflow
from metainformant.networks.config.config import NetworkWorkflowConfig

wf = NetworkWorkflow(config=NetworkWorkflowConfig())
wf.build_network(edges=edge_list).detect_communities().enrich_pathways()
results = wf.results
```
