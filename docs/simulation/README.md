# SIMULATION

## Overview
Evolutionary and Population Genetics Simulation module for METAINFORMANT.

## 📦 Contents
- **[models/](models/)**
- **[visualization/](visualization/)**
- **[workflow/](workflow/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Simulation Module"
        M[models/] --> |sequences.py| SQ[Sequence Evolution]
        M --> |popgen.py| PG[Population Genetics]
        M --> |rna.py| RN[RNA-seq Counts]
        M --> |agents.py| AG[Agent-Based Ecosystems]

        W[workflow/] --> |workflow.py| WF[SimulationConfig & Orchestration]
        V[visualization/] --> |visualization.py| VZ[Simulation Plots]
        B[benchmark/] --> BM[Performance Benchmarks]
    end

    WF --> SQ
    WF --> PG
    WF --> RN
    WF --> AG
    SQ --> VZ
    PG --> VZ
    AG --> VZ
```

## Usage
Import module:
```python
from metainformant.simulation import ...
```
