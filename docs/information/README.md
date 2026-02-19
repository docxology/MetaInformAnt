# INFORMATION

## Overview
Information theory analysis module for METAINFORMANT.

## 📦 Contents
- **[integration/](integration/)** — Cross-module information integration
- **[metrics/](metrics/)** — Information metrics (`core/`, `advanced/`, `analysis/` subpackages)
- **[network_info/](network_info/)** — Network information analysis
- **[workflow/](workflow/)** — Information theory workflows

## 📊 Structure

```mermaid
graph TD
    subgraph "Information Module"
        MC[metrics/core/] --> |syntactic.py| SH[Shannon Entropy, MI, KL]
        MC --> |continuous.py| DE[Differential Entropy]
        MC --> |estimation.py| EST[Estimators]

        MA[metrics/advanced/] --> |geometry.py| FG[Fisher-Rao, Hellinger]
        MA --> |channel.py| CH[Channel Capacity]
        MA --> |decomposition.py| DC[PID, Synergy]
        MA --> |semantic.py| SEM[Semantic Similarity]

        MN[metrics/analysis/] --> |analysis.py| PRO[Information Profiles]

        NI[network_info/] --> |information_flow.py| TF[Transfer Entropy, Granger]

        INT[integration/] --> |integration.py| BI[Cross-Omic Integration]

        WF[workflow/] --> |workflows.py| BAT[Batch Pipelines]
    end
```

## Usage
Import module:
```python
from metainformant.information import ...
```
