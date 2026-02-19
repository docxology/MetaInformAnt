# ML

## Overview
Machine learning module for METAINFORMANT.

## 📦 Contents
- **[evaluation/](evaluation/)**
- **[features/](features/)**
- **[models/](models/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "ML Module"
        F[features/] --> |features.py| FS[Feature Selection]
        F --> |dimensionality.py| DR[Dimensionality Reduction]

        M[models/] --> |classification.py| CL[Classification Models]
        M --> |regression.py| RG[Regression Models]

        E[evaluation/] --> EV[Model Evaluation & Scoring]
        I[interpretability/] --> IN[SHAP, Feature Importance]
        A[automl/] --> AM[Automated Model Selection]

        L[llm/] --> |ollama/client.py| OC[OllamaClient]
    end

    FS --> CL
    FS --> RG
    DR --> EV
    CL --> EV
    RG --> EV
    EV --> IN
```

## Usage
Import module:
```python
from metainformant.ml import ...
```
