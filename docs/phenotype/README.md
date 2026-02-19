# PHENOTYPE

## Overview
Phenotype module for MetaInformAnt.

## 📦 Contents
- **[analysis/](analysis/)**
- **[behavior/](behavior/)**
- **[chemical/](chemical/)**
- **[data/](data/)**
- **[electronic/](electronic/)**
- **[morphological/](morphological/)**
- **[sonic/](sonic/)**
- **[visualization/](visualization/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Phenotype Module"
        D[data/] --> |scraper.py| SC[AntWiki Web Scraping]
        D --> |antwiki.py| AW[AntWiki Record Loading]

        A[analysis/] --> |life_course.py| LC[Life Course Trajectories]

        B[behavior/] --> BH[Behavioral Phenotypes]
        C[chemical/] --> CH[Chemical Profiles]
        MO[morphological/] --> MR[Morphometric Measurements]
        SO[sonic/] --> SN[Acoustic Phenotypes]
        EL[electronic/] --> EE[Electronic Phenotyping]

        V[visualization/] --> |visualization.py| VZ[Trait Plots & Networks]
        W[workflow/] --> |pipeline.py| PL[PhenotypePipeline]

        GI[gwas_integration/] --> GW[GWAS-Phenotype Linking]
        IN[integration/] --> IT[Cross-Domain Integration]
    end

    D --> A
    A --> V
    A --> W
    GI --> IT
```

## Usage
Import module:
```python
from metainformant.phenotype import ...
```
