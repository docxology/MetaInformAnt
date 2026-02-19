# ECOLOGY

## Overview
Ecology and community analysis module for METAINFORMANT.

## 📦 Contents
- **[analysis/](analysis/)**
- **[visualization/](visualization/)**
- `[__init__.py](__init__.py)`

## 📊 Structure

```mermaid
graph TD
    subgraph "Ecology Module"
        AC[analysis/community.py] --> |calculate_diversity| DIV[Shannon, Simpson, Richness]

        AO[analysis/ordination.py] --> |pcoa| PCO[Principal Coordinates]
        AO --> |nmds| NM[Non-metric MDS]
        AO --> |cca| CC[Canonical Correspondence]

        AI[analysis/indicators.py] --> |indval| IV[Indicator Species]
        AI --> |permanova| PA[PERMANOVA]
        AI --> |anosim| AN[ANOSIM]
        AI --> |simper| SP[SIMPER]

        AF[analysis/functional.py] --> |functional_richness| FR[FRic, FEve, FDiv]
        AF --> |raos_quadratic_entropy| RQ[Rao's Q]

        AM[analysis/macroecology.py] --> |fit_logseries| SAD[SAD Models]
        AM --> |species_area_power| SAR[Species-Area]

        PH[phylogenetic/diversity.py] --> PD[Phylogenetic Diversity]

        VZ[visualization/] --> VP[Community Plots]
    end
```

## Usage
Import module:
```python
from metainformant.ecology import ...
```
