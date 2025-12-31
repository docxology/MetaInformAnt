# Ecology Module

The `ecology` module provides tools for ecological metadata management, community analysis, and environmental data integration.

## Overview

This module handles ecological data including species diversity, community composition, and environmental parameters that influence biological systems.

### Module Architecture

```mermaid
graph TB
    subgraph "Ecology Module"
        Community[community<br/>Community Analysis]
        Environmental[environmental<br/>Environmental Data]
        Interactions[interactions<br/>Species Interactions]
        Workflow[workflow<br/>Workflow Orchestration]
    end
    
    subgraph "Input Data"
        Abundance[Abundance Tables]
        EnvData[Environmental Data]
        InteractionData[Interaction Data]
    end
    
    subgraph "Other Modules"
        Networks_Mod[networks]
        DNA_Mod[dna]
        Info_Mod[information]
    end
    
    Abundance --> Community
    EnvData --> Environmental
    InteractionData --> Interactions
    Community --> Workflow
    Environmental --> Workflow
    Interactions --> Workflow
    Workflow --> Networks_Mod
    Workflow --> DNA_Mod
    Workflow --> Info_Mod
```

### Community Diversity Analysis

```mermaid
graph TD
    A[Community Data] --> B[Abundance Matrix]
    B --> C[Species Richness]

    C --> D{Index Type}
    D -->|Alpha Diversity| E[Shannon Index]
    D -->|Species Richness| F[Observed Species]
    D -->|Evenness| G[Simpson Index]
    D -->|Dominance| H[Berger-Parker]

    E --> I[Diversity Metrics]
    F --> I
    G --> I
    H --> I

    I --> J[Rarefaction Curves]
    J --> K[Sample Size Standardization]

    K --> L[Statistical Comparison]
    L --> M[Significance Testing]

    M --> N[Diversity Results]
    N --> O[Community Interpretation]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style I fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style O fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Data Types"
        P[Count Data] -.-> B
        Q[Presence/Absence] -.-> B
        R[Biomass Data] -.-> B
        S[Cover Data] -.-> B
    end

    subgraph "Diversity Measures"
        T[Richness] -.-> C
        U[Evenness] -.-> I
        V[Dominance] -.-> I
    end

    subgraph "Statistical Tests"
        W[t-test] -.-> M
        X[ANOVA] -.-> M
        Y[PERMANOVA] -.-> M
    end
```

### Beta Diversity and Community Comparison

```mermaid
graph TD
    A[Multiple Communities] --> B[Distance Matrix]
    B --> C{Dissimilarity Measure}

    C -->|Bray-Curtis| D[Abundance-based]
    C -->|Jaccard| E[Presence/Absence]
    C -->|Unifrac| F[Phylogenetic]
    C -->|Sorensen| G[Quantitative Sorensen]

    D --> H[Dissimilarity Values]
    E --> H
    F --> H
    G --> H

    H --> I[Ordination Analysis]
    I --> J[NMDS]
    I --> K[PCoA]
    I --> L[MDS]

    J --> M[Visualization]
    K --> M
    L --> M

    M --> N[Cluster Analysis]
    N --> O[Hierarchical Clustering]
    N --> P[K-means Clustering]

    O --> Q[Community Groups]
    P --> Q

    Q --> R[Indicator Species]
    R --> S[Characteristic Taxa]

    S --> T[Community Classification]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style H fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style T fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Distance Measures"
        U[Euclidean] -.-> B
        V[Manhattan] -.-> B
        W[Canberra] -.-> B
    end

    subgraph "Ordination Methods"
        X[PCA] -.-> I
        Y[CA] -.-> I
        Z[DCA] -.-> I
    end

    subgraph "Classification"
        AA[Supervised] -.-> Q
        BB[Unsupervised] -.-> Q
        CC[Fuzzy] -.-> Q
    end
```

### Species-Environment Relationships

```mermaid
graph TD
    A[Species Data] --> B[Environmental Data]
    B --> C[Variable Selection]

    A --> D[Community Matrix]
    D --> E[Transformation]

    C --> F[Collinearity Check]
    F --> G[Variable Reduction]

    E --> H[Statistical Modeling]
    G --> H

    H --> I{Method}
    I -->|CCA| J[Canonical Correlation]
    I -->|RDA| K[Redundancy Analysis]
    I -->|PLS| L[Partial Least Squares]
    I -->|GLM| M[Generalized Linear Models]

    J --> N[Environment-Community Relations]
    K --> N
    L --> N
    M --> N

    N --> O[Significance Testing]
    O --> P[Permutation Tests]
    O --> Q[ANOVA]

    P --> R[Significant Relationships]
    Q --> R

    R --> S[Variable Importance]
    S --> T[Key Environmental Drivers]

    T --> U[Management Implications]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style H fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Environmental Variables"
        V[Climate] -.-> B
        W[Soil] -.-> B
        X[Topography] -.-> B
        Y[Disturbance] -.-> B
    end

    subgraph "Data Transformations"
        Z[Hellinger] -.-> E
        AA[Chord] -.-> E
        BB[Log] -.-> E
    end

    subgraph "Model Evaluation"
        CC[R²] -.-> O
        DD[AIC] -.-> O
        EE[BIC] -.-> O
    end
```

### Ecological Network Analysis

```mermaid
graph TD
    A[Species Interaction Data] --> B[Network Construction]
    B --> C[Interaction Matrix]

    C --> D[Network Properties]
    D --> E[Connectivity]
    D --> F[Modularity]
    D --> G[Centrality Measures]

    E --> H[Network Metrics]
    F --> H
    G --> H

    H --> I[Community Detection]
    I --> J[Modules Identification]

    J --> K[Keystone Species]
    K --> L[Network Hubs]

    L --> M[Trophic Levels]
    M --> N[Food Web Analysis]

    N --> O[Stability Analysis]
    O --> P[Resilience Measures]

    P --> Q[Conservation Priorities]
    Q --> R[Management Strategies]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style H fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style R fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Interaction Types"
        S[Predation] -.-> C
        T[Competition] -.-> C
        U[Mutualism] -.-> C
        V[Parasitism] -.-> C
    end

    subgraph "Network Metrics"
        W[Degree Distribution] -.-> H
        X[Clustering Coefficient] -.-> H
        Y[Betweenness] -.-> H
        Z[Closeness] -.-> H
    end

    subgraph "Stability Measures"
        AA[Robustness] -.-> P
        BB[Invariance] -.-> P
        CC[Resistance] -.-> P
    end
```

### Temporal Community Dynamics

```mermaid
graph TD
    A[Temporal Community Data] --> B[Time Series Analysis]
    B --> C[Succession Patterns]

    C --> D[Community Trajectories]
    D --> E[Directional Changes]
    E --> F[Convergence Analysis]

    F --> G[Successional Stages]
    G --> H[Transition Probabilities]

    H --> I[Markov Chain Models]
    I --> J[Stable States]

    J --> K[Disturbance Response]
    K --> L[Recovery Trajectories]

    L --> M[Resilience Metrics]
    M --> N[Vulnerability Assessment]

    N --> O[Conservation Strategies]
    O --> P[Adaptive Management]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style F fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style P fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Temporal Patterns"
        Q[Seasonal Changes] -.-> C
        R[Succession] -.-> C
        S[Cyclic Dynamics] -.-> C
        T[Irregular Changes] -.-> C
    end

    subgraph "Analytical Methods"
        U[Time Series Decomposition] -.-> B
        V[Change Point Detection] -.-> E
        W[Trend Analysis] -.-> F
    end

    subgraph "Management Applications"
        X[Restoration] -.-> O
        Y[Monitoring] -.-> O
        Z[Impact Assessment] -.-> O
    end
```

### Biodiversity Assessment Framework

```mermaid
graph TD
    A[Biodiversity Data] --> B[Scale Definition]
    B --> C{Assessment Level}

    C -->|Alpha| D[Local Diversity]
    C -->|Beta| E[Between-Habitat]
    C -->|Gamma| F[Regional Diversity]

    D --> G[Diversity Indices]
    E --> H[Turnover Measures]
    F --> I[Species Pool]

    G --> J[Richness Estimation]
    H --> K[Dissimilarity Metrics]
    I --> L[Species-Area Relationships]

    J --> M[Rarefaction]
    K --> N[Nestedness]
    L --> O[Endemism Analysis]

    M --> P[Sampling Completeness]
    N --> Q[Community Structure]
    O --> R[Conservation Priority]

    P --> S[Assessment Results]
    Q --> S
    R --> S

    S --> T[Biodiversity Status]
    T --> U[Management Recommendations]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style J fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Taxonomic Groups"
        V[Plants] -.-> A
        W[Animals] -.-> A
        X[Microorganisms] -.-> A
        Y[Fungi] -.-> A
    end

    subgraph "Spatial Scales"
        Z[Plot] -.-> D
        AA[Habitat] -.-> E
        BB[Landscape] -.-> F
        CC[Biome] -.-> F
    end

    subgraph "Conservation Metrics"
        DD[Hotspots] -.-> R
        EE[Gaps] -.-> R
        FF[Threats] -.-> R
    end
```

## Submodules

### Community Analysis (`community.py`)
Tools for analyzing ecological communities and species interactions.

**Key Features:**
- Species diversity metrics (Shannon, Simpson indices)
- Community composition analysis
- Species abundance distributions
- Beta diversity calculations

**Usage:**
```python
from metainformant.ecology.community import (
    shannon_diversity,
    simpson_diversity,
    species_richness,
    pielou_evenness,
    chao1_estimator
)

# Community analysis
abundance_data = [10, 8, 6, 4, 2, 1]
diversity = shannon_diversity(abundance_data)
simpson = simpson_diversity(abundance_data)
richness = species_richness(abundance_data)
evenness = pielou_evenness(abundance_data)
richness_est = chao1_estimator([int(x) for x in abundance_data])
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import population
from metainformant.ecology.community import shannon_diversity

# Ecological context for genetic diversity
genetic_diversity = population.nucleotide_diversity(sequences)
# Ecological diversity from species abundances
species_abundances = [10, 8, 6, 4, 2]
ecological_diversity = shannon_diversity(species_abundances)
```

### With Visualization Module
```python
from metainformant.ecology.community import shannon_diversity, simpson_diversity
from metainformant.visualization import lineplot, scatter_plot
import numpy as np

# Visualize diversity metrics across sites
sites = ["Site1", "Site2", "Site3", "Site4"]
abundance_matrices = [
    [10, 8, 6, 4, 2],
    [5, 7, 9, 3, 1],
    [12, 4, 8, 6, 3],
    [8, 6, 5, 4, 2]
]

# Calculate diversity for each site
shannon_values = [shannon_diversity(abund) for abund in abundance_matrices]
simpson_values = [simpson_diversity(abund) for abund in abundance_matrices]

# Visualize diversity patterns
ax = scatter_plot(shannon_values, simpson_values, 
                  xlabel="Shannon Diversity", ylabel="Simpson Diversity",
                  title="Diversity Metrics Across Sites")
```

### With Information Theory Module
```python
from metainformant.ecology.community import shannon_diversity
from metainformant.information import shannon_entropy

# Compare ecological diversity with information-theoretic entropy
species_abundances = [10, 8, 6, 4, 2, 1]

# Ecological diversity (normalized by total abundance)
eco_diversity = shannon_diversity(species_abundances)

# Information-theoretic entropy (from relative abundances)
total = sum(species_abundances)
proportions = [a / total for a in species_abundances]
info_entropy = shannon_entropy(proportions)

# Both metrics measure diversity/information content
```

## Data Sources

- Species occurrence databases
- Environmental monitoring data
- Biodiversity surveys and inventories
- Ecological metadata repositories

## Performance Features

- Efficient processing of large ecological datasets
- Memory-optimized community calculations
- Support for sparse ecological matrices

## Testing

Comprehensive tests cover:
- Diversity metric calculations
- Community analysis algorithms
- Integration with environmental data

## Dependencies

- NumPy for numerical computations
- Pandas for data manipulation
- Optional: specialized ecological packages

## See Also

- **[AGENTS.md](AGENTS.md)**: AI agent contributions and development details for the ecology module

## Related Modules

The Ecology module integrates with several other METAINFORMANT modules:

- **Networks Module**: Ecological interaction networks, food web analysis, and community structure modeling
- **DNA Module**: DNA barcoding for species identification; biodiversity assessment using sequence data
- **Information Module**: Information-theoretic analysis of ecological complexity and biodiversity patterns
- **Multi-omics Module**: Ecological multi-omics data integration; environmental metagenomics
- **Visualization Module**: Ecological data visualization, diversity plots, and community structure diagrams
- **ML Module**: Machine learning analysis of ecological data; species distribution modeling
- **Quality Module**: Ecological data validation and quality control; biodiversity assessment standards
- **Simulation Module**: Ecosystem simulation and modeling; virtual ecological experiments
- **Math Module**: Mathematical modeling of population dynamics and ecological processes
- **Phenotype Module**: Ecological trait analysis and morphological biodiversity studies

This module provides essential tools for ecological data analysis and biodiversity research.
