# Ecology Module

The `ecology` module provides tools for ecological metadata management, community analysis, and environmental data integration.

## Overview

This module handles ecological data including species diversity, community composition, and environmental parameters that influence biological systems.

### Module Architecture

```mermaid
graph TB
    subgraph "Ecology Module"
        CommunitycommunityCommunityAnalysis[community_Community Analysis]
        EnvironmentalenvironmentalEnvironmentalData[environmental_Environmental Data]
        InteractionsinteractionsSpeciesInteractions[interactions_Species Interactions]
        WorkflowworkflowWorkflowOrchestration[workflow_Workflow Orchestration]
    end
    
    subgraph "Input Data"
        AbundanceabundanceTables[Abundance Tables]
        EnvDataenvironmentalData[Environmental Data]
        InteractionDatainteractionData[Interaction Data]
    end
    
    subgraph "Other Modules"
        networks[networks]
        dna[dna]
        information[information]
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
    AcommunityData[Community Data] --> BabundanceMatrix[Abundance Matrix]
    B --> CspeciesRichness[Species Richness]

    C --> D{Index Type}
    D -->|Alpha Diversity| EshannonIndex[Shannon Index]
    D -->|Species Richness| FobservedSpecies[Observed Species]
    D -->|Evenness| GsimpsonIndex[Simpson Index]
    D -->|Dominance| H[Berger-Parker]

    E --> IdiversityMetrics[Diversity Metrics]
    F --> I
    G --> I
    H --> I

    I --> JrarefactionCurves[Rarefaction Curves]
    J --> KsampleSizeStandardization[Sample Size Standardization]

    K --> LstatisticalComparison[Statistical Comparison]
    L --> MsignificanceTesting[Significance Testing]

    M --> NdiversityResults[Diversity Results]
    N --> OcommunityInterpretation[Community Interpretation]


    subgraph "Data Types"
        PcountData[Count Data] -.-> B
        Q[Presence/Absence] -.-> B
        RbiomassData[Biomass Data] -.-> B
        ScoverData[Cover Data] -.-> B
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
    AmultipleCommunities[Multiple Communities] --> BdistanceMatrix[Distance Matrix]
    B --> C{Dissimilarity Measure}

    C -->|Bray-Curtis| D[Abundance-based]
    C -->|Jaccard| E[Presence/Absence]
    C -->|Unifrac| F[Phylogenetic]
    C -->|Sorensen| GquantitativeSorensen[Quantitative Sorensen]

    D --> HdissimilarityValues[Dissimilarity Values]
    E --> H
    F --> H
    G --> H

    H --> IordinationAnalysis[Ordination Analysis]
    I --> J[NMDS]
    I --> K[PCoA]
    I --> L[MDS]

    J --> M[Visualization]
    K --> M
    L --> M

    M --> NclusterAnalysis[Cluster Analysis]
    N --> OhierarchicalClustering[Hierarchical Clustering]
    N --> Pk-meansClustering[K-means Clustering]

    O --> QcommunityGroups[Community Groups]
    P --> Q

    Q --> RindicatorSpecies[Indicator Species]
    R --> ScharacteristicTaxa[Characteristic Taxa]

    S --> TcommunityClassification[Community Classification]


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
    AspeciesData[Species Data] --> BenvironmentalData[Environmental Data]
    B --> CvariableSelection[Variable Selection]

    A --> DcommunityMatrix[Community Matrix]
    D --> E[Transformation]

    C --> FcollinearityCheck[Collinearity Check]
    F --> GvariableReduction[Variable Reduction]

    E --> HstatisticalModeling[Statistical Modeling]
    G --> H

    H --> I{Method}
    I -->|CCA| JcanonicalCorrelation[Canonical Correlation]
    I -->|RDA| KredundancyAnalysis[Redundancy Analysis]
    I -->|PLS| LpartialLeastSquares[Partial Least Squares]
    I -->|GLM| MgeneralizedLinearModels[Generalized Linear Models]

    J --> Nenvironment-communityRelations[Environment-Community Relations]
    K --> N
    L --> N
    M --> N

    N --> OsignificanceTesting[Significance Testing]
    O --> PpermutationTests[Permutation Tests]
    O --> Q[ANOVA]

    P --> RsignificantRelationships[Significant Relationships]
    Q --> R

    R --> SvariableImportance[Variable Importance]
    S --> TkeyEnvironmentalDrivers[Key Environmental Drivers]

    T --> UmanagementImplications[Management Implications]


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
    AspeciesInteractionData[Species Interaction Data] --> BnetworkConstruction[Network Construction]
    B --> CinteractionMatrix[Interaction Matrix]

    C --> DnetworkProperties[Network Properties]
    D --> E[Connectivity]
    D --> F[Modularity]
    D --> GcentralityMeasures[Centrality Measures]

    E --> HnetworkMetrics[Network Metrics]
    F --> H
    G --> H

    H --> IcommunityDetection[Community Detection]
    I --> JmodulesIdentification[Modules Identification]

    J --> KkeystoneSpecies[Keystone Species]
    K --> LnetworkHubs[Network Hubs]

    L --> MtrophicLevels[Trophic Levels]
    M --> NfoodWebAnalysis[Food Web Analysis]

    N --> OstabilityAnalysis[Stability Analysis]
    O --> PresilienceMeasures[Resilience Measures]

    P --> QconservationPriorities[Conservation Priorities]
    Q --> RmanagementStrategies[Management Strategies]


    subgraph "Interaction Types"
        S[Predation] -.-> C
        T[Competition] -.-> C
        U[Mutualism] -.-> C
        V[Parasitism] -.-> C
    end

    subgraph "Network Metrics"
        WdegreeDistribution[Degree Distribution] -.-> H
        XclusteringCoefficient[Clustering Coefficient] -.-> H
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
    AtemporalCommunityData[Temporal Community Data] --> BtimeSeriesAnalysis[Time Series Analysis]
    B --> CsuccessionPatterns[Succession Patterns]

    C --> DcommunityTrajectories[Community Trajectories]
    D --> EdirectionalChanges[Directional Changes]
    E --> FconvergenceAnalysis[Convergence Analysis]

    F --> GsuccessionalStages[Successional Stages]
    G --> HtransitionProbabilities[Transition Probabilities]

    H --> ImarkovChainModels[Markov Chain Models]
    I --> JstableStates[Stable States]

    J --> KdisturbanceResponse[Disturbance Response]
    K --> LrecoveryTrajectories[Recovery Trajectories]

    L --> MresilienceMetrics[Resilience Metrics]
    M --> NvulnerabilityAssessment[Vulnerability Assessment]

    N --> OconservationStrategies[Conservation Strategies]
    O --> PadaptiveManagement[Adaptive Management]


    subgraph "Temporal Patterns"
        QseasonalChanges[Seasonal Changes] -.-> C
        R[Succession] -.-> C
        ScyclicDynamics[Cyclic Dynamics] -.-> C
        TirregularChanges[Irregular Changes] -.-> C
    end

    subgraph "Analytical Methods"
        UtimeSeriesDecomposition[Time Series Decomposition] -.-> B
        VchangePointDetection[Change Point Detection] -.-> E
        WtrendAnalysis[Trend Analysis] -.-> F
    end

    subgraph "Management Applications"
        X[Restoration] -.-> O
        Y[Monitoring] -.-> O
        ZimpactAssessment[Impact Assessment] -.-> O
    end
```

### Biodiversity Assessment Framework

```mermaid
graph TD
    AbiodiversityData[Biodiversity Data] --> BscaleDefinition[Scale Definition]
    B --> C{Assessment Level}

    C -->|Alpha| DlocalDiversity[Local Diversity]
    C -->|Beta| E[Between-Habitat]
    C -->|Gamma| FregionalDiversity[Regional Diversity]

    D --> GdiversityIndices[Diversity Indices]
    E --> HturnoverMeasures[Turnover Measures]
    F --> IspeciesPool[Species Pool]

    G --> JrichnessEstimation[Richness Estimation]
    H --> KdissimilarityMetrics[Dissimilarity Metrics]
    I --> Lspecies-areaRelationships[Species-Area Relationships]

    J --> M[Rarefaction]
    K --> N[Nestedness]
    L --> OendemismAnalysis[Endemism Analysis]

    M --> PsamplingCompleteness[Sampling Completeness]
    N --> QcommunityStructure[Community Structure]
    O --> RconservationPriority[Conservation Priority]

    P --> SassessmentResults[Assessment Results]
    Q --> S
    R --> S

    S --> TbiodiversityStatus[Biodiversity Status]
    T --> UmanagementRecommendations[Management Recommendations]


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
