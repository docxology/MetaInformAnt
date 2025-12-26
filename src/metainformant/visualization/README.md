# Visualization Module

The `visualization` module provides unified plotting and animation utilities for creating publication-quality figures and interactive visualizations across all METAINFORMANT domains.

## Overview

This module offers a cohesive API for creating various types of plots, from statistical charts to phylogenetic trees and animations. It integrates multiple plotting backends while providing a consistent interface.

### Module Architecture

```mermaid
graph TB
    subgraph "Visualization Module"
        Basic[basic<br/>Basic Plots]
        Statistical[statistical<br/>Statistical Plots]
        Trees[trees<br/>Phylogenetic Trees]
        Networks[networks<br/>Network Visualization]
        Animation[animation<br/>Animations]
        Style[style<br/>Styling]
    end
    
    subgraph "Input Data"
        Data[Data Arrays]
        Trees_Data[Tree Data]
        Network_Data[Network Data]
        TimeSeries[Time Series]
    end
    
    subgraph "All Modules"
        All[All Domain Modules]
    end
    
    Data --> Basic
    Data --> Statistical
    Trees_Data --> Trees
    Network_Data --> Networks
    TimeSeries --> Animation
    All --> Basic
    All --> Statistical
    All --> Trees
    All --> Networks
    All --> Animation
```

### Visualization Pipeline Framework

```mermaid
graph TD
    A[Biological Data] --> B[Data Processing]
    B --> C[Format Conversion]

    C --> D{Visualization Category}
    D -->|Basic| E[Line/Scatter/Bar Plots]
    D -->|Statistical| F[Distribution/Comparison Plots]
    D -->|Genomic| G[Manhattan/Volcano Plots]
    D -->|Expression| H[Heatmaps/Enrichment Plots]
    D -->|Dimension Reduction| I[PCA/UMAP/t-SNE]
    D -->|Network| J[Graph/Network Visualization]
    D -->|Time Series| K[Time Series Analysis]
    D -->|Multi-dimensional| L[Parallel Coordinates/3D]
    D -->|Quality Control| M[QC Metrics Plots]
    D -->|Information Theory| N[Entropy/MI Networks]

    E --> O[Plot Generation]
    F --> O
    G --> O
    H --> O
    I --> O
    J --> O
    K --> O
    L --> O
    M --> O
    N --> O

    O --> P[Styling & Customization]
    P --> Q[Color Schemes]
    P --> R[Layout & Annotations]
    P --> S[Font & Size Settings]

    Q --> T[Figure Assembly]
    R --> T
    S --> T

    T --> U{Output Format}
    U -->|Static| V[PNG/PDF/SVG]
    U -->|Interactive| W[HTML/JavaScript]
    U -->|Animation| X[GIF/MP4]

    V --> Y[Publication Quality]
    W --> Z[Web Integration]
    X --> AA[Dynamic Visualization]

    Y --> BB[Final Output]
    Z --> BB
    AA --> BB

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style O fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style BB fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Data Types"
        CC[Genomic] -.-> A
        DD[Transcriptomic] -.-> A
        EE[Proteomic] -.-> A
        FF[Phenotypic] -.-> A
        GG[Network] -.-> A
    end

    subgraph "Visualization Categories"
        HH[20+ Specialized Modules] -.-> D
        II[100+ Plot Functions] -.-> O
        JJ[Domain Integrations] -.-> O
    end
```

### GWAS Visualization Suite

```mermaid
graph TD
    A[GWAS Results] --> B[Manhattan Plot]
    B --> C[Genome-wide Significance]

    A --> D[Q-Q Plot]
    D --> E[P-value Distribution]

    A --> F[Regional Plot]
    F --> G[Locus Zoom]

    A --> H[PCA Plot]
    H --> I[Population Structure]

    A --> J[Kinship Heatmap]
    J --> K[Relatedness Matrix]

    C --> L[Visualization Suite]
    E --> L
    G --> L
    I --> L
    K --> L

    L --> M{Additional Analysis}
    M -->|Yes| N[Forest Plot]
    M -->|No| O[Complete Suite]

    N --> P[Meta-analysis]
    P --> O

    O --> Q[Publication Figures]
    Q --> R[Interactive Exploration]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style L fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style R fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Core GWAS Plots"
        S[Manhattan] -.-> B
        T[Q-Q] -.-> D
        U[Regional/LocusZoom] -.-> F
        V[Population Structure] -.-> H
    end

    subgraph "Advanced Plots"
        W[Forest] -.-> N
        X[Locus Compare] -.-> M
        Y[Functional Enrichment] -.-> M
        Z[Replication] -.-> M
    end
```

### Single-Cell Visualization Framework

```mermaid
graph TD
    A[Single-Cell Data] --> B[Quality Control Plots]
    B --> C[Library Size Distribution]

    A --> D[Preprocessing Visualization]
    D --> E[Gene Detection vs UMIs]

    A --> F[Dimensionality Reduction]
    F --> G[PCA/UMAP/t-SNE Plots]

    A --> H[Clustering Results]
    H --> I[Cluster Visualization]

    A --> J[Differential Expression]
    J --> K[Volcano Plots]

    A --> L[Trajectory Analysis]
    L --> M[Pseudotime Plots]

    C --> N[Visualization Dashboard]
    E --> N
    G --> N
    I --> N
    K --> N
    M --> N

    N --> O[Interactive Exploration]
    O --> P[Cell Type Annotation]
    P --> Q[Marker Gene Expression]

    Q --> R[Publication Figures]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style N fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style R fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "QC Visualizations"
        S[Library Size] -.-> B
        T[Gene Detection] -.-> B
        U[MT Content] -.-> B
        V[Complexity] -.-> B
    end

    subgraph "Analysis Plots"
        W[DR Embeddings] -.-> F
        X[Cluster Markers] -.-> H
        Y[DE Volcano] -.-> J
        Z[Trajectory] -.-> L
    end
```

### Information Theory Visualization

```mermaid
graph TD
    A[Information Data] --> B[Entropy Profiles]
    B --> C[Shannon Entropy Plot]

    A --> D[Mutual Information]
    D --> E[MI Heatmap/Network]

    A --> F[Complexity Measures]
    F --> G[Rényi Spectra]

    A --> H[Sequence Analysis]
    H --> I[Information Landscapes]

    A --> J[Network Information]
    J --> K[Information Flow Networks]

    C --> L[Information Visualization Suite]
    E --> L
    G --> L
    I --> L
    K --> L

    L --> M[Statistical Significance]
    M --> N[Information Landscapes]

    N --> O[Publication Quality]
    O --> P[Interactive Exploration]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style L fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style P fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Syntactic Information"
        Q[Shannon Entropy] -.-> B
        R[Joint Entropy] -.-> B
        S[Mutual Information] -.-> D
        T[Transfer Entropy] -.-> D
    end

    subgraph "Semantic Information"
        U[Information Content] -.-> F
        V[Semantic Similarity] -.-> F
        W[Semantic Entropy] -.-> F
    end

    subgraph "Visualization Types"
        X[Heatmaps] -.-> L
        Y[Networks] -.-> L
        Z[Profiles] -.-> L
        AA[Landscapes] -.-> N
    end
```

### Animation and Time Series Framework

```mermaid
graph TD
    A[Time Series Data] --> B[Data Preparation]
    B --> C[Temporal Ordering]

    C --> D{Animation Type}
    D -->|Evolution| E[Progressive Changes]
    D -->|Clustering| F[Cluster Dynamics]
    D -->|Network| G[Network Evolution]
    D -->|Trajectory| H[Cell Trajectories]

    E --> I[Frame Generation]
    F --> I
    G --> I
    H --> I

    I --> J[Animation Parameters]
    J --> K[Frame Rate]
    J --> L[Transition Effects]
    J --> M[Color Schemes]

    K --> N[Animation Assembly]
    L --> N
    M --> N

    N --> O{Output Format}
    O -->|GIF| P[GIF Animation]
    O -->|MP4| Q[Video Format]
    O -->|HTML| R[Interactive Animation]

    P --> S[Static Sharing]
    Q --> S
    R --> T[Web Integration]

    S --> U[Complete Animation]
    T --> U

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style I fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Time Series Types"
        V[Genomic Evolution] -.-> A
        W[Gene Expression] -.-> A
        X[Cell Differentiation] -.-> A
        Y[Population Dynamics] -.-> A
    end

    subgraph "Animation Categories"
        Z[Progressive] -.-> D
        AA[Dynamic] -.-> D
        BB[Temporal] -.-> D
        CC[Spatial] -.-> D
    end
```

### Multi-Panel Figure Layout System

```mermaid
graph TD
    A[Multiple Plots] --> B[Layout Design]
    B --> C[Grid Specification]

    A --> D[Plot Integration]
    D --> E[Shared Axes]
    D --> F[Common Legends]

    C --> G[Figure Assembly]
    E --> G
    F --> G

    G --> H[Styling Consistency]
    H --> I[Font Sizes]
    H --> J[Color Palettes]
    H --> K[Line Styles]

    I --> L[Final Layout]
    J --> L
    K --> L

    L --> M{Export Format}
    M -->|Single Figure| N[Publication Layout]
    M -->|Multi-Panel| O[Complex Visualization]
    M -->|Interactive| P[Dashboard Style]

    N --> Q[Journal Submission]
    O --> Q
    P --> R[Web Presentation]

    Q --> S[Complete Visualization]
    R --> S

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style G fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style S fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Layout Types"
        T[Grid Layout] -.-> B
        U[Free-form] -.-> B
        V[Subplot Mosaic] -.-> B
        W[Nested Layouts] -.-> B
    end

    subgraph "Integration Features"
        X[Shared Legends] -.-> D
        Y[Common Scales] -.-> D
        Z[Coordinated Zoom] -.-> D
        AA[Linked Brushing] -.-> P
    end
```

## Module Organization

The visualization package is organized into category-specific modules for clear organization and maintainability:

### Core Plotting Modules

#### Basic Plots (`basic.py`)
Fundamental plotting functions for simple visualizations.

**Functions:**
- `lineplot`: Simple line plots
- `scatter_plot`: Scatter plots
- `bar_plot`: Bar charts
- `pie_chart`: Pie charts
- `area_plot`: Area plots (filled line plots)
- `step_plot`: Step plots
- `heatmap`: Basic heatmaps

**Usage:**
```python
from metainformant.visualization import lineplot, scatter_plot, bar_plot, heatmap

ax = lineplot(None, [1, 4, 2, 8, 5], label="Data")
ax = scatter_plot([1, 2, 3], [4, 5, 6], xlabel="X", ylabel="Y")
ax = bar_plot(["A", "B", "C"], [10, 20, 15])
```

#### Statistical Plots (`statistical.py`)
Statistical visualization functions for data analysis.

**Functions:**
- `histogram`: Histograms
- `box_plot`: Box plots
- `violin_plot`: Violin plots
- `qq_plot`: Q-Q plots for p-value analysis
- `correlation_heatmap`: Correlation matrices
- `density_plot`: Kernel density estimation
- `ridge_plot`: Overlapping density plots
- `roc_curve`: ROC curves
- `precision_recall_curve`: Precision-recall curves
- `residual_plot`: Regression residuals
- `leverage_plot`: Regression leverage

**Usage:**
```python
from metainformant.visualization import histogram, box_plot, qq_plot
import numpy as np

data = np.random.normal(0, 1, 1000)
ax = histogram(data, bins=30)
ax = box_plot([np.random.normal(0, 1, 100) for _ in range(3)], labels=["A", "B", "C"])
```

#### Genomics Plots (`genomics.py`)
Genomic visualization functions for GWAS and sequence analysis.

**Functions:**
- `manhattan_plot`: Manhattan plots for GWAS
- `volcano_plot`: Volcano plots for differential expression
- `regional_plot`: Regional plots for specific genomic regions
- `circular_manhattan_plot`: Circular genome-wide views
- `chromosome_ideogram`: Chromosome maps with markers
- `coverage_plot`: Sequencing coverage visualization
- `variant_plot`: Variant visualization

**Usage:**
```python
from metainformant.visualization import manhattan_plot, volcano_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'position': range(1000, 10000, 100),
    'pvalue': np.random.uniform(1e-9, 1, 90),
    'chromosome': ['chr1'] * 45 + ['chr2'] * 45
})
data['neg_log10_p'] = -np.log10(data['pvalue'])
ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
```

#### Expression Plots (`expression.py`)
Expression analysis visualization functions.

**Functions:**
- `expression_heatmap`: Expression heatmaps with clustering
- `enrichment_plot`: Pathway/gene set enrichment
- `gene_expression_plot`: Single gene expression visualization
- `differential_expression_plot`: Top differentially expressed genes
- `log_fold_change_plot`: Log fold change distributions

#### Dimensionality Reduction (`dimred.py`)
Visualization functions for dimensionality reduction techniques.

**Functions:**
- `pca_plot`: PCA scatter plots
- `umap_plot`: UMAP visualizations
- `tsne_plot`: t-SNE visualizations
- `pca_scree_plot`: Variance explained plots
- `pca_loadings_plot`: PCA loadings visualization
- `biplot`: PCA biplots (samples and loadings)

**Usage:**
```python
from metainformant.visualization import pca_plot, umap_plot
import pandas as pd
import numpy as np

data = pd.DataFrame({
    'PC1': np.random.normal(0, 1, 100),
    'PC2': np.random.normal(0, 1, 100),
    'group': ['A'] * 50 + ['B'] * 50
})
ax = pca_plot(data, hue='group')
```

#### Network Plots (`networks.py`)
Network visualization functions.

**Functions:**
- `network_plot`: Basic network graphs
- `circular_network_plot`: Circular layouts
- `hierarchical_network_plot`: Hierarchical layouts
- `force_directed_plot`: Force-directed layouts
- `community_network_plot`: Community-based coloring

#### Time Series (`timeseries.py`)
Time series visualization functions.

**Functions:**
- `time_series_plot`: Basic time series plots
- `autocorrelation_plot`: Autocorrelation analysis
- `seasonal_decomposition_plot`: Trend, seasonal, residual decomposition
- `forecast_plot`: Forecasts with confidence intervals
- `trend_plot`: Trend analysis with fitted lines

#### Multi-dimensional (`multidim.py`)
Multi-dimensional visualization functions.

**Functions:**
- `pairplot_dataframe`: Pair plots for DataFrames
- `parallel_coordinates_plot`: Parallel coordinates
- `radar_chart`: Radar charts (spider charts)
- `splom_plot`: Scatter plot matrices
- `scatter_3d`: 3D scatter plots

#### Quality Control (`quality.py`)
Quality control visualization functions.

**Functions:**
- `qc_metrics_plot`: Multiple QC metrics
- `quality_score_plot`: Quality score distributions
- `per_base_quality_plot`: Position-specific quality
- `adapter_content_plot`: Adapter content visualization
- `sequence_length_distribution`: Sequence length distributions

#### Information Theory (`information.py`)
Information theory visualization functions.

**Functions:**
- `entropy_plot`: Entropy across positions
- `mutual_information_plot`: MI matrix heatmaps
- `information_profile_plot`: Comprehensive information profiles
- `renyi_spectrum_plot`: Rényi entropy spectra
- `information_network_plot`: MI-based networks

### Specialized Modules

#### Phylogenetic Trees (`trees.py`)
Phylogenetic tree visualization with multiple layouts.

**Functions:**
- `plot_phylo_tree`: Standard tree plots
- `circular_tree_plot`: Circular layouts
- `unrooted_tree_plot`: Unrooted trees
- `tree_comparison_plot`: Side-by-side tree comparison
- `tree_annotation_plot`: Trees with annotations

**Usage:**
```python
from metainformant.visualization import plot_phylo_tree, circular_tree_plot
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")
ax = plot_phylo_tree(tree)
ax = circular_tree_plot(tree)
```

#### Animations (`animations.py`)
Dynamic visualization animations.

**Functions:**
- `animate_time_series`: Time series animations
- `animate_evolution`: Evolutionary process animations
- `animate_clustering`: Clustering iteration animations
- `animate_network`: Network evolution animations
- `animate_trajectory`: Trajectory inference animations

**Usage:**
```python
from metainformant.visualization import animate_time_series

data = [1, 2, 3, 2, 4, 5, 3, 6]
fig, anim = animate_time_series(data, interval_ms=500)
```

### Utility Modules

#### Styling (`style.py`)
Publication-quality style presets and color palettes.

**Functions:**
- `apply_publication_style`: Apply publication-quality settings
- `get_color_palette`: Get color palettes (primary, accessible, colorblind, viridis_like)
- `get_figure_size`: Get predefined figure sizes
- `set_font_family`: Set font family
- `set_font_size`: Set font size
- `reset_style`: Reset to defaults

#### Layout (`layout.py`)
Multi-panel figure creation and layout management.

**Functions:**
- `create_subplot_grid`: Create subplot grids
- `create_multi_panel`: Create multi-panel layouts
- `add_shared_axis_labels`: Add shared labels
- `hide_unused_subplots`: Hide unused subplots
- `adjust_spacing`: Adjust subplot spacing

#### Export (`export.py`)
High-resolution figure export utilities.

**Functions:**
- `save_figure`: Save figure with high-resolution settings
- `save_figure_multiformat`: Save in multiple formats
- `batch_export_figures`: Batch export multiple figures
- `get_supported_formats`: Get list of supported formats

#### Interactive (`interactive.py`)
Plotly integration for interactive plots (optional dependency).

**Functions:**
- `create_interactive_scatter`: Interactive scatter plots
- `create_interactive_heatmap`: Interactive heatmaps
- `convert_matplotlib_to_plotly`: Convert matplotlib to Plotly
- `is_plotly_available`: Check Plotly availability

### Domain Integration Modules

#### GWAS Integration (`gwas_integration.py`)
Unified access to GWAS visualization functions.

Re-exports functions from `gwas.visualization_*` modules including:
- Manhattan plots, circular Manhattan plots
- Q-Q plots, lambda GC plots
- Regional plots, gene annotation plots
- PCA plots, kinship heatmaps
- Variant property plots
- Effect size plots

#### Single-Cell Integration (`singlecell_integration.py`)
Unified access to single-cell visualization functions.

Re-exports functions from `singlecell.visualization` including:
- QC metrics plots
- Embedding plots (UMAP, t-SNE, PCA)
- Gene expression plots

#### Information Theory Integration (`information_integration.py`)
Unified access to information theory visualization functions.

Re-exports functions from `information.visualization` including:
- Entropy plots
- Mutual information matrices
- Information profiles
- Rényi spectra

#### Life Events Integration (`life_events_integration.py`)
Unified access to life events visualization functions.

Re-exports functions from `life_events.visualization` including:
- Event timeline plots
- Event embeddings visualization
- Attention heatmaps
- Prediction importance plots

## Backward Compatibility

The module maintains backward compatibility through:

1. **`plots.py` wrapper**: Re-exports functions from new modules for existing code
2. **Unchanged function names**: All existing function names remain the same
3. **Import paths**: Both old and new import paths work

**Old imports (still work):**
```python
from metainformant.visualization.plots import lineplot, scatter_plot
```

**New imports (recommended):**
```python
from metainformant.visualization import lineplot, scatter_plot
# or
from metainformant.visualization.basic import lineplot, scatter_plot
```

## Integration Features

### Backend Support
- **Matplotlib**: Primary backend for static plots
- **Seaborn**: Enhanced statistical visualizations (optional)
- **Plotly**: Interactive web-based visualizations (optional)
- **NetworkX**: Network graph visualization (optional)

### Consistent API
All plotting functions follow a consistent pattern:
```python
ax = plot_function(data, **kwargs)
ax.set_xlabel("X Label")
ax.set_ylabel("Y Label")
ax.set_title("Plot Title")
```

## Performance Considerations

- **Memory Management**: Figures are properly closed after saving
- **Vector Graphics**: SVG/PDF export for scalable graphics
- **Batch Processing**: Efficient handling of multiple plots
- **Large Dataset Optimization**: Downsampling for large datasets

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import phylogeny
from metainformant.visualization import plot_phylo_tree

tree = phylogeny.neighbor_joining_tree(dna_sequences)
ax = plot_phylo_tree(tree)
```

### With GWAS Module
```python
from metainformant.gwas import association_test_linear
from metainformant.visualization import manhattan_plot, qq_plot

results = association_test_linear(genotype_matrix, phenotypes)
ax = manhattan_plot(results, 'position', 'neg_log10_p', 'chromosome')
```

### With Single-Cell Module
```python
from metainformant.singlecell import compute_umap
from metainformant.visualization import umap_plot

data = compute_umap(data, min_dist=0.1)
ax = umap_plot(data, x_col='UMAP1', y_col='UMAP2', hue='cluster')
```

## Testing

Comprehensive tests ensure:
- Plot rendering correctness
- Data accuracy in visualizations
- Animation functionality
- Export format validation
- Performance benchmarks

## Dependencies

- **Required**: matplotlib, numpy
- **Optional**: 
  - seaborn (enhanced styling)
  - plotly (interactive plots)
  - networkx (network visualization)
  - scipy (statistical functions)
  - scikit-learn (machine learning plots)

This module provides a visualization toolkit for biological data analysis and figure generation.
