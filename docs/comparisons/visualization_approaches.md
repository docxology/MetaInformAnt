# Visualization Approaches: Comprehensive Module Comparison

## Overview

The `visualization` module in METAINFORMANT is a comprehensive plotting system with **70+ plot types** organized into specialized submodules. This guide helps you choose the right visualization approach for your biological data, audience, and interactivity needs.

## Visualization Module Organization

```
visualization/
├── plots/                  # Core plotting functions (basic, general, specialized, multidim, animations)
├── genomics/               # Genomics-specific visualizations (genomics, expression, networks, trees)
├── analysis/               # Statistical, dimensionality reduction, time-series, information, quality plots
├── dashboards/             # Composite multi-panel figures and Plotly interactive figures
├── interactive_dashboards/ # JSON-serializable interactive plot data and HTML export
└── config/                 # Color palettes, themes, styling conventions
```

The sections below group related plots under short names (`basic/`, `statistical/`, `dimred/`, `trees/`, ...);
each example shows the real import path.

## Quick Decision Matrix

| Data Domain | Primary Plot Types | Key Submodule | Example Functions |
|-------------|-------------------|---------------|-------------------|
| **Genomics & Variants** | Manhattan, QQ, Volcano, Ideogram | `genomics/` | `manhattan_plot()`, `volcano_plot()`, `chromosome_ideogram()` |
| **Gene Expression** | Heatmaps, Volcano, DE plots | `expression/` | `plot_expression_heatmap()`, `plot_differential_expression()` |
| **Phylogenetics** | Trees (various layouts) | `trees/` | `plot_phylo_tree()`, `circular_tree_plot()` |
| **Networks & Pathways** | Network graphs, hierarchical layouts | `networks/` | `plot_network_basic()`, `plot_network_hierarchical()` |
| **Statistical Analysis** | Histograms, Boxplots, Violin, ROC | `statistical/` | `histogram()`, `box_plot()`, `violin_plot()`, `roc_curve()` |
| **Dimensionality Reduction** | PCA, UMAP, t-SNE plots | `dimred/` | `pca_plot()`, `plot_umap()`, `plot_tsne()` |
| **Time Series** | Line plots, Forecasts, Autocorrelation | `timeseries/` | `plot_time_series()`, `plot_forecast()`, `plot_autocorrelation()` |
| **Multi-Dimensional** | Pair plots, Parallel coords, Radar | `multidim/` | `plot_pairwise_relationships()`, `plot_parallel_coordinates()`, `plot_radar_chart()` |
| **Quality Control** | Quality scores, Adapter content, Length dist | `quality/` | `plot_per_base_quality_boxplot()`, `plot_adapter_content()` |
| **Information Theory** | Entropy, Mutual info profiles | `information/` | `plot_entropy_profile()`, `plot_mutual_information_matrix()` |
| **Simulations** | Animation of evolutionary/genetic processes | `animations/` | `animate_evolution()`, `animate_trajectory()` |
| **General Purpose** | Line, Scatter, Bar, Pie, Heatmap | `basic/` | `lineplot()`, `scatter_plot()`, `bar_plot()`, `heatmap()` |

## 1. Data Type → Plot Type Decision Tree

```mermaid
flowchart TD
    A[Start: What data do you have?] --> B{Data Type}
    
    B -->|Genomic positions<br/>SNPs, genes, peaks| Genomics[Genomics Plots<br/>manhattan/volcano/ideogram]
    B -->|Expression values<br/>genes × samples| Expression[Expression Plots<br/>heatmap/volcano/MA]
    B -->|Tree/phylogeny<br/>Newick format| Tree[Phylogenetic Trees<br/>rectangular/circular/radial]
    B -->|Network edges<br/>PPI, pathways| Network[Network Plots<br/>graph/circular/hierarchical]
    B -->|High-dimensional<br/>multi-omic features| DimRed[Dim. Reduction<br/>PCA/UMAP/t-SNE]
    B -->|Time points<br/>longitudinal| TimeSeries[Time Series<br/>line/forecast/ACF]
    B -->|Multiple variables<br/>per sample| MultiDim[Multi-Dim Plots<br/>pairplot/parallel coords]
    B -->|Quality metrics<br/>Phred scores, lengths| QC[QC Plots<br/>quality scores/adapters]
    B -->|Prob. distributions<br/>entropy/mutual info| Info[Information Theory<br/>entropy/MI profiles]
    B -->|Simple numeric<br/>vectors, matrices| Basic[Basic Plots<br/>line/scatter/bar/pie]
    
    Genomics --> C{Goal?}
    Expression --> D{Goal?}
    Tree --> E{Tree Type?}
    Network --> F{Network Size?}
    DimRed --> G{Method?}
    TimeSeries --> H{Temporal Pattern?}
    MultiDim --> I{Relationships?}
    QC --> J{QC Aspect?}
    Info --> K{Measure?}
    Basic --> L{Chart Type?}
    
    C -->|Identify peaks| G1[manhattan_plot]
    C -->|Compare effect vs sig| G2[volcano_plot]
    C -->|Show chromosomal dist| G3[chromosome_ideogram]
    
    D -->|Differential expression| E1[volcano_plot]
    D -->|Sample similarity| E2[expression_heatmap]
    D -->|Treat effect size| E3[plot_differential_expression]
    
    E -->|Small tree (<100 taxa)| T1[plot_phylo_tree]
    E -->|Large tree with collapse| T2[circular_tree_plot]
    E -->|Cladogram only| T3[plot_phylo_tree]
    
    F -->|Small network (<50 nodes)| N1[plot_network_basic]
    F -->|Pathway/diagram| N2[plot_network_hierarchical]
    F -->|Large network (>1000)| N3[plot_network_circular]
    
    G -->|Linear structure| D1[pca_plot]
    G -->|Nonlinear manifold| D2[plot_umap / plot_tsne]
    G -->|Show variance| D3[pca_scree_plot]
    
    H -->|Trend over time| TS1[plot_time_series]
    H -->|Seasonal/decompose| TS2[plot_seasonal_decomposition]
    H -->|Autocorrelation| TS3[plot_autocorrelation]
    
    I -->|Pairwise relationships| M1[plot_pairwise_relationships]
    I -->|Multi-var trajectories| M2[plot_parallel_coordinates]
    I -->|Radar/skill profile| M3[plot_radar_chart]
    
    J -->|Per-base quality| QC1[plot_per_base_quality_boxplot]
    J -->|Adapter contamination| QC2[plot_adapter_content]
    J -->|Read length dist| QC3[plot_length_distribution]
    
    K -->|Entropy by k-mer size| I1[plot_entropy_profile]
    K -->|Feature dependence| I2[plot_mutual_information_matrix]
    
    L -->|Trend over x| B1[lineplot]
    L -->|Correlation| B2[scatter_plot]
    L -->|Comparison| B3[bar_plot]
    L -->|Proportion| B4[pie_chart]
    L -->|2D density| B5[heatmap]
```

## 2. Comprehensive Plot Comparison Table

| Plot Type | Input Data | Use Case | Interactivity | Module | Complexity |
|-----------|-------------|-----------|---------------|--------|------------|
| **Manhattan** | GWAS summary stats (chr, pos, p) | Genome-wide association peaks | Low (static) | `genomics` | Medium |
| **Volcano** | Effect size vs p-value | Differential expression/association | Low | `genomics` | Low |
| **QQ** | Observed vs expected p-values | Inflation diagnostic | Low | `statistical` | Low |
| **Ideogram** | Chromosome ideogram with bands | Genomic context overview | Low | `genomics` | Medium |
| **Heatmap** | Matrix (samples × features) | Expression clusters, correlation | Medium (zoom) | `basic` | Low |
| **Expression Heatmap** | Expression matrix (genes × samples) | Gene expression patterns | High (gene search) | `expression` | Medium |
| **PCA Plot** | PCA scores matrix | Sample ordination, batch detection | Medium | `dimred` | Low |
| **UMAP/t-SNE** | High-dim embedding | Non-linear sample structure | Medium | `dimred` | Medium |
| **Scree Plot** | Variance per component | Dimensionality choice | Low | `dimred` | Low |
| **Phylogenetic Tree** | Newick tree file | Evolutionary relationships | Low | `trees` | High |
| **Network Graph** | Edge list, adjacency | PPI, regulatory networks | High (layout) | `networks` | High |
| **Pathway Diagram** | Pathway + expression data | Pathway activity visualization | Low | `networks` | High |
| **Line Plot** | x, y vectors | Time series, trends | Low | `basic` | Low |
| **Time Series Forecast** | Time-indexed series | Future predictions | Medium | `timeseries` | High |
| **Autocorrelation** | Time series | Seasonality detection | Low | `timeseries` | Medium |
| **Histogram** | Numeric vector | Distribution visualization | Low | `statistical` | Low |
| **Boxplot** | Numeric by group | Group comparison | Low | `statistical` | Low |
| **Violin Plot** | Numeric by group | Distribution + density | Low | `statistical` | Medium |
| **ROC Curve** | Pred/true labels | Classifier performance | Low | `statistical` | Low |
| **Pair Plot** | DataFrame (numeric cols) | Multivariate relationships | Medium | `multidim` | Medium |
| **Parallel Coord.** | DataFrame (all cols) | High-dim pattern discovery | High (brushing) | `multidim` | High |
| **Radar Chart** | Multi-axis variables | Skill/profile comparison | Low | `multidim` | Low |
| **Bar Plot** | Categorical + value | Group comparisons | Low | `basic` | Low |
| **Pie Chart** | Categorical proportions | Composition (limited use) | Low | `basic` | Low |
| **Quality Scores** | FASTQ quality matrix | Per-base quality | Low | `quality` | Medium |
| **Adapter Content** | Adapter positions | Contamination detection | Low | `quality` | Low |
| **Entropy Plot** | Sequence entropy | Complexity, conservation | Low | `information` | Medium |
| **Mutual Information** | Variable pairs | Dependence visualization | Medium | `information` | High |
| **Animation** | Time/lambda series | Dynamic processes | Very High | `animations` | Very High |

### Plot Selection Matrix by Domain

| Domain → Plot | Genomics | Transcriptomics | Proteomics | Networks | Phylogeny | Single-Cell | Spatial |
|---------------|----------|-----------------|------------|----------|-----------|-------------|---------|
| Manhattan ✓ | ✓ (GWAS) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| Volcano ✓ | ✓ | ✓ (DE) | ✓ (DE) | ✗ | ✗ | ✓ | ✓ |
| Heatmap ✓ | ✓ (peaks) | ✓ (expr) | ✓ (prot) | ✓ (adjacency) | ✗ | ✓ (markers) | ✓ (spots) |
| PCA ✓ | ✓ (pop struct) | ✓ (batch) | ✓ (batch) | ✓ (embedding) | ✗ | ✓ | ✓ |
| Network ✓ | ✗ | ✗ | ✗ | ✓ (PPI) | ✗ | ✗ | ✗ |
| Tree ✓ | ✗ | ✗ | ✗ | ✗ | ✓ (phylo) | ✗ | ✗ |
| TimeSeries ✓ | ✗ | ✓ (kinetics) | ✗ | ✗ | ✗ | ✗ | ✗ |
| DimRed ✓ | ✓ (snp-pca) | ✓ (scRNA UMAP) | ✓ | ✓ | ✗ | ✓ (cell clusters) | ✓ (tissue domains) |
| QC ✓ | ✓ (variant) | ✓ (fastq/expr) | ✓ (mass-spec) | ✗ | ✗ | ✓ (mito%, ribo%) | ✓ (spot quality) |

## 3. Interactivity & Output Format Guide

### Static vs Interactive Decision

| Need | Recommended Approach | Implementation |
|------|---------------------|----------------|
| **Publication figure** | Static (PNG, TIFF, PDF, SVG) | `savefig("fig.pdf")` |
| **Exploratory analysis** | Interactive (Plotly) | `interactive_scatter()`, `interactive_heatmap()` |
| **Dashboard / report** | Composite layout | `multi_panel([panel1, panel2])` |
| **Animation for talks** | GIF animation | `animate_trajectory(..., output_path="trajectory.gif")` |
| **Web app / Shiny** | Interactive Plotly | Export to JSON, embed in Streamlit |

### Output Format Matrix

| Plot Type | Static Formats | Interactive Formats | Resolution | File Size |
|-----------|----------------|---------------------|------------|-----------|
| Manhattan | PDF, PNG (300 DPI) | Plotly HTML | 10–20 MB (PDF) | N/A |
| Heatmap (large) | TIFF (600 DPI) | None (too heavy) | 50–100 MB | N/A |
| Network (small) | SVG (vector) | Cytoscape JSON | <1 MB | Small |
| Tree | Newick (tree only), PDF | N/A | <5 MB | N/A |
| Animation | GIF | HTML5 Canvas | 50–200 MB | Medium |
| Dashboard | Multi-page PDF | Bokeh/Streamlit app | N/A | Large |

## 4. Module-Specific Visualization Submodules

### genomics/ Submodule

Genomics-specific plots for variant and expression analysis.

```python-snippet
from metainformant.visualization.genomics.genomics import (
    manhattan_plot,
    volcano_plot,
    regional_plot,
    chromosome_ideogram,
    circular_manhattan_plot
)

# Manhattan plot: GWAS results
manhattan_plot(
    gwas_df,
    chr_col='chrom',
    pos_col='pos',
    pval_col='p_value'
)

# Regional plot: Zoom into locus (expects CHR, BP and P columns)
regional_plot(
    gwas_df,
    chr='1',
    start=10500000,
    end=10600000
)
```

**When to use**:
- GWAS results → Manhattan, QQ
- Differential expression → Volcano
- Genomic context → Ideogram, circular Manhattan

### expression/ Submodule

Gene expression visualization across samples/conditions.

```python-snippet
from metainformant.visualization.genomics.expression import (
    plot_expression_heatmap,
    plot_differential_expression,
    plot_enrichment_barplot
)

# Heatmap of top DE genes (z-score normalized)
plot_expression_heatmap(expr_df.loc[top_de_genes])

# Differential expression plot: log-FC vs -log10 adjusted p-value,
# point size scaled by mean expression (log2FoldChange, padj and baseMean columns)
plot_differential_expression(de_df)

# Pathway enrichment results (Term, p.adjust and Count columns)
plot_enrichment_barplot(enrichment_df)
```

**When to use**:
- Bulk RNA-seq DE results → Volcano, DE plot, Heatmap
- Expression patterns → Clustered heatmaps
- Sample relationships → PCA, correlation matrix

### trees/ Submodule

Phylogenetic tree visualization in multiple layouts.

```python-snippet
from metainformant.visualization.genomics.trees import (
    plot_phylo_tree,
    circular_tree_plot,
    unrooted_tree_plot
)
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")

# Standard hierarchical tree
plot_phylo_tree(tree, with_labels=True)

# Circular layout for large trees
circular_tree_plot(tree)

# Unrooted (spring-layout) view
unrooted_tree_plot(tree)
```

**When to use**:
- Evolutionary relationships → Rectangular/circular trees
- Large trees (>50 taxa) → Circular layout with collapse
- Cladograms only → Rectangular without branch lengths

### networks/ Submodule

Biological network visualization (PPI, regulatory, pathways).

```python-snippet
from metainformant.visualization.genomics.networks import (
    plot_network_basic,
    plot_network_circular,
    plot_network_hierarchical,
    plot_community_network
)

# Network graph with spring layout (NetworkX graph)
plot_network_basic(
    networkx_graph,
    node_color='steelblue'
)

# Hierarchical layout for pathway-like graphs
plot_network_hierarchical(pathway_graph)

# Color nodes by community membership ({node: community_id})
plot_community_network(
    networkx_graph,
    communities
)
```

**When to use**:
- Small-medium PPI/network → Spring layout
- Pathway diagrams → Hierarchical layout
- Large networks (>1000 nodes) → Circular layout

### statistical/ Submodule

Statistical diagnostic and summary plots.

```python-snippet
from metainformant.visualization.analysis.statistical import (
    histogram,
    box_plot,
    violin_plot,
    qq_plot,
    roc_curve,
    density_plot
)

# Q-Q plot for GWAS p-values (uniform) or model residuals ('norm')
qq_plot(
    gwas_df['p_value'].to_numpy(),
    distribution='uniform'
)

# ROC curve for classifier (the AUC is shown in the legend)
roc_curve(
    true_labels,
    predictions
)

# Violin plot for expression by group (one array per group)
violin_plot([expr_control, expr_treated])
```

**When to use**:
- Distribution assessment → Histogram, density
- Group comparison → Boxplot, violin
- Model diagnostics → QQ, residual plots
- Classifier evaluation → ROC, precision-recall

### dimred/ Submodule

Dimensionality reduction visualizations.

```python-snippet
from metainformant.visualization.plots.general import pca_plot
from metainformant.visualization.analysis.dimred import plot_umap, plot_tsne
from metainformant.gwas.visualization.population import pca_scree_plot

# PCA scores plot (DataFrame with PC1, PC2, PC1_variance, PC2_variance and a 'batch' column)
pca_plot(
    pca_df,
    hue='batch'
)

# UMAP for non-linear embedding of a feature matrix (requires umap-learn)
plot_umap(
    feature_matrix,
    c=cluster_codes,
    s=5
)

# Scree plot: variance per PC (cumulative curve included)
pca_scree_plot(
    pca.explained_variance_ratio_,
    output_file='scree.png'
)
```

**When to use**:
- Linear structure → PCA (faster, interpretable)
- Non-linear manifolds → UMAP (preserves local) or t-SNE (probabilistic)
- Variance selection → Scree plot
- Batch detection → Color samples by batch in PCA

### quality/ Submodule

Quality control visualizations for sequencing data.

```python-snippet
from metainformant.visualization.analysis.quality import (
    plot_per_base_quality_boxplot,
    plot_adapter_content,
    plot_length_distribution,
    plot_gc_distribution
)

# Per-base quality scores (FastQC-style): position range -> list of quality scores
plot_per_base_quality_boxplot(per_base_qualities)

# Adapter contamination: adapter name -> content (%) at each position
plot_adapter_content(adapter_content_by_name)

# GC content distribution (list of per-read GC percentages)
plot_gc_distribution(gc_percent)
```

**When to use**:
- FASTQ QC → Quality scores, adapter content, length dist
- BAM QC → Coverage depth, insert size
- Assembly QC → N50 plot, BUSCO scores

### information/ Submodule

Information theory visualizations.

```python-snippet
from metainformant.visualization.analysis.information import (
    plot_entropy_profile,
    plot_mutual_information_matrix
)

# Shannon entropy of each sequence at several k-mer sizes
plot_entropy_profile(
    entropy_by_sequence  # {sequence name: [H(k=1), H(k=2), ...]}
)

# Mutual information between positions (compares correlated sites)
plot_mutual_information_matrix(
    mi_scores,
    labels=position_labels
)
```

**When to use**:
- Sequence complexity → Entropy profile
- Co-evolving residues → MI matrix

### animations/ Submodule

Animated visualizations for dynamic processes.

```python-snippet
from metainformant.visualization.plots.animations import (
    animate_time_series,
    animate_evolution,
    animate_clustering,
    animate_network
)

# Animate time-series (array with one row per series; saved as GIF)
fig, anim = animate_time_series(
    y,
    interval=50,
    output_path="time_series.gif"
)

# Sequence evolution animation (one sequence per generation)
animate_evolution(
    sequence_history,
    output_path="evolution.gif"
)

# Clustering progression: 2D points plus the cluster labels at each step
animate_clustering(
    points_2d,
    cluster_labels_over_time,
    output_path="trajectory.gif"
)
```

**When to use**:
- Temporal processes → `animate_time_series`
- Evolutionary algorithms → `animate_evolution`
- Trajectory inference → `animate_clustering`
- Network dynamics → `animate_network`

## 5. Audience & Complexity Guide

### ForPublications (High-Quality Static)

**Recommended**: `plots/` + `config/` with publication theme

```python-snippet
from metainformant.visualization.config.themes import apply_theme
apply_theme("publication")  # Sans-serif 10pt; DPI: 300

# Create figure
from metainformant.visualization.plots import manhattan_plot
ax = manhattan_plot(gwas_df, output_path="manhattan.pdf")
```

**Key modules**: `basic/`, `genomics/`, `statistical/`, `trees/`

---

### For Exploratory Analysis (Fast, Simple)

**Recommended**: Default `basic/` functions, low-code

```python
from metainformant.visualization.plots.basic import scatter_plot
from metainformant.visualization.analysis.statistical import histogram
scatter_plot(df['PC1'].to_numpy(), df['PC2'].to_numpy(), c=df['batch'].astype('category').cat.codes)
histogram(df['expression'].to_numpy(), bins=50)
```

**Key modules**: `basic/`, `dimred/`, `statistical/`

---

### For Interactive Exploration (Web-based)

**Recommended**: Plotly figures from `dashboards.interactive`

```python-snippet
from metainformant.visualization.dashboards.interactive import interactive_scatter

fig = interactive_scatter(embedding_df, x="PC1", y="PC2", color="batch")
fig.show()  # Opens in browser (requires plotly)
```

**Key modules**: `dashboards/` (`interactive_scatter`, `interactive_heatmap`, `interactive_volcano`, `interactive_manhattan`)

---

### For Dashboards (Composite Figures)

**Recommended**: `dashboards/` module

```python-snippet
from metainformant.visualization.dashboards.composite import multi_panel
from metainformant.visualization.plots.general import manhattan_plot, qq_plot, volcano_plot

fig = multi_panel(
    [
        {"func": manhattan_plot, "kwargs": {"data": gwas_df}, "title": "Manhattan"},
        {"func": qq_plot, "kwargs": {"data": gwas_df['P'].to_numpy(), "distribution": "uniform"}, "title": "Q-Q"},
        {"func": volcano_plot, "kwargs": {"data": de_df}, "title": "Volcano"},
    ],
    ncols=2,
    output_path="gwas_dashboard.pdf",
)
```

**Key modules**: `dashboards/`, `interactive_dashboards/`

---

### For Presentations (Animated / Large Print)

**Recommended**: `animations/` + the presentation theme

```python-snippet
from metainformant.visualization.config.themes import apply_theme
from metainformant.visualization.plots.animations import animate_time_series

apply_theme("presentation")  # Larger fonts, thicker lines

animate_time_series(data, output_path="trend.gif")
```

**Key modules**: `animations/`, `basic/`, `genomics/`

## 6. Performance & Scalability

| Plot Type | Time Complexity | Memory (n points) | Max Reasonable Points | Notes |
|-----------|----------------|-------------------|-----------------------|-------|
| Manhattan | O(n_snps) | O(n_snps) | 10⁷ SNPs (genome-wide) | Fast, uses binning |
| Heatmap | O(n×m) | O(n×m) | 10⁴ × 10⁴ (100M cells) | Large: use sparse or downsampling |
| Network (spring) | O(V³) worst-case | O(V²) | 1,000 nodes | Large: use `plot_network_circular` |
| PCA/UMAP | O(n × p × n_components) | O(n × p) | 100K samples × 10K features | Use incremental PCA for larger |
| Tree | O(n_taxa) | O(n_taxa²) | 1,000 taxa (uncluttered) | Large: use circular with collapse |
| Animation | O(n_frames × n_points) | O(n_points) | Depends on frames real-time | Pre-render for large |
| Pairplot | O(n² × p²) | O(n × p) | 5K samples × 10 vars | Large: use `sns.pairplot(sample=1000)` |

**Optimization tips**:
- Downsample before plotting if n > 10,000 points (scatter)
- Use hexbin or 2D histogram for dense scatter plots
- For large heatmaps, downsample or clip values
- Network: use `plot_network_circular` for large graphs
- Animation: pre-render frames, don't real-time compute

## 7. Styling & Customization

### Theme Selection

| Theme | Use Case | Colors | Font |
|-------|----------|--------|------|
| `publication` | Journal figures | Matplotlib default cycle, light grid, no top/right spines | Sans-serif 10pt, 300 DPI |
| `presentation` | Slides | Matplotlib default cycle, thicker lines | Sans-serif 14pt |
| `poster` | Posters | Matplotlib default cycle, thickest lines | Sans-serif 18pt |
| `dark` | Dark background talks | Dark grey background, light text | Matplotlib default, light grey |
| `minimal` | Clean figures without grid | Matplotlib default cycle | Sans-serif 10pt |
| `colorblind` | Universal accessibility | Okabe-Ito color cycle | Matplotlib default, 10pt |

`reset_theme()` restores the matplotlib defaults.

```python-snippet
from metainformant.visualization.config.themes import apply_theme
from metainformant.visualization.config.palettes import categorical

apply_theme('publication')
colors = categorical(8, palette='wong')  # Okabe-Ito scale
```

### Color Palette Guide

| Palette Type | Best For | Example |
|--------------|----------|---------|
| `categorical` | Group differences (≤12 groups) | Sex, population, cluster |
| `continuous` | Quantitative values (heatmap) | p-values, expression Z-score |
| `diverging` | Values centered at 0 | log-FC, Z-score |
| `sequential` | Monotonic increase | Quality scores, depth |
| `colorblind` | Universal accessibility | Red-green safe (Okabe-Ito) |

## 8. Example: Multi-Panel Figure Assembly

```python-snippet
from metainformant.visualization.dashboards.composite import multi_panel
from metainformant.visualization.genomics.genomics import manhattan_plot, volcano_plot
from metainformant.visualization.analysis.statistical import histogram, qq_plot

# Create 2×2 figure grid (panels fill the grid row by row)
fig = multi_panel(
    [
        {"func": manhattan_plot, "kwargs": {"data": gwas_df}},                                        # top-left
        {"func": qq_plot, "kwargs": {"data": gwas_df['P'].to_numpy(), "distribution": "uniform"}},  # top-right
        {"func": histogram, "kwargs": {"data": gwas_df['beta'].to_numpy()}},                        # bottom-left
        {"func": volcano_plot, "kwargs": {"data": de_df}},                                            # bottom-right
    ],
    ncols=2,
    figsize=(12, 10),
    output_path="gwas_manuscript_figure.pdf",
)
```

## 9. Cross-Module Visualization Usage Examples

### Equivalence: GWAS Plotting Across Modules

```python-snippet
# Option 1: Direct from gwas module (has built-in viz)
# gwas_results: list of dicts with chrom, pos and p_value keys
from metainformant.gwas.visualization import manhattan_plot, qq_plot
manhattan_plot(gwas_results)
qq_plot([r['p_value'] for r in gwas_results])

# Option 2: Using unified visualization module (DataFrame input)
from metainformant.visualization.genomics.genomics import manhattan_plot
manhattan_plot(
    gwas_df,
    chr_col='chrom',
    pos_col='pos',
    pval_col='p_value'
)
```

### Equivalence: Expression Visualization

```python-snippet
# z-score normalized heatmap (genomics submodule)
from metainformant.visualization.genomics.expression import plot_expression_heatmap
plot_expression_heatmap(expression_matrix)

# General-purpose heatmap of the values as given (plots submodule)
from metainformant.visualization.plots.general import expression_heatmap
expression_heatmap(expression_matrix, cmap="RdBu_r")
```

## 10. Advanced: Custom Plot Development

If existing plots don't fit your needs, write a plain function that follows the conventions of `visualization.plots`
(optional `ax` and `output_path` arguments, returns the `Axes`), so it also works as a `multi_panel` panel:

```python-snippet
import matplotlib.pyplot as plt
import numpy as np
from metainformant.visualization.config.conventions import save_figure_deterministic

def custom_manhattan_plot(data, highlight_genes=None, *, ax=None, output_path=None):
    highlight_genes = highlight_genes or []
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))
    # Custom logic
    for chrom in data['chrom'].unique():
        chrom_data = data[data['chrom'] == chrom]
        ax.scatter(chrom_data['pos'], -np.log10(chrom_data['p']))
    # Highlight specific genes
    for gene in highlight_genes:
        gene_pos = data[data['gene'] == gene]['pos'].iloc[0]
        ax.axvline(gene_pos, color='red', linestyle='--')
    if output_path:
        save_figure_deterministic(ax.figure, output_path)
    return ax

# Use it
ax = custom_manhattan_plot(gwas_df, highlight_genes=['APOE', 'CLU'])
```

## 11. Summary Table: All Visualization Submodules

| Submodule | Plot Count | Primary Domains | Key Functions | Best For |
|-----------|------------|-----------------|---------------|----------|
| `basic/` | 12 | General | line, scatter, bar, pie, area, step, heatmap | Quick exploratory plots |
| `statistical/` | 10 | Stats | histogram, boxplot, violin, qq, roc, density | Statistical diagnostics |
| `genomics/` | 8 | GWAS, genomics | manhattan, volcano, regional, ideogram, circular Manhattan | Variant/expression studies |
| `expression/` | 6 | RNA-seq, proteomics | plot_expression_heatmap, plot_differential_expression, plot_enrichment_barplot | Gene/protein expression |
| `dimred/` | 6 | ML, scRNA | pca, umap, tsne, scree | Dimensionality reduction viz |
| `networks/` | 5 | PPI, pathways | plot_network_basic, plot_community_network | Biological networks |
| `trees/` | 4 | Phylogeny | plot_phylo_tree, circular_tree_plot | Evolutionary trees |
| `timeseries/` | 5 | Temporal | line, forecast, autocorrelation, seasonal | Time-series data |
| `multidim/` | 4 | High-dim | pairplot, parallel_coords, radar | Multi-variable relationships |
| `quality/` | 4 | QC | plot_per_base_quality_boxplot, plot_adapter_content, plot_length_distribution | Sequencing QC |
| `information/` | 3 | Info theory | plot_entropy_profile, plot_mutual_information_matrix | Info theoretic analysis |
| `animations/` | 6 | Dynamic | animate_time_series, animate_evolution | Animated figures |
| `dashboards/` | N/A | Composite | multi_panel, genomic_overview, qc_summary | Multi-plot figures |
| `config/` | N/A | Theming | apply_theme, categorical | Figure aesthetics |

**Total**: 70+ individual plot functions across all submodules.

## 12. Related Documentation

- **[visualization/index.md](../visualization/index.md)** — Visualization module overview
- **[visualization/README.md](../visualization/README.md)** — Full function reference
- **[visualization/genomics.md](../visualization/genomics.md)** — Genomics-specific plots
- **[visualization/statistical.md](../visualization/statistical.md)** — Statistical plots
- **[visualization/trees.md](../visualization/trees.md)** — Phylogenetic tree visualization
- **[visualization/animations.md](../visualization/animations.md)** — Animation guide
- **[COMPARISON_GUIDES.md](../COMPARISON_GUIDES.md)** — Master comparison index
