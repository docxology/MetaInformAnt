# Visualization Approaches: Comprehensive Module Comparison

## Overview

The `visualization` module in METAINFORMANT is a comprehensive plotting system with **70+ plot types** organized into specialized submodules. This guide helps you choose the right visualization approach for your biological data, audience, and interactivity needs.

## Visualization Module Organization

```
visualization/
├── plots/           # Core plotting functions (basic, general, specialized)
├── genomics/        # Genomics-specific visualizations
├── analysis/        # Statistical, dimensionality reduction, quality plots
├── dashboards/      # Composite multi-panel figures
├── config/          # Color palettes, themes, styling
├── animations/      # Animated time-series, evolutionary, network plots
└── integration/     # Domain-specific integration (GWAS, single-cell, etc.)
```

## Quick Decision Matrix

| Data Domain | Primary Plot Types | Key Submodule | Example Functions |
|-------------|-------------------|---------------|-------------------|
| **Genomics & Variants** | Manhattan, QQ, Volcano, Ideogram | `genomics/` | `manhattan_plot()`, `volcano_plot()`, `ideogram()` |
| **Gene Expression** | Heatmaps, Volcano, MA plots | `expression/` | `expression_heatmap()`, `ma_plot()` |
| **Phylogenetics** | Trees (various layouts) | `trees/` | `plot_phylo_tree()`, `plot_tree_circular()` |
| **Networks & Pathways** | Network graphs, pathway diagrams | `networks/` | `plot_network()`, `plot_pathway()` |
| **Statistical Analysis** | Histograms, Boxplots, Violin, ROC | `statistical/` | `histogram()`, `boxplot()`, `violin_plot()`, `roc_curve()` |
| **Dimensionality Reduction** | PCA, UMAP, t-SNE plots | `dimred/` | `pca_plot()`, `umap_plot()`, `tsne_plot()` |
| **Time Series** | Line plots, Forecasts, Autocorrelation | `timeseries/` | `lineplot()`, `timeseries_forecast()`, `autocorrelation()` |
| **Multi-Dimensional** | Pair plots, Parallel coords, Radar | `multidim/` | `pairplot()`, `parallel_coordinates()`, `radar_chart()` |
| **Quality Control** | Quality scores, Adapter content, Length dist | `quality/` | `quality_scores()`, `adapter_content()` |
| **Information Theory** | Entropy, Mutual info profiles | `information/` | `entropy_plot()`, `mutual_information_heatmap()` |
| **Simulations** | Animation of evolutionary/genetic processes | `animations/` | `animate_evolution()`, `animate_trajectory()` |
| **General Purpose** | Line, Scatter, Bar, Pie, Heatmap | `basic/` | `lineplot()`, `scatter_plot()`, `barplot()`, `heatmap()` |

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
    C -->|Show chromosomal dist| G3[ideogram]
    
    D -->|Differential expression| E1[volcano_plot]
    D -->|Sample similarity| E2[expression_heatmap]
    D -->|Treat effect size| E3[ma_plot]
    
    E -->|Small tree (<100 taxa)| T1[plot_phylo_tree]
    E -->|Large tree with collapse| T2[plot_tree_circular]
    E -->|Cladogram only| T3[plot_tree_rectangular]
    
    F -->|Small network (<50 nodes)| N1[plot_network spring]
    F -->|Pathway/diagram| N2[plot_network hierarchical]
    F -->|Large network (>1000)| N3[plot_network_circular]
    
    G -->|Linear structure| D1[pca_plot]
    G -->|Nonlinear manifold| D2[umap_plot / tsne_plot]
    G -->|Show variance| D3[scree_plot]
    
    H -->|Trend over time| TS1[lineplot with CI]
    H -->|Seasonal/decompose| TS2[seasonal_decompose]
    H -->|Autocorrelation| TS3[autocorrelation_plot]
    
    I -->|Pairwise relationships| M1[pairplot]
    I -->|Multi-var trajectories| M2[parallel_coordinates]
    I -->|Radar/skill profile| M3[radar_chart]
    
    J -->|Per-base quality| QC1[quality_scores]
    J -->|Adapter contamination| QC2[adapter_content]
    J -->|Read length dist| QC3[length_distribution]
    
    K -->|Entropy per position| I1[entropy_plot]
    K -->|Feature dependence| I2[mutual_information_heatmap]
    K -->|Info profile| I3[information_content_profile]
    
    L -->|Trend over x| B1[lineplot]
    L -->|Correlation| B2[scatter_plot]
    L -->|Comparison| B3[barplot / grouped_bar]
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
| **Exploratory analysis** | Interactive (Plotly) | `interactive=True` parameter |
| **Dashboard / report** | Composite layout | `create_dashboard([plot1, plot2])` |
| **Animation for talks** | MP4/GIF animation | `animate_trajectory(..., save_as="mp4")` |
| **Web app / Shiny** | Interactive Plotly | Export to JSON, embed in Streamlit |

### Output Format Matrix

| Plot Type | Static Formats | Interactive Formats | Resolution | File Size |
|-----------|----------------|---------------------|------------|-----------|
| Manhattan | PDF, PNG (300 DPI) | Plotly HTML | 10–20 MB (PDF) | N/A |
| Heatmap (large) | TIFF (600 DPI) | None (too heavy) | 50–100 MB | N/A |
| Network (small) | SVG (vector) | Cytoscape JSON | <1 MB | Small |
| Tree | Newick (tree only), PDF | N/A | <5 MB | N/A |
| Animation | MP4 (1080p) | HTML5 Canvas | 50–200 MB | Medium |
| Dashboard | Multi-page PDF | Bokeh/Streamlit app | N/A | Large |

## 4. Module-Specific Visualization Submodules

### genomics/ Submodule

Genomics-specific plots for variant and expression analysis.

```python-snippet
from metainformant.visualization.genomics import (
    manhattan_plot,
    volcano_plot,
    regional_plot,
    ideogram,
    circos_plot
)

# Manhattan plot: GWAS results
manhattan_plot(
    data=gwas_df,
    position_col='pos',
    pval_col='p_value',
    chromosome_col='chrom',
    highlight_peaks=True,
    genome_build='GRCh38'
)

# Regional plot: Zoom into locus
regional_plot(
    summary_stats=gwas_df,
    locus='chr1:10500000-10600000',
    annotation='genes.gtf'
)
```

**When to use**:
- GWAS results → Manhattan, QQ
- Differential expression → Volcano, MA
- Genomic context → Ideogram, Circos

### expression/ Submodule

Gene expression visualization across samples/conditions.

```python-snippet
from metainformant.visualization.expression import (
    expression_heatmap,
    ma_plot,
    enrichment_plot
)

# Heatmap of top DE genes
expression_heatmap(
    expression_matrix=expr_df,
    gene_list=top_de_genes,
    sample_annotation=metadata_df
)

# MA plot: log-FC vs mean expression
ma_plot(
    de_results=de_df,
    mean_col='base_mean',
    logfc_col='log2FoldChange',
    sig_col='padj'
)
```

**When to use**:
- Bulk RNA-seq DE results → Volcano, MA, Heatmap
- Expression patterns → Clustered heatmaps
- Sample relationships → PCA, correlation matrix

### trees/ Submodule

Phylogenetic tree visualization in multiple layouts.

```python-snippet
from metainformant.visualization.trees import (
    plot_phylo_tree,
    plot_tree_rectangular,
    plot_tree_circular,
    plot_tree_radial
)
from Bio import Phylo

tree = Phylo.read("tree.nwk", "newick")

# Standard rectangular tree
plot_phylo_tree(tree, branch_labels=True)

# Circular layout for large trees
plot_tree_circular(tree, label_angle=45)

# Radial cladogram
plot_tree_radial(tree, show_confidence=True)
```

**When to use**:
- Evolutionary relationships → Rectangular/circular trees
- Large trees (>50 taxa) → Circular layout with collapse
- Cladograms only → Rectangular without branch lengths

### networks/ Submodule

Biological network visualization (PPI, regulatory, pathways).

```python-snippet
from metainformant.visualization.networks import (
    plot_network,
    plot_network_circular,
    plot_network_hierarchical,
    plot_pathway
)

# Network graph with spring layout
plot_network(
    G=networkx_graph,
    node_colors=node_attributes,
    edge_weights=edge_list,
    layout='spring'
)

# Pathway diagram with expression overlay
plot_pathway(
    pathway_graph=graphml_file,
    node_values=expression_zscores,
    highlight_nodes=significant_genes
)
```

**When to use**:
- Small-medium PPI/network → Spring layout
- Pathway diagrams → Hierarchical layout
- Large networks (>1000 nodes) → Circular/fa layout

### statistical/ Submodule

Statistical diagnostic and summary plots.

```python-snippet
from metainformant.visualization.statistical import (
    histogram,
    boxplot,
    violin_plot,
    qq_plot,
    roc_curve,
    density_plot
)

# Q-Q plot for GWAS or model residuals
qq_plot(
    observed_pvalues=gwas_df['p_value'],
    theoretical=stats.uniform(0,1)
)

# ROC curve for classifier
roc_curve(
    y_true=true_labels,
    y_score=predictions,
    auc=True
)

# Violin plot for expression by group
violin_plot(
    data=expr_long_df,
    x='group',
    y='expression',
    hue='treatment'
)
```

**When to use**:
- Distribution assessment → Histogram, density
- Group comparison → Boxplot, violin
- Model diagnostics → QQ, residual plots
- Classifier evaluation → ROC, precision-recall

### dimred/ Submodule

Dimensionality reduction visualizations.

```python-snippet
from metainformant.visualization.dimred import (
    pca_plot,
    umap_plot,
    tsne_plot,
    scree_plot
)

# PCA scores plot
pca_plot(
    pca_result=pca_scores,
    variance_explained=pca.explained_variance_ratio_,
    color_by=metadata['batch']
)

# UMAP for non-linear embedding
umap_plot(
    embedding=umap_coords,
    labels=cell_types,
    point_size=5
)

# Scree plot: variance per PC
scree_plot(
    variance_ratio=pca.explained_variance_ratio_,
    cumulative=True
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
from metainformant.visualization.quality import (
    quality_scores,
    adapter_content,
    length_distribution,
    gc_content
)

# Per-base quality scores (FastQC-style)
quality_scores(
    quality_matrix=fastq_quality_array,
    per_base=True
)

# Adapter contamination across read length
adapter_content(
    adapter_counts=adapter_df,
    read_lengths=lengths
)

# GC content distribution
gc_content(
    gc_percent=gc_array,
    theoretical=expected_gc
)
```

**When to use**:
- FASTQ QC → Quality scores, adapter content, length dist
- BAM QC → Coverage depth, insert size
- Assembly QC → N50 plot, BUSCO scores

### information/ Submodule

Information theory visualizations.

```python-snippet
from metainformant.visualization.information import (
    entropy_plot,
    mutual_information_heatmap,
    information_content_profile
)

# Shannon entropy along sequence alignment
entropy_plot(
    alignment_matrix=msa_array,  # sequences × positions
    alphabet='DNA'
)

# Mutual information between positions (compares correlated sites)
mutual_information_heatmap(
    mi_matrix=mi_scores,
    annotation=structure_contacts
)

# Information content profile (sequence logo alternative)
information_content_profile(
    pfm=position_frequency_matrix,
    background=bg_freq
)
```

**When to use**:
- Sequence conservation → Entropy plot
- Co-evolving residues → MI heatmap
- Motif strength → Information content profile

### animations/ Submodule

Animated visualizations for dynamic processes.

```python-snippet
from metainformant.visualization.animations import (
    animate_time_series,
    animate_evolution,
    animate_clustering,
    animate_network
)

# Animate time-series with confidence bands
fig, anim = animate_time_series(
    timepoints=t,
    values=y,
    confidence=ci,
    interval=50,
    save_as="time_series.mp4"
)

# Evolutionary simulation animation
animate_evolution(
    population=pop_history,
    fitness_landscape=fitness_func,
    generations=100,
    save_as="evolution.mp4"
)

# Clustering progression (e.g., cell trajectory)
animate_clustering(
    embeddings=umap_per_timepoint,
    labels=cell_labels,
    save_as="trajectory.mp4"
)
```

**When to use**:
- Temporal processes → `animate_time_series`
- Evolutionary algorithms → `animate_evolution`
- Trajectory inference → `animate_clustering`
- Network dynamics → `animate_network`

## 5. Audience & Complexity Guide

### ForPublications (High-Quality Static)

**Recommended**: `plots/` + `styling/` with publication theme

```python-snippet
from metainformant.visualization import set_publication_style
set_publication_style()  # Font: Helvetica, 8pt; DPI: 300

# Create figure
from metainformant.visualization.plots import manhattan_plot
ax = manhattan_plot(gwas_df, save="manhattan.pdf")
```

**Key modules**: `basic/`, `genomics/`, `statistical/`, `trees/`

---

### For Exploratory Analysis (Fast, Simple)

**Recommended**: Default `basic/` functions, low-code

```python
from metainformant.visualization import scatter_plot, histogram
scatter_plot(x=df['PC1'], y=df['PC2'], hue=df['batch'])
histogram(df['expression'], bins=50)
```

**Key modules**: `basic/`, `dimred/`, `statistical/`

---

### For Interactive Exploration (Web-based)

**Recommended**: Enable Plotly backend

```python-snippet
from metainformant.visualization import enable_interactive
enable_interactive()  # Switches to Plotly

fig = pca_plot(embedding, labels=metadata, interactive=True)
fig.show()  # Opens in browser
```

**Key modules**: `basic/`, `dimred/`, `multidim/` (interactive=True)

---

### For Dashboards (Composite Figures)

**Recommended**: `dashboards/` module

```python-snippet
from metainformant.visualization.dashboards import (
    create_dashboard,
    multi_panel_figure
)

fig = create_dashboard([
    (manhattan_plot, {"data": gwas_df}),
    (qq_plot, {"pvalues": gwas_df['p']}),
    (volcano_plot, {"de_results": de_df})
], layout="2x2")
fig.save("gwas_dashboard.pdf")
```

**Key modules**: `dashboards/`, `integration/`

---

### For Presentations (Animated / Large Print)

**Recommended**: `animations/` + high-contrast styling

```python-snippet
from metainformant.visualization import set_presentation_style
set_presentation_style()  # Larger fonts, bold colors

animate_time_series(data, save_as="trend.mp4")
```

**Key modules**: `animations/`, `basic/`, `genomics/`

## 6. Performance & Scalability

| Plot Type | Time Complexity | Memory (n points) | Max Reasonable Points | Notes |
|-----------|----------------|-------------------|-----------------------|-------|
| Manhattan | O(n_snps) | O(n_snps) | 10⁷ SNPs (genome-wide) | Fast, uses binning |
| Heatmap | O(n×m) | O(n×m) | 10⁴ × 10⁴ (100M cells) | Large: use sparse or downsampling |
| Network (spring) | O(V³) worst-case | O(V²) | 1,000 nodes | Large: use `fa` or `circular` |
| PCA/UMAP | O(n × p × n_components) | O(n × p) | 100K samples × 10K features | Use incremental PCA for larger |
| Tree | O(n_taxa) | O(n_taxa²) | 1,000 taxa (uncluttered) | Large: use circular with collapse |
| Animation | O(n_frames × n_points) | O(n_points) | Depends on frames real-time | Pre-render for large |
| Pairplot | O(n² × p²) | O(n × p) | 5K samples × 10 vars | Large: use `sns.pairplot(sample=1000)` |

**Optimization tips**:
- Downsample before plotting if n > 10,000 points (scatter)
- Use hexbin or 2D histogram for dense scatter plots
- For large heatmaps, use `sparse=True` or clip values
- Network: use `layout='fa'` (forceatlas2) for large graphs
- Animation: pre-render frames, don't real-time compute

## 7. Styling & Customization

### Theme Selection

| Theme | Use Case | Colors | Font |
|-------|----------|--------|------|
| `publication` | Journal figures | Colorblind-safe (Set2) | Helvetica 8pt |
| `presentation` | Slides/poster | High-contrast (tab10) | Arial 14pt bold |
| `default` | General use | Matplotlib default | DejaVu Sans |
| `dark` | Dark background talks | Neon/bright | Sans-serif white |

```python-snippet
from metainformant.visualization.styling import (
    set_theme,
    set_color_palette,
    get_palette
)

set_theme('publication')
set_color_palette('colorblind')  # Okabe-Ito scale
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
from metainformant.visualization.dashboards import create_figure_grid
from metainformant.visualization.genomics import manhattan_plot, qq_plot
from metainformant.visualization.statistical import histogram

# Create 2×2 figure grid
fig = create_figure_grid(
    plots=[
        (manhattan_plot, {"data": gwas_df}, (0, 0)),      # top-left
        (qq_plot, {"pvalues": gwas_df['p']}, (0, 1)),    # top-right
        (histogram, {"data": gwas_df['beta']}, (1, 0)),  # bottom-left
        (volcano_plot, {"de": de_df}, (1, 1)),           # bottom-right
    ],
    shared_legend=False,
    figsize=(12, 10),
    hspace=0.3,
    wspace=0.2
)

fig.save("gwas_manuscript_figure.pdf", dpi=300)
```

## 9. Cross-Module Visualization Usage Examples

### Equivalence: GWAS Plotting Across Modules

```python-snippet
# Option 1: Direct from gwas module (has built-in viz)
from metainformant.gwas.visualization import plot_manhattan, plot_qq
plot_manhattan(gwas_results)
plot_qq(gwas_results['p_value'])

# Option 2: Using unified visualization module
from metainformant.visualization.genomics import manhattan_plot, qq_plot
manhattan_plot(
    data=gwas_results,
    position_col='pos',
    pval_col='p_value',
    chromosome_col='chrom'
)
```

### Equivalence: Expression Visualization

```python-snippet
# From rna module
from metainformant.rna import analysis as rna_analysis
rna_analysis.plot_expression_heatmap(expression_matrix)

# From visualization module (more general)
from metainformant.visualization.expression import expression_heatmap
expression_heatmap(
    expression_matrix,
    annotation=metadata,
    cluster_rows=True,
    cluster_cols=True
)
```

## 10. Advanced: Custom Plot Development

If existing plots don't fit your needs, extend via `visualization.plots`:

```python-snippet
from metainformant.visualization.plots.base import BasePlotter
import matplotlib.pyplot as plt

class CustomManhattanPlotter(BasePlotter):
    def __init__(self, highlight_genes=None):
        super().__init__()
        self.highlight_genes = highlight_genes or []
    
    def plot(self, data, ax=None):
        ax = ax or plt.gca()
        # Custom logic
        for chrom in data['chrom'].unique():
            chrom_data = data[data['chrom'] == chrom]
            ax.scatter(chrom_data['pos'], -np.log10(chrom_data['p']))
        # Highlight specific genes
        for gene in self.highlight_genes:
            gene_pos = data[data['gene'] == gene]['pos'].iloc[0]
            ax.axvline(gene_pos, color='red', linestyle='--')
        return ax

# Use it
plotter = CustomManhattanPlotter(highlight_genes=['APOE', 'CLU'])
plotter.plot(gwas_df)
```

## 11. Summary Table: All Visualization Submodules

| Submodule | Plot Count | Primary Domains | Key Functions | Best For |
|-----------|------------|-----------------|---------------|----------|
| `basic/` | 12 | General | line, scatter, bar, pie, area, step, heatmap | Quick exploratory plots |
| `statistical/` | 10 | Stats | histogram, boxplot, violin, qq, roc, density | Statistical diagnostics |
| `genomics/` | 8 | GWAS, genomics | manhattan, volcano, regional, ideogram, circos | Variant/expression studies |
| `expression/` | 6 | RNA-seq, proteomics | expression_heatmap, ma_plot, enrichment | Gene/protein expression |
| `dimred/` | 6 | ML, scRNA | pca, umap, tsne, scree | Dimensionality reduction viz |
| `networks/` | 5 | PPI, pathways | plot_network, plot_pathway | Biological networks |
| `trees/` | 4 | Phylogeny | plot_phylo_tree, plot_tree_circular | Evolutionary trees |
| `timeseries/` | 5 | Temporal | line, forecast, autocorrelation, seasonal | Time-series data |
| `multidim/` | 4 | High-dim | pairplot, parallel_coords, radar | Multi-variable relationships |
| `quality/` | 4 | QC | quality_scores, adapter_content, length_dist | Sequencing QC |
| `information/` | 3 | Info theory | entropy_plot, mutual_information | Info theoretic analysis |
| `animations/` | 6 | Dynamic | animate_time_series, animate_evolution | Animated figures |
| `dashboards/` | N/A | Composite | create_dashboard, multi_panel | Multi-plot figures |
| `styling/` | N/A | Theming | set_theme, set_palette | Figure aesthetics |

**Total**: 70+ individual plot functions across all submodules.

## 12. Related Documentation

- **[visualization/index.md](../visualization/index.md)** — Visualization module overview
- **[visualization/README.md](../visualization/README.md)** — Full function reference
- **[visualization/genomics.md](../visualization/genomics.md)** — Genomics-specific plots
- **[visualization/statistical.md](../visualization/statistical.md)** — Statistical plots
- **[visualization/trees.md](../visualization/trees.md)** — Phylogenetic tree visualization
- **[visualization/animations.md](../visualization/animations.md)** — Animation guide
- **[COMPARISON_GUIDES.md](../COMPARISON_GUIDES.md)** — Master comparison index
