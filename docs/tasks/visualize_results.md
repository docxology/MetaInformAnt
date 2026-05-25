# Visualization Quick Reference

Common plots and how to generate them across all METAINFORMANT modules.

## When to Use

Use `visualize_results` for publication-ready figures and exploratory data analysis—not for interactive dashboards (use `dashboards/` module) or static report generation (use `multiqc/`).

## Table of Contents

- [Manhattan Plot (GWAS)](#manhattan-plot-gwas)
- [PCA Plot (Multi-omics)](#pca-plot-multi-omics)
- [Heatmap (DEGs)](#heatmap-degs)
- [Phylogenetic Tree](#phylogenetic-tree)
- [Network Graph](#network-graph)
- [Advanced Examples](#advanced-examples)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

```python
from metainformant.visualization import plots
import pandas as pd

gwas = pd.read_csv("gwas_results.tsv", sep='	')
fig = plots.manhattan(
    gwas,
    genome="amellifera",  # or provide chr_lengths dict
    highlight_genes=["Amel_004123", "Amel_007890"],
    suggestiveline=-np.log10(1e-5),
    genomewideline=-np.log10(5e-8)
)
fig.savefig("manhattan.png", dpi=300, bbox_inches='tight')
```

## PCA Plot (Multi-omics)

```python-snippet
from metainformant.visualization import pca_plots

# From expression matrix
expr = pd.read_csv("counts.tsv", sep='	', index_col=0)
fig = pca_plots.samples(expr, color_by="tissue", shape_by="individual")
fig.savefig("pca_by_tissue.png")
```

## Heatmap (DEGs)

```python-snippet
from metainformant.visualization import heatmaps

deseq = pd.read_csv("deseq2_results.tsv", sep='	')
top100 = deseq.nsmallest(100, "padj")

counts = pd.read_csv("normalized_counts.tsv", sep='	', index_col=0)
counts_top100 = counts.loc[top100["gene_id"]]

fig = heatmaps.expression(
    counts_top100,
    z_score="row",
    cluster_rows=True,
    cluster_cols=True,
    annotation=pd.DataFrame({"tissue": sample_meta["tissue"]})
)
fig.savefig("heatmap_top100_degs.png")
```

## Phylogenetic Tree

```python
from metainformant.dna import phylogeny

tree = phylogeny.read_newick("tree.nwk")
fig = phylogeny.plot(
    tree,
    show_branch_support=True,
    label_leaves=True,
    tip_labels_replace={"_": " "}
)
fig.savefig("phylogeny.pdf")
```

## Network Graph

```python-snippet
from metainformant.networks import graph, visualize

ppi = graph.read_edgelist("string_ppi.tsv")
fig = visualize.spring_layout(
    ppi,
    node_size="degree",
    node_color="community",
    highlight_nodes=["TP53", "BRCA1"]
)
fig.savefig("ppi_network.png")
```

## Advanced Examples

### Circos plot for multi-omics integration
```python-snippet
from metainformant.visualization import circos

fig = circos.genome_view(
    genome="amellifera",
    tracks={
        "GWAS": {"file": "gwas_hits.tsv", "type": "manhattan"},
        "RNA": {"file": "deseq2_results.tsv", "type": "volcano"},
        "Methylation": {"file": "meth_delta.tsv", "type": "heatmap"}
    },
    highlight_genes=["Amel_004123", "Amel_007890"]
)
fig.savefig("multiomics_circos.pdf", dpi=300)
```

### Interactive Plotly dashboard (standalone HTML)
```python
from metainformant.visualization.dashboards import interactive

dashboard = interactive.create_dashboard(
    data={
        "expression": "counts.tsv",
        "gwas": "gwas_results.tsv",
        "metadata": "sample_info.tsv"
    },
    layout="grid",
    widgets=["manhattan", "pca", "heatmap"]
)
dashboard.save("explorer.html")
# Opens in browser: file:///path/to/explorer.html
```

### Animated trajectory plot (single-cell pseudotime)
```python-snippet
from metainformant.visualization import animation

# Requires monotonic ordering of cells (e.g., from Slingshot/palantir)
cells = pd.read_csv("pseudotime.tsv")

fig = animation.trajectory(
    embedding=cells[["UMAP1", "UMAP2"]].values,
    pseudotime=cells["pseudotime"],
    color_by=cells["cluster"],
    fps=30,
    duration_secs=10
)
fig.save("trajectory.mp4")
```

### 3D protein structure with PyMOL rendering
```python-snippet
from metainformant.protein import visualize as prot_viz

# Load AlphaFold structure
structure = prot_viz.load_pdb("AF-P12345-F1-model_v4.pdb")

# Highlight domains and mutations
fig = prot_viz.render(
    structure,
    highlight_residues=[45, 128, 256],
    show_surface=True,
    color_by="B-factor"  # Flexibility
)
fig.save("protein_structure.png", width=1920, height=1080)
```

## Expected Output

### Manhattan plot (file)
```
File: manhattan.png
Size: 2.4 MB
Dimensions: 2400×1200 px (300 DPI)
Genome: Apis mellifera (chr1-chr16, chrX)
Significant hits (p < 5e-8): 23 loci
Suggestive hits (p < 1e-5): 142 SNPs
Annotated genes: AMEL_004123 (chr6:28.5 Mb), AMEL_007890 (chr11:12.1 Mb)
```

### PCA plot caption-ready data preview
```
> head pca.tsv
sample         PC1      PC2      PC3      tissue    individual
S1_R1          12.34    -5.67    2.89     brain     bee_001
S2_R1          11.89    -5.23    3.01     brain     bee_001
S3_R1         -12.45     8.91   -2.34     ovary     bee_045
...
Variance explained: PC1=18.2%, PC2=8.7%, PC3=4.1%
Color scheme: tissue (brain=blue, ovary=red, leg=green)
```

### Heatmap with dendrogram (PDF output)
```
File: heatmap_top100_degs.pdf
Genes: 100 (top by padj)
Samples: 48 (ordered by tissue cluster)
Z-score: row-scaling (mean=0, sd=1)
Clustering: Ward's method, Euclidean distance
Dendrogram leaves: brain samples cluster tightly (bootstrap >95%)
```

### Network graph statistics
```
Graph: STRING PPI (score > 0.7)
Nodes: 1,234 proteins
Edges: 4,567 interactions
Diameter: 8 hops (TP53 → BRCA1)
Clustering coefficient: 0.42
Top 5 hubs: TP53 (deg=156), BRCA1 (deg=142), EGFR (deg=134), ...
Community structure: 7 modules (Louvain, resolution=1.0)
```

## Common Pitfalls

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| `MemoryError` during heatmap with 50k genes | Distance matrix O(n²) for clustering | Subset genes first (`deseq.top_n(200)`), or use `cluster_rows=False` |
| Manhattan plot x-axis labels overlap | Many small chromosomes or custom genome | Use `rotate_xlabels=45` or `chrom_label_step=5` (label every 5th) |
| Network graph layout takes forever | >10k nodes without fast layout algorithm | Switch to `kamada_kawai` (O(n³) but better for small) → actually use `sfdp` or `sparse` from graph-tool; or sample subgraph |
| `ValueError: color mapping size mismatch` | `color_by` column has more unique values than palette | Set `palette="tab20"` for >10 categories, or provide `palette={...}` dict |
| Plot saves blank/empty file | Matplotlib backend not set for non-interactive | Add `import matplotlib; matplotlib.use('Agg')` at top of script |
| Font missing warnings (CJK/Russian) | System lacks Unicode fonts | Install `fonts-noto` (Linux) or set `fontfamily="DejaVu Sans"` |

---

**Related:** [Full visualization docs](../visualization/index.md) | [Plot library API](../visualization/plots.md) | [Gallery](../visualization/gallery.md) | [Quality control plots](../visualization/quality.md)
