# Visualization Quick Reference

Common plots and how to generate them across all METAINFORMANT modules.

## When to Use

Use `visualize_results` for publication-ready figures and exploratory data analysis—not for interactive dashboards (use `dashboards/` module) or static report generation (use `multiqc/`).

## Table of Contents

- [Manhattan Plot](#manhattan-plot)
- [Expected Manhattan Plot Output](#manhattan-plot-file)
- [PCA Plot (Multi-omics)](#pca-plot-multi-omics)
- [Heatmap (DEGs)](#heatmap-degs)
- [Phylogenetic Tree](#phylogenetic-tree)
- [Network Graph](#network-graph)
- [Advanced Examples](#advanced-examples)
- [Expected Output](#expected-output)
- [Common Pitfalls](#common-pitfalls)

---

## Manhattan Plot

```python
import pandas as pd
from metainformant.gwas.visualization import manhattan_plot

# TSV with one row per variant and chrom, pos and p_value columns
gwas = pd.read_csv("gwas_results.tsv", sep="\t")
fig = manhattan_plot(
    gwas.to_dict("records"),  # the plot functions take a list of per-variant dictionaries
    output_path="manhattan.png",  # the output directory must already exist
    significance_threshold=5e-8,  # genome-wide line
    suggestive_threshold=1e-5,
    gene_annotations=[  # an arrow and label at the nearest variant
        {"chrom": "6", "pos": 28_500_000, "gene_name": "Amel_004123"},
        {"chrom": "11", "pos": 12_100_000, "gene_name": "Amel_007890"},
    ],
)
```

## PCA Plot (Multi-omics)

```python-snippet
import pandas as pd
from sklearn.decomposition import PCA
from metainformant.visualization.plots.general import pca_plot

# From expression matrix (samples x genes) and sample metadata
expr = pd.read_csv("counts.tsv", sep='\t', index_col=0)
sample_meta = pd.read_csv("sample_info.tsv", sep='\t', index_col=0)

pca = PCA(n_components=2).fit(expr)
scores = pd.DataFrame(pca.transform(expr), index=expr.index, columns=["PC1", "PC2"])
scores["PC1_variance"], scores["PC2_variance"] = pca.explained_variance_ratio_ * 100
scores["tissue"] = sample_meta["tissue"]

ax = pca_plot(scores, hue="tissue", output_path="pca_by_tissue.png")
```

## Heatmap (DEGs)

```python-snippet
import pandas as pd
from metainformant.visualization.plots.general import expression_heatmap

deseq = pd.read_csv("deseq2_results.tsv", sep='\t')
top100 = deseq.nsmallest(100, "padj")

counts = pd.read_csv("normalized_counts.tsv", sep='\t', index_col=0)
counts_top100 = counts.loc[top100["gene_id"]]

# Row z-score (mean=0, sd=1 per gene)
z = counts_top100.sub(counts_top100.mean(axis=1), axis=0).div(counts_top100.std(axis=1), axis=0)

ax = expression_heatmap(z, vmin=-3, vmax=3, cbar_label="Z-score", output_path="heatmap_top100_degs.png")
```

## Phylogenetic Tree

```python
from Bio import Phylo
from metainformant.visualization.genomics.trees import plot_phylo_tree

# Read a Newick file with Bio.Phylo. A tree built with metainformant.dna.phylogeny converts the same way:
# Phylo.read(StringIO(phylogeny.to_newick(tree)), "newick")
tree = Phylo.read("tree.nwk", "newick")
ax = plot_phylo_tree(tree, with_labels=True, output_path="phylogeny.pdf")
# Other layouts: circular_tree_plot, unrooted_tree_plot, tree_annotation_plot, tree_comparison_plot
```

## Network Graph

```python-snippet
from metainformant.networks.analysis.community import detect_communities
from metainformant.networks.analysis.graph import load_network
from metainformant.visualization.genomics.networks import plot_community_network

ppi = load_network("string_ppi.tsv", format="edgelist")
communities = detect_communities(ppi)  # node -> community id (Louvain by default)
ax = plot_community_network(ppi, communities, output_path="ppi_network.png")
```

## Advanced Examples

### Interactive Plotly plots (standalone HTML)
```python
import pandas as pd
from metainformant.visualization.dashboards import interactive

# With Plotly, each call returns a Plotly figure and writes a standalone HTML file (a suffix other than .html is replaced).
# Without Plotly an .html output_path raises ImportError; an image path such as .png returns a static matplotlib Axes.
gwas = pd.read_csv("gwas_results.tsv", sep="\t")  # chrom, pos and p_value columns
interactive.interactive_manhattan(
    gwas, chrom_col="chrom", pos_col="pos", pvalue_col="p_value", output_path="manhattan.html"
)

deseq = pd.read_csv("deseq2_results.tsv", sep="\t")  # gene, log2FoldChange and pvalue columns
interactive.interactive_volcano(deseq, output_path="volcano.html")

counts = pd.read_csv("normalized_counts.tsv", sep="\t", index_col=0)
interactive.interactive_heatmap(counts.iloc[:100], output_path="heatmap.html")
# Opens in browser: file:///path/to/manhattan.html
```

### Animated trajectory plot (single-cell pseudotime)
```python-snippet
import pandas as pd
from metainformant.visualization.plots.animations import animate_trajectory

# Requires monotonic ordering of cells (e.g., from Slingshot/palantir)
cells = pd.read_csv("pseudotime.tsv", sep='\t').sort_values("pseudotime")

# One trajectory = the cells' UMAP coordinates ordered along pseudotime
fig, anim = animate_trajectory(
    [cells[["UMAP1", "UMAP2"]].values],
    interval=100,             # ms per frame
    output_path="trajectory.gif",
)
```

### 3D protein structure (interactive Plotly viewer)
```python-snippet
from pathlib import Path
from metainformant.protein.structure.pdb import load_pdb_file
from metainformant.protein.visualization.general import create_interactive_structure_viewer

# Load AlphaFold structure
structure = load_pdb_file(Path("AF-P12345-F1-model_v4.pdb"))

# Interactive 3D atom plot saved as HTML (requires plotly)
fig = create_interactive_structure_viewer(structure, output_path="protein_structure.html")
```

## Expected Output

### Manhattan plot (file)
```
File: manhattan.png
Dimensions: about 2900×1700 px (300 DPI; the default 12×6 inch figure, cropped to its contents)
Title: Manhattan Plot (n=600 variants, 1 GW-significant)
Lines: genome-wide significance at p = 5e-8, suggestive at p = 1e-5
Annotated genes: Amel_004123 (chr6:28.5 Mb), Amel_007890 (chr11:12.1 Mb)
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

### Heatmap (PNG output)
```
File: heatmap_top100_degs.png
Genes: 100 (top by padj)
Samples: 48
Z-score: row-scaling (mean=0, sd=1)
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
| `MemoryError` during heatmap with 50k genes | Distance matrix O(n²) for clustering | Subset genes first (`deseq.nsmallest(200, "padj")`) |
| Manhattan plot x-axis labels overlap | Many small chromosomes or contigs | Pass `figsize=(16, 6)`, or draw on your own axes (`manhattan_plot(results, ax=ax)`) and rotate the labels with `ax.tick_params(axis="x", labelrotation=45)` |
| Network graph layout takes forever | More than ~10k nodes | Plot a subgraph instead (`metainformant.networks.analysis.graph.create_subgraph(graph, nodes)`), for example the largest communities |
| `ValueError: Missing required columns for PCA plot` | The score table lacks `PC1`, `PC2`, `PC1_variance` or `PC2_variance` | Add the columns as in the PCA example (variance as percentages) |
| `FileNotFoundError` when saving a GWAS figure | `manhattan_plot`, `qq_plot` and `regional_plot` do not create missing directories | Create the output directory first (`metainformant.core.io.ensure_directory(...)`) |
| `ImportError: Plotly is required to write the interactive HTML file ...` from an `interactive_*` function in `metainformant.visualization.dashboards.interactive` | Plotly is not installed and `output_path` ends in `.html` or `.htm`; matplotlib cannot write HTML | Install Plotly (`uv sync --extra visualization`), or pass an image path such as `.png` to get the static fallback |
| Plot saves blank/empty file | Matplotlib backend not set for non-interactive | Add `import matplotlib; matplotlib.use('Agg')` at top of script |
| Font missing warnings (CJK/Russian) | System lacks Unicode fonts | Install `fonts-noto` (Linux) or set `plt.rcParams["font.family"] = "DejaVu Sans"` |

---

**Related:** [Full visualization docs](../visualization/index.md) | [Plot library API](../visualization/plots.md) | [Gallery](../visualization/gallery.md) | [Quality control plots](../visualization/quality.md)
