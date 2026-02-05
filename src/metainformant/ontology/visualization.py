"""Ontology visualization functions.

This module provides comprehensive visualization capabilities for ontologies,
including gene ontology (GO) trees, semantic similarity matrices, enrichment plots,
and functional annotation networks.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

from metainformant.core import logging, paths, validation

# Optional scientific dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    nx = None

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_go_dag(
    go_graph: Any,  # nx.DiGraph
    terms: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot Gene Ontology directed acyclic graph (DAG).

    Args:
        go_graph: NetworkX DiGraph representing GO hierarchy
        terms: Optional list of GO terms to highlight
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx is required for GO DAG plotting. " "Install with: uv pip install networkx")

    validation.validate_type(go_graph, nx.DiGraph, "go_graph")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Use hierarchical layout
    try:
        pos = nx.spring_layout(go_graph, k=1, iterations=50)
    except (ValueError, nx.NetworkXError):
        pos = nx.random_layout(go_graph)

    # Draw nodes
    node_colors = ["red" if terms and node in terms else "lightblue" for node in go_graph.nodes()]

    nx.draw_networkx_nodes(go_graph, pos, node_color=node_colors, node_size=300, alpha=0.8, ax=ax)

    # Draw edges
    nx.draw_networkx_edges(go_graph, pos, edge_color="gray", arrows=True, arrowsize=10, ax=ax)

    # Draw labels (only for highlighted terms or small graphs)
    if terms or len(go_graph) < 20:
        labels = {node: node for node in (terms if terms else go_graph.nodes())}
        nx.draw_networkx_labels(go_graph, pos, labels, font_size=8, ax=ax)

    ax.set_title("Gene Ontology DAG")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"GO DAG plot saved to {output_path}")

    return ax


def plot_semantic_similarity_matrix(
    similarity_matrix: np.ndarray,
    term_labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot semantic similarity matrix heatmap.

    Args:
        similarity_matrix: Square matrix of semantic similarities
        term_labels: Optional labels for matrix axes
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(similarity_matrix, np.ndarray, "similarity_matrix")

    if similarity_matrix.shape[0] != similarity_matrix.shape[1]:
        raise ValueError("Similarity matrix must be square")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    if HAS_SEABORN:
        mask = np.triu(np.ones_like(similarity_matrix, dtype=bool))
        sns.heatmap(
            similarity_matrix,
            mask=mask,
            annot=len(similarity_matrix) <= 10,
            fmt=".2f",
            cmap="YlOrRd",
            square=True,
            ax=ax,
            xticklabels=term_labels,
            yticklabels=term_labels,
            **kwargs,
        )
    else:
        im = ax.imshow(similarity_matrix, cmap="YlOrRd", aspect="equal")
        ax.set_xticks(range(len(similarity_matrix)))
        ax.set_yticks(range(len(similarity_matrix)))
        if term_labels:
            ax.set_xticklabels(term_labels, rotation=45, ha="right")
            ax.set_yticklabels(term_labels)
        plt.colorbar(im, ax=ax)

    ax.set_title("Semantic Similarity Matrix")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Semantic similarity matrix saved to {output_path}")

    return ax


def plot_go_enrichment_barplot(
    enrichment_results: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot GO enrichment results as a bar plot.

    Args:
        enrichment_results: List of enrichment result dictionaries
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(enrichment_results, list, "enrichment_results")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract data
    terms = [result.get("term", f"Term {i}") for i, result in enumerate(enrichment_results)]
    pvalues = [-np.log10(result.get("pvalue", 1.0)) for result in enrichment_results]
    fold_changes = [result.get("fold_change", 1.0) for result in enrichment_results]

    # Sort by significance
    sorted_idx = np.argsort(pvalues)[::-1]
    terms = [terms[i] for i in sorted_idx]
    pvalues = [pvalues[i] for i in sorted_idx]

    # Plot horizontal bar chart
    bars = ax.barh(range(len(terms)), pvalues, color="skyblue", alpha=0.8)

    # Add term labels
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms)
    ax.set_xlabel("-log₁₀(p-value)")
    ax.set_title("GO Enrichment Analysis")

    # Add value labels on bars
    for i, (bar, pval) in enumerate(zip(bars, pvalues)):
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2, ".2f", ha="left", va="center", fontsize=8)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"GO enrichment barplot saved to {output_path}")

    return ax


def plot_go_enrichment_dotplot(
    enrichment_results: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot GO enrichment results as a dot plot (bubble plot).

    Args:
        enrichment_results: List of enrichment result dictionaries
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(enrichment_results, list, "enrichment_results")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract data
    terms = [result.get("term", f"Term {i}") for i, result in enumerate(enrichment_results)]
    pvalues = [-np.log10(result.get("pvalue", 1.0)) for result in enrichment_results]
    gene_counts = [result.get("gene_count", 1) for result in enrichment_results]

    # Create scatter plot with bubble sizes
    scatter = ax.scatter(
        pvalues, range(len(terms)), s=np.array(gene_counts) * 20, c=pvalues, cmap="Reds", alpha=0.7, edgecolors="black"
    )

    # Add term labels
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms)
    ax.set_xlabel("-log₁₀(p-value)")
    ax.set_title("GO Enrichment Analysis (Dot Plot)")

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("-log₁₀(p-value)")

    # Add legend for bubble sizes
    sizes = [5, 10, 20, 30]
    labels = [f"{s//20} genes" for s in sizes]
    legend_elements = [
        plt.scatter([], [], s=s, c="gray", alpha=0.5, edgecolors="black", label=label)
        for s, label in zip(sizes, labels)
    ]
    ax.legend(handles=legend_elements, title="Gene Count", bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"GO enrichment dotplot saved to {output_path}")

    return ax


def plot_ontology_network(
    ontology_graph: Any,  # nx.Graph
    node_colors: Dict[str, str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot ontology relationship network.

    Args:
        ontology_graph: NetworkX graph of ontology relationships
        node_colors: Optional mapping of node names to colors
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError(
            "networkx is required for ontology network plotting. " "Install with: uv pip install networkx"
        )

    validation.validate_type(ontology_graph, nx.Graph, "ontology_graph")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Layout
    pos = nx.spring_layout(ontology_graph, k=1, iterations=50)

    # Node colors
    if node_colors:
        colors = [node_colors.get(node, "lightblue") for node in ontology_graph.nodes()]
    else:
        colors = "lightblue"

    # Draw network
    nx.draw_networkx_nodes(ontology_graph, pos, node_color=colors, node_size=300, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(ontology_graph, pos, edge_color="gray", alpha=0.5, ax=ax)
    nx.draw_networkx_labels(ontology_graph, pos, font_size=8, ax=ax)

    ax.set_title("Ontology Relationship Network")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Ontology network plot saved to {output_path}")

    return ax


def plot_information_content_profile(
    term_ic: Dict[str, float],
    sorted_terms: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot information content profile for ontology terms.

    Args:
        term_ic: Dictionary mapping terms to information content values
        sorted_terms: Optional ordered list of terms to plot
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(term_ic, dict, "term_ic")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Sort terms by information content if not provided
    if sorted_terms is None:
        sorted_terms = sorted(term_ic.keys(), key=lambda x: term_ic[x], reverse=True)

    ic_values = [term_ic[term] for term in sorted_terms]

    # Plot
    ax.plot(range(len(sorted_terms)), ic_values, "b-", linewidth=2, marker="o", markersize=4)
    ax.fill_between(range(len(sorted_terms)), ic_values, alpha=0.3, color="blue")

    # Add term labels for top terms
    n_labels = min(10, len(sorted_terms))
    step = len(sorted_terms) // n_labels
    for i in range(0, len(sorted_terms), step):
        ax.annotate(
            sorted_terms[i], (i, ic_values[i]), xytext=(5, 5), textcoords="offset points", fontsize=8, rotation=45
        )

    ax.set_xlabel("Terms (sorted by IC)")
    ax.set_ylabel("Information Content")
    ax.set_title("Ontology Term Information Content")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Information content profile saved to {output_path}")

    return ax


def plot_go_term_hierarchy(
    go_graph: Any,  # nx.DiGraph
    root_term: str,
    max_depth: int = 3,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot GO term hierarchy from a root term.

    Args:
        go_graph: NetworkX DiGraph representing GO hierarchy
        root_term: Root GO term ID
        max_depth: Maximum depth to traverse
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError(
            "networkx is required for GO term hierarchy plotting. " "Install with: uv pip install networkx"
        )

    validation.validate_type(go_graph, nx.DiGraph, "go_graph")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract subgraph up to max_depth
    nodes_to_include = {root_term}
    current_level = {root_term}

    for depth in range(max_depth):
        next_level = set()
        for node in current_level:
            next_level.update(go_graph.successors(node))
        nodes_to_include.update(next_level)
        current_level = next_level

    subgraph = go_graph.subgraph(nodes_to_include)

    # Hierarchical layout
    def hierarchy_pos(G, root, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5):
        pos = {root: (xcenter, vert_loc)}
        if len(G) == 1:
            return pos

        children = list(G.successors(root))
        if len(children) != 0:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos.update(
                    hierarchy_pos(G, child, width=dx, vert_gap=vert_gap, vert_loc=vert_loc - vert_gap, xcenter=nextx)
                )
        return pos

    pos = hierarchy_pos(subgraph, root_term)

    # Draw
    nx.draw_networkx_nodes(subgraph, pos, node_color="lightblue", node_size=500, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(subgraph, pos, edge_color="gray", arrows=True, arrowsize=15, ax=ax)
    nx.draw_networkx_labels(subgraph, pos, font_size=8, ax=ax)

    ax.set_title(f"GO Term Hierarchy from {root_term}")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"GO term hierarchy plot saved to {output_path}")

    return ax


def plot_functional_annotation_heatmap(
    annotation_matrix: np.ndarray,
    term_labels: List[str],
    gene_labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot functional annotation matrix as heatmap.

    Args:
        annotation_matrix: Binary matrix of gene-term annotations
        term_labels: Labels for GO terms
        gene_labels: Optional labels for genes
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(annotation_matrix, np.ndarray, "annotation_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot binary heatmap
    im = ax.imshow(annotation_matrix, cmap="Blues", aspect="auto", origin="lower")

    ax.set_xlabel("GO Terms")
    ax.set_ylabel("Genes")
    ax.set_title("Functional Annotation Matrix")

    # Set tick labels if provided
    if gene_labels and len(gene_labels) <= 50:  # Only show labels for reasonable sizes
        ax.set_yticks(range(len(gene_labels)))
        ax.set_yticklabels(gene_labels)
    else:
        ax.set_ylabel(f"Genes (n={annotation_matrix.shape[0]})")

    if len(term_labels) <= 50:
        ax.set_xticks(range(len(term_labels)))
        ax.set_xticklabels(term_labels, rotation=45, ha="right")
    else:
        ax.set_xlabel(f"GO Terms (n={annotation_matrix.shape[1]})")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Annotated (1) / Not Annotated (0)")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Functional annotation heatmap saved to {output_path}")

    return ax


def plot_semantic_similarity_clustermap(
    similarity_matrix: np.ndarray,
    term_labels: List[str] | None = None,
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Any:
    """Plot semantic similarity clustermap.

    Args:
        similarity_matrix: Square matrix of semantic similarities
        term_labels: Optional labels for matrix axes
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        seaborn ClusterGrid object (or matplotlib figure)
    """
    if not HAS_SEABORN:
        raise ImportError("seaborn required for clustermap visualization")

    validation.validate_type(similarity_matrix, np.ndarray, "similarity_matrix")

    if similarity_matrix.shape[0] != similarity_matrix.shape[1]:
        raise ValueError("Similarity matrix must be square")

    # Create clustermap
    if term_labels:
        g = sns.clustermap(
            similarity_matrix,
            xticklabels=term_labels,
            yticklabels=term_labels,
            cmap="YlOrRd",
            figsize=figsize,
            **kwargs,
        )
    else:
        g = sns.clustermap(similarity_matrix, cmap="YlOrRd", figsize=figsize, **kwargs)

    plt.suptitle("Semantic Similarity Clustermap", y=1.02)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Semantic similarity clustermap saved to {output_path}")

    return g


def create_interactive_go_network(
    go_graph: nx.DiGraph,
    term_annotations: Dict[str, Any] | None = None,
    *,
    output_path: str | Path | None = None,
    **kwargs,
) -> Any:
    """Create an interactive GO network visualization using Plotly.

    Args:
        go_graph: NetworkX DiGraph representing GO hierarchy
        term_annotations: Optional annotations for terms
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive GO network")

    validation.validate_type(go_graph, nx.DiGraph, "go_graph")

    # Create positions
    pos = nx.spring_layout(go_graph, dim=3, k=1, iterations=50)

    # Extract node and edge data
    node_x = []
    node_y = []
    node_z = []
    node_text = []

    for node in go_graph.nodes():
        x, y, z = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)

        # Create hover text
        text = f"GO Term: {node}"
        if term_annotations and node in term_annotations:
            text += f"<br>{term_annotations[node]}"
        node_text.append(text)

    # Create edges
    edge_x = []
    edge_y = []
    edge_z = []

    for edge in go_graph.edges():
        x0, y0, z0 = pos[edge[0]]
        x1, y1, z1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])

    # Create figure
    fig = go.Figure()

    # Add edges
    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="gray", width=1),
            hoverinfo="none",
            name="Relationships",
        )
    )

    # Add nodes
    fig.add_trace(
        go.Scatter3d(
            x=node_x,
            y=node_y,
            z=node_z,
            mode="markers",
            marker=dict(size=6, color="lightblue", line=dict(width=1, color="black")),
            text=node_text,
            hovertemplate="%{text}<extra></extra>",
            name="GO Terms",
        )
    )

    # Update layout
    fig.update_layout(
        title="Interactive GO Network",
        scene=dict(
            xaxis=dict(showticklabels=False),
            yaxis=dict(showticklabels=False),
            zaxis=dict(showticklabels=False),
        ),
        **kwargs,
    )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive GO network saved to {html_path}")

    return fig
