"""Multi-omics data integration and visualization functions.

This module provides comprehensive visualization capabilities for integrated
multi-omics data analysis, including cross-platform correlation plots,
integrated heatmaps, pathway enrichment visualizations, and systems-level
analysis results.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, Rectangle

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

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
    from plotly.subplots import make_subplots

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None
    make_subplots = None


def plot_multiomics_correlation_heatmap(
    correlation_matrix: np.ndarray,
    omics_names: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 10),
    **kwargs,
) -> Axes:
    """Plot correlation heatmap across multiple omics layers.

    Args:
        correlation_matrix: Correlation matrix between all omics features
        omics_names: Names of omics types/layers
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(correlation_matrix, np.ndarray, "correlation_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot correlation heatmap
    if HAS_SEABORN:
        sns.heatmap(
            correlation_matrix,
            annot=len(correlation_matrix) <= 20,
            fmt=".2f",
            cmap="coolwarm",
            center=0,
            ax=ax,
            xticklabels=omics_names,
            yticklabels=omics_names,
            **kwargs,
        )
    else:
        im = ax.imshow(correlation_matrix, cmap="coolwarm", aspect="equal", vmin=-1, vmax=1)
        ax.set_xticks(range(len(omics_names)))
        ax.set_yticks(range(len(omics_names)))
        ax.set_xticklabels(omics_names, rotation=45, ha="right")
        ax.set_yticklabels(omics_names)
        plt.colorbar(im, ax=ax)

    ax.set_title("Multi-Omics Correlation Matrix")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics correlation heatmap saved to {output_path}")

    return ax


def plot_integrated_omics_heatmap(
    omics_data: Dict[str, np.ndarray],
    feature_names: List[str] | None = None,
    sample_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (14, 10),
    **kwargs,
) -> Axes:
    """Plot integrated heatmap showing multiple omics layers.

    Args:
        omics_data: Dictionary mapping omics types to data matrices
        feature_names: Optional names for features
        sample_names: Optional names for samples
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(omics_data, dict, "omics_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Combine all omics data
    combined_data = []
    omics_labels = []
    feature_labels = []

    for omics_type, data in omics_data.items():
        combined_data.append(data)
        omics_labels.extend([omics_type] * data.shape[0])
        if feature_names:
            feature_labels.extend(feature_names[: data.shape[0]])

    combined_matrix = np.vstack(combined_data) if combined_data else np.array([])

    # Plot integrated heatmap
    if HAS_SEABORN:
        ax = sns.heatmap(
            combined_matrix, cmap="RdYlBu_r", ax=ax, xticklabels=sample_names, yticklabels=feature_labels, **kwargs
        )
    else:
        im = ax.imshow(combined_matrix, cmap="RdYlBu_r", aspect="auto")
        if sample_names:
            ax.set_xticks(range(len(sample_names)))
            ax.set_xticklabels(sample_names, rotation=45, ha="right")
        plt.colorbar(im, ax=ax)

    ax.set_title("Integrated Multi-Omics Heatmap")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Integrated omics heatmap saved to {output_path}")

    return ax


def plot_omics_layer_comparison(
    omics_datasets: Dict[str, np.ndarray],
    comparison_metric: str = "correlation",
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot comparison between different omics layers.

    Args:
        omics_datasets: Dictionary mapping omics types to data matrices
        comparison_metric: Metric to use for comparison ('correlation', 'euclidean', 'cosine')
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(omics_datasets, dict, "omics_datasets")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    omics_types = list(omics_datasets.keys())
    n_omics = len(omics_types)

    # Calculate pairwise comparisons
    comparison_matrix = np.zeros((n_omics, n_omics))

    for i in range(n_omics):
        for j in range(n_omics):
            data1 = omics_datasets[omics_types[i]]
            data2 = omics_datasets[omics_types[j]]

            if comparison_metric == "correlation":
                # Calculate mean correlation between features
                corr_matrix = np.corrcoef(data1.T, data2.T)
                comparison_matrix[i, j] = np.mean(corr_matrix[: data1.shape[1], data1.shape[1] :])
            elif comparison_metric == "euclidean":
                # Mean Euclidean distance
                distances = []
                for k in range(min(data1.shape[0], data2.shape[0])):
                    dist = np.linalg.norm(data1[k] - data2[k])
                    distances.append(dist)
                comparison_matrix[i, j] = np.mean(distances)
            elif comparison_metric == "cosine":
                # Mean cosine similarity
                similarities = []
                for k in range(min(data1.shape[0], data2.shape[0])):
                    cos_sim = np.dot(data1[k], data2[k]) / (np.linalg.norm(data1[k]) * np.linalg.norm(data2[k]))
                    similarities.append(cos_sim)
                comparison_matrix[i, j] = np.mean(similarities)

    # Plot comparison matrix
    im = ax.imshow(comparison_matrix, cmap="viridis", aspect="equal")
    ax.set_xticks(range(n_omics))
    ax.set_yticks(range(n_omics))
    ax.set_xticklabels(omics_types, rotation=45, ha="right")
    ax.set_yticklabels(omics_types)
    ax.set_title(f"Omics Layer Comparison ({comparison_metric})")

    # Add values on heatmap
    for i in range(n_omics):
        for j in range(n_omics):
            ax.text(
                j,
                i,
                ".2f",
                ha="center",
                va="center",
                color="white" if comparison_matrix[i, j] > np.mean(comparison_matrix) else "black",
            )

    plt.colorbar(im, ax=ax, label=comparison_metric.capitalize())

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Omics layer comparison saved to {output_path}")

    return ax


def plot_multiomics_pca(
    integrated_data: np.ndarray,
    omics_membership: np.ndarray | None = None,
    sample_labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot PCA of integrated multi-omics data.

    Args:
        integrated_data: Integrated feature matrix (samples x features)
        omics_membership: Optional array indicating which omics each feature belongs to
        sample_labels: Optional sample labels
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(integrated_data, np.ndarray, "integrated_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Perform PCA
    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(integrated_data)

    # Plot PCA
    if omics_membership is not None:
        unique_omics = np.unique(omics_membership)
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_omics)))

        for i, omics_type in enumerate(unique_omics):
            mask = omics_membership == omics_type
            ax.scatter(
                pca_coords[mask, 0], pca_coords[mask, 1], c=[colors[i]], label=f"Omics {omics_type}", alpha=0.7, s=50
            )
    else:
        scatter = ax.scatter(
            pca_coords[:, 0], pca_coords[:, 1], c=range(len(pca_coords)), cmap="viridis", alpha=0.7, s=50
        )

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title("Multi-Omics PCA")

    if omics_membership is not None:
        ax.legend()

    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics PCA saved to {output_path}")

    return ax


def plot_pathway_enrichment_integration(
    enrichment_results: Dict[str, List[Dict[str, Any]]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot integrated pathway enrichment across multiple omics.

    Args:
        enrichment_results: Dictionary mapping omics types to enrichment results
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(enrichment_results, dict, "enrichment_results")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract common pathways across omics
    all_pathways = set()
    for results in enrichment_results.values():
        all_pathways.update([r.get("pathway", r.get("term", "")) for r in results])

    all_pathways = list(all_pathways)
    n_pathways = len(all_pathways)
    omics_types = list(enrichment_results.keys())

    # Create enrichment matrix
    enrichment_matrix = np.full((n_pathways, len(omics_types)), np.nan)

    for i, pathway in enumerate(all_pathways):
        for j, omics_type in enumerate(omics_types):
            results = enrichment_results[omics_type]
            for result in results:
                if result.get("pathway", result.get("term", "")) == pathway:
                    enrichment_matrix[i, j] = -np.log10(result.get("pvalue", 1.0))
                    break

    # Plot heatmap
    if HAS_SEABORN:
        sns.heatmap(
            enrichment_matrix,
            cmap="Reds",
            ax=ax,
            xticklabels=omics_types,
            yticklabels=all_pathways,
            mask=np.isnan(enrichment_matrix),
            **kwargs,
        )
    else:
        im = ax.imshow(enrichment_matrix, cmap="Reds", aspect="auto")
        ax.set_xticks(range(len(omics_types)))
        ax.set_yticks(range(n_pathways))
        ax.set_xticklabels(omics_types, rotation=45, ha="right")
        ax.set_yticklabels(all_pathways)
        plt.colorbar(im, ax=ax)

    ax.set_title("Integrated Pathway Enrichment")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Pathway enrichment integration saved to {output_path}")

    return ax


def plot_multiomics_network(
    integration_network: Any,
    node_omics_types: Dict[str, str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 10),
    **kwargs,
) -> Axes:
    """Plot multi-omics integration network.

    Args:
        integration_network: NetworkX graph of integrated features
        node_omics_types: Mapping of node IDs to omics types
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required for multi-omics network visualization")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Get unique omics types and assign colors
    omics_types = list(set(node_omics_types.values()))
    colors = plt.cm.Set1(np.linspace(0, 1, len(omics_types)))
    omics_color_map = dict(zip(omics_types, colors))

    # Node colors based on omics type
    node_colors = [omics_color_map[node_omics_types.get(node, "unknown")] for node in integration_network.nodes()]

    # Layout
    pos = nx.spring_layout(integration_network, k=1, iterations=50)

    # Draw network
    nx.draw_networkx_nodes(integration_network, pos, node_color=node_colors, node_size=200, alpha=0.8, ax=ax)
    nx.draw_networkx_edges(integration_network, pos, edge_color="gray", alpha=0.3, ax=ax)

    # Add labels for important nodes only
    if len(integration_network.nodes()) <= 30:
        nx.draw_networkx_labels(integration_network, pos, font_size=8, ax=ax)

    ax.set_title("Multi-Omics Integration Network")
    ax.axis("off")

    # Add legend
    legend_elements = [
        plt.scatter([], [], c=color, label=omics_type, s=100) for omics_type, color in omics_color_map.items()
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics network saved to {output_path}")

    return ax


def plot_omics_factor_analysis(
    factor_loadings: np.ndarray,
    feature_names: List[str] | None = None,
    factor_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot factor analysis results for multi-omics data.

    Args:
        factor_loadings: Factor loading matrix (features x factors)
        feature_names: Optional names for features
        factor_names: Optional names for factors
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(factor_loadings, np.ndarray, "factor_loadings")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot factor loadings as heatmap
    if HAS_SEABORN:
        if factor_names:
            xticklabels = factor_names
        else:
            xticklabels = [f"Factor {i+1}" for i in range(factor_loadings.shape[1])]

        sns.heatmap(
            factor_loadings,
            cmap="RdBu_r",
            center=0,
            ax=ax,
            xticklabels=xticklabels,
            yticklabels=feature_names,
            **kwargs,
        )
    else:
        im = ax.imshow(factor_loadings, cmap="RdBu_r", aspect="auto")
        if factor_names:
            ax.set_xticks(range(len(factor_names)))
            ax.set_xticklabels(factor_names, rotation=45, ha="right")
        plt.colorbar(im, ax=ax)

    ax.set_title("Multi-Omics Factor Loadings")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Omics factor analysis saved to {output_path}")

    return ax


def plot_multiomics_upset(
    feature_sets: Dict[str, set],
    *,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Any:
    """Create an UpSet plot for multi-omics feature overlaps.

    Args:
        feature_sets: Dictionary mapping omics types to sets of features
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Figure object
    """
    try:
        import upsetplot
    except ImportError:
        raise ImportError("upsetplot required for UpSet plot visualization")

    # Create UpSet plot data
    from upsetplot import from_contents

    upset_data = from_contents(feature_sets)

    # Create plot
    fig = plt.figure(figsize=figsize)
    upsetplot.plot(upset_data, fig=fig, **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Multi-omics UpSet plot saved to {output_path}")

    return fig


def create_interactive_multiomics_dashboard(
    multiomics_data: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive multi-omics dashboard using Plotly.

    Args:
        multiomics_data: Dictionary containing multi-omics datasets and analysis results
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY or not make_subplots:
        raise ImportError("plotly required for interactive multi-omics dashboard")

    # Create subplot figure
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("Correlation Matrix", "PCA", "Integration Network", "Enrichment"),
        specs=[[{"type": "heatmap"}, {"type": "scatter"}], [{"type": "scatter"}, {"type": "bar"}]],
    )

    # Add correlation heatmap
    if "correlation_matrix" in multiomics_data:
        corr_matrix = multiomics_data["correlation_matrix"]
        fig.add_trace(go.Heatmap(z=corr_matrix, colorscale="RdBu", zmid=0), row=1, col=1)

    # Add PCA plot
    if "pca_coords" in multiomics_data:
        pca_coords = multiomics_data["pca_coords"]
        fig.add_trace(go.Scatter(x=pca_coords[:, 0], y=pca_coords[:, 1], mode="markers"), row=1, col=2)

    # Add network plot (simplified)
    if "network_data" in multiomics_data:
        network_data = multiomics_data["network_data"]
        # Simplified network visualization
        fig.add_trace(go.Scatter(x=[0, 1, 2], y=[0, 1, 0], mode="markers+lines"), row=2, col=1)

    fig.update_layout(title="Interactive Multi-Omics Dashboard", **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive multi-omics dashboard saved to {html_path}")

    return fig
