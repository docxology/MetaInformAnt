"""Population structure visualization for GWAS.

This module provides plots for visualizing population stratification and structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def pca_plot(
    pca_data: Dict[str, Any],
    output_file: Optional[str | Path] = None,
    title: str = "PCA Plot",
    pc1: int = 0,
    pc2: int = 1,
) -> Optional[Any]:
    """Create a PCA plot for population structure visualization.

    Args:
        pca_data: Dictionary containing PCA data with keys:
                 'pcs': 2D array of principal components
                 'explained_variance': array of explained variance ratios
                 'sample_names': list of sample names (optional)
                 'populations': list of population labels (optional)
        output_file: Optional output file path
        title: Plot title
        pc1: Index of first principal component to plot (0-based)
        pc2: Index of second principal component to plot (0-based)

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> data = {'pcs': pcs_array, 'explained_variance': var_ratios}
        >>> plot = pca_plot(data)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for PCA plot")
        return None

    # Validate input data
    required_keys = ["pcs", "explained_variance"]
    missing_keys = [key for key in required_keys if key not in pca_data]
    if missing_keys:
        logger.error(f"Missing required keys in pca_data: {missing_keys}")
        return None

    pcs = pca_data["pcs"]
    explained_var = pca_data["explained_variance"]

    if len(pcs.shape) != 2:
        logger.error("PCA data 'pcs' must be a 2D array")
        return None

    if pc1 >= pcs.shape[1] or pc2 >= pcs.shape[1]:
        logger.error(f"PC indices {pc1}, {pc2} out of range for {pcs.shape[1]} components")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Extract PC data
    x_data = pcs[:, pc1]
    y_data = pcs[:, pc2]

    # Check for population labels
    populations = pca_data.get("populations")
    sample_names = pca_data.get("sample_names")

    if populations is not None and len(populations) == len(x_data):
        # Color by population
        unique_pops = list(set(populations))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_pops)))

        for i, pop in enumerate(unique_pops):
            mask = [p == pop for p in populations]
            ax.scatter(x_data[mask], y_data[mask], c=[colors[i]], label=pop, alpha=0.7, s=50)
        ax.legend(title="Population", bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        # Simple scatter plot
        ax.scatter(x_data, y_data, alpha=0.7, s=50, c="blue")

    # Add sample labels if provided and not too many samples
    if sample_names is not None and len(sample_names) <= 50:
        for i, name in enumerate(sample_names):
            ax.annotate(name, (x_data[i], y_data[i]), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8)

    # Labels and title
    var1 = explained_var[pc1] * 100
    var2 = explained_var[pc2] * 100

    ax.set_xlabel(f"PC{pc1+1} ({var1:.1f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC{pc2+1} ({var2:.1f}% variance)", fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.grid(True, alpha=0.3)

    # Add center lines
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved PCA plot to {output_file}")

    return plt.gcf()


def pca_scree_plot(
    explained_variance: List[float], output_file: Optional[str | Path] = None, title: str = "PCA Scree Plot"
) -> Optional[Any]:
    """Create a scree plot showing explained variance by principal component.

    Args:
        explained_variance: List of explained variance ratios (0-1)
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for scree plot")
        return None

    if not explained_variance:
        logger.error("No explained variance data provided")
        return None

    # Convert to percentages if not already
    if max(explained_variance) <= 1.0:
        explained_variance_pct = [v * 100 for v in explained_variance]
    else:
        explained_variance_pct = explained_variance

    fig, ax = plt.subplots(figsize=(10, 6))

    # Create scree plot
    components = list(range(1, len(explained_variance) + 1))
    bars = ax.bar(components, explained_variance_pct, color="skyblue", edgecolor="navy", alpha=0.7)

    # Add cumulative variance line
    cumulative = np.cumsum(explained_variance_pct)
    ax.plot(components, cumulative, "r-o", linewidth=2, markersize=4, label="Cumulative")

    # Add value labels on bars
    for bar, pct in zip(bars, explained_variance_pct):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, height + 0.5, f"{pct:.1f}%", ha="center", va="bottom", fontsize=8)

    ax.set_xlabel("Principal Component", fontsize=12)
    ax.set_ylabel("Explained Variance (%)", fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xticks(components)
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend()

    # Add elbow detection heuristic (find where additional variance drops below 5%)
    if len(explained_variance_pct) > 1:
        for i in range(1, len(explained_variance_pct)):
            if explained_variance_pct[i] < 5.0:
                ax.axvline(x=i + 0.5, color="red", linestyle="--", alpha=0.7, label=f"Elbow at PC{i+1}")
                break

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved PCA scree plot to {output_file}")

    return fig


def kinship_heatmap(kinship_matrix: np.ndarray | List[List[float]], output_path: Optional[str | Path] = None) -> Any:
    """Create kinship matrix heatmap.

    Args:
        kinship_matrix: Kinship matrix
        output_path: Path to save the plot (optional)

    Returns:
        matplotlib Figure object
    """
    try:
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create kinship heatmap")
        return None

    logger.info("Creating kinship heatmap")

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot heatmap
        im = ax.imshow(kinship_matrix, cmap="viridis", aspect="equal")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Kinship coefficient")

        # Labels and title
        ax.set_title("Kinship Matrix Heatmap")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Sample")

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved kinship heatmap to {output_path}")

        return fig

    except Exception as e:
        logger.error(f"Error creating kinship heatmap: {e}")
        return None


def admixture_plot(
    admixture_data: Dict[str, Any], output_path: Optional[str | Path] = None, figsize: tuple[int, int] = (12, 6)
) -> dict[str, Any]:
    """Create an admixture plot showing population ancestry proportions.

    Args:
        admixture_data: Dictionary containing admixture proportions and sample info
        output_path: Path to save the plot
        figsize: Figure size

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "error", "message": "matplotlib not available"}

    try:
        # Extract data from admixture_data
        admixture_proportions = admixture_data.get("admixture_proportions", [])
        sample_names = admixture_data.get("sample_names", [])
        population_labels = admixture_data.get("population_labels", [])
        ancestry_names = admixture_data.get("ancestry_names", [])

        if not admixture_proportions:
            return {"status": "error", "message": "No admixture proportions provided"}

        admixture_matrix = np.array(admixture_proportions)

        if admixture_matrix.ndim != 2:
            return {"status": "error", "message": "Admixture proportions must be 2D array"}

        n_samples, n_ancestries = admixture_matrix.shape

        # Sort samples by population if population labels are provided
        if population_labels and len(population_labels) == n_samples:
            # Group by population
            pop_indices = {}
            for i, pop in enumerate(population_labels):
                if pop not in pop_indices:
                    pop_indices[pop] = []
                pop_indices[pop].append(i)

            # Sort within populations by ancestry proportion
            sorted_indices = []
            for pop in sorted(pop_indices.keys()):
                indices = pop_indices[pop]
                # Sort by the proportion of the first ancestry component
                sorted_pop_indices = sorted(indices, key=lambda x: admixture_matrix[x, 0], reverse=True)
                sorted_indices.extend(sorted_pop_indices)
        else:
            # Sort by first ancestry component
            sorted_indices = list(range(n_samples))
            sorted_indices.sort(key=lambda x: admixture_matrix[x, 0], reverse=True)

        # Reorder data
        sorted_proportions = admixture_matrix[sorted_indices]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        # Create stacked bar plot
        bottom = np.zeros(n_samples)

        # Use a color map for ancestries
        colors = plt.cm.tab10(np.linspace(0, 1, n_ancestries))

        ancestry_labels = ancestry_names if ancestry_names else [f"Ancestry {i+1}" for i in range(n_ancestries)]

        for i in range(n_ancestries):
            ax.bar(
                range(n_samples),
                sorted_proportions[:, i],
                bottom=bottom,
                color=colors[i],
                width=1.0,
                edgecolor="none",
                alpha=0.8,
            )
            bottom += sorted_proportions[:, i]

        # Add population separators if population labels exist
        if population_labels and len(population_labels) == n_samples:
            sorted_pops = [population_labels[i] for i in sorted_indices]
            pop_changes = [i for i in range(1, len(sorted_pops)) if sorted_pops[i] != sorted_pops[i - 1]]

            for change_point in pop_changes:
                ax.axvline(x=change_point - 0.5, color="black", linewidth=1, alpha=0.7)

            # Add population labels
            unique_pops = sorted(set(population_labels))
            pop_positions = []
            for pop in unique_pops:
                indices = [i for i, p in enumerate(sorted_pops) if p == pop]
                if indices:
                    pop_positions.append((pop, sum(indices) / len(indices)))

            for pop, pos in pop_positions:
                ax.text(
                    pos,
                    -0.05,
                    pop,
                    ha="center",
                    va="top",
                    transform=ax.get_xaxis_transform(),
                    fontsize=10,
                    fontweight="bold",
                )

        ax.set_xlabel("Samples")
        ax.set_ylabel("Ancestry Proportion")
        ax.set_title("Population Admixture Proportions")
        ax.set_xlim(-0.5, n_samples - 0.5)
        ax.set_ylim(0, 1)

        # Create legend
        legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i], alpha=0.8) for i in range(n_ancestries)]
        ax.legend(legend_elements, ancestry_labels, loc="center left", bbox_to_anchor=(1.02, 0.5))

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        return {
            "status": "success",
            "n_samples": n_samples,
            "n_ancestries": n_ancestries,
            "populations": list(set(population_labels)) if population_labels else None,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "error", "message": str(e)}
