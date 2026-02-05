"""Population structure visualization for GWAS.

This module provides plots for visualizing population stratification and structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def _load_pca_from_file(pca_file: str | Path) -> Dict[str, Any]:
    """Load PCA data from a TSV file.

    Expects a TSV with columns: sample_id, PC1, PC2, ... PCn.

    Args:
        pca_file: Path to TSV file containing PCA data.

    Returns:
        Dictionary with 'pcs', 'explained_variance', 'sample_names'.
    """
    pca_file = Path(pca_file)
    lines = pca_file.read_text().strip().split("\n")
    header = lines[0].split("\t")
    pc_cols = [c for c in header if c.upper().startswith("PC")]
    pc_indices = [header.index(c) for c in pc_cols]

    sample_names = []
    pcs_data = []
    for line in lines[1:]:
        parts = line.split("\t")
        sample_names.append(parts[0])
        pcs_data.append([float(parts[i]) for i in pc_indices])

    pcs_array = np.array(pcs_data)
    # Compute approximate explained variance from variance of each PC column
    variances = np.var(pcs_array, axis=0)
    total_var = np.sum(variances)
    explained_variance = variances / total_var if total_var > 0 else np.zeros(len(pc_cols))

    return {
        "pcs": pcs_array,
        "explained_variance": explained_variance,
        "sample_names": sample_names,
    }


def pca_plot(
    pca_data: Dict[str, Any] | str | Path,
    output_file: Optional[str | Path] = None,
    title: str = "PCA Plot",
    pc1: int = 0,
    pc2: int = 1,
) -> Dict[str, Any]:
    """Create a PCA plot for population structure visualization.

    Args:
        pca_data: Dictionary containing PCA data with keys:
                 'pcs': 2D array of principal components
                 'explained_variance': array of explained variance ratios
                 'sample_names': list of sample names (optional)
                 'populations': list of population labels (optional)
                 OR a file path (str/Path) to a TSV file with PCA data.
        output_file: Optional output file path
        title: Plot title
        pc1: Index of first principal component to plot (0-based)
        pc2: Index of second principal component to plot (0-based)

    Returns:
        Dictionary with status and metadata.

    Example:
        >>> data = {'pcs': pcs_array, 'explained_variance': var_ratios}
        >>> result = pca_plot(data)
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib or seaborn not available for PCA plot")
        return {"status": "failed", "error": "matplotlib/seaborn not available"}

    try:
        # If pca_data is a file path, load it
        if isinstance(pca_data, (str, Path)):
            pca_data = _load_pca_from_file(pca_data)

        # Validate input data
        required_keys = ["pcs", "explained_variance"]
        missing_keys = [key for key in required_keys if key not in pca_data]
        if missing_keys:
            logger.error(f"Missing required keys in pca_data: {missing_keys}")
            return {"status": "failed", "error": f"Missing required keys: {missing_keys}"}

        pcs = pca_data["pcs"]
        explained_var = pca_data["explained_variance"]

        if len(pcs.shape) != 2:
            logger.error("PCA data 'pcs' must be a 2D array")
            return {"status": "failed", "error": "PCA data 'pcs' must be a 2D array"}

        if pc1 >= pcs.shape[1] or pc2 >= pcs.shape[1]:
            logger.error(f"PC indices {pc1}, {pc2} out of range for {pcs.shape[1]} components")
            return {"status": "failed", "error": f"PC indices out of range"}

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
                ax.annotate(
                    name, (x_data[i], y_data[i]), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8
                )

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

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": pcs.shape[0],
            "n_components": pcs.shape[1],
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating PCA plot: {e}")
        return {"status": "failed", "error": str(e)}


def pca_scree_plot(
    explained_variance: List[float] | str | Path,
    output_file: Optional[str | Path] = None,
    title: str = "PCA Scree Plot",
) -> Dict[str, Any]:
    """Create a scree plot showing explained variance by principal component.

    Args:
        explained_variance: List of explained variance ratios (0-1),
            OR a file path to a text file with one eigenvalue per line.
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dictionary with status and metadata.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for scree plot")
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # If explained_variance is a file path, load eigenvalues from it
        if isinstance(explained_variance, (str, Path)):
            eigenval_path = Path(explained_variance)
            lines = eigenval_path.read_text().strip().split("\n")
            explained_variance = [float(line.strip()) for line in lines if line.strip()]

        if not explained_variance:
            logger.error("No explained variance data provided")
            return {"status": "failed", "error": "No explained variance data provided"}

        # Convert to percentages if not already
        if max(explained_variance) <= 1.0:
            explained_variance_pct = [v * 100 for v in explained_variance]
        else:
            explained_variance_pct = list(explained_variance)

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
            ax.text(
                bar.get_x() + bar.get_width() / 2.0, height + 0.5, f"{pct:.1f}%", ha="center", va="bottom", fontsize=8
            )

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

        plt.close(fig)
        return {
            "status": "success",
            "n_components": len(explained_variance),
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating scree plot: {e}")
        return {"status": "failed", "error": str(e)}


def _load_kinship_from_file(kinship_file: str | Path) -> tuple[np.ndarray, List[str]]:
    """Load a kinship matrix from a TSV file.

    Expects a TSV with an optional header row of sample IDs, then numeric rows.

    Args:
        kinship_file: Path to TSV file containing kinship matrix.

    Returns:
        Tuple of (kinship_matrix as float ndarray, sample_labels list).
    """
    kinship_file = Path(kinship_file)
    lines = kinship_file.read_text().strip().split("\n")

    # Try to detect whether the first line is a header (non-numeric first field)
    first_fields = lines[0].split("\t")
    sample_labels: List[str] = []
    data_start = 0
    try:
        float(first_fields[0])
        # First field is numeric -- no header
    except ValueError:
        # First field is not numeric -- it is a header
        sample_labels = first_fields
        data_start = 1

    rows = []
    for line in lines[data_start:]:
        parts = line.split("\t")
        rows.append([float(v) for v in parts])

    matrix = np.array(rows, dtype=np.float64)
    return matrix, sample_labels


def kinship_heatmap(
    kinship_matrix: np.ndarray | List[List[float]] | str | Path,
    output_path_or_labels: Optional[str | Path | List[str]] = None,
    output_path_if_labels: Optional[str | Path] = None,
    title: Optional[str] = None,
    population_labels: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Create kinship matrix heatmap.

    This function supports multiple call signatures:
      - kinship_heatmap(file_path, output_path, ...)
      - kinship_heatmap(matrix, output_path, ...)
      - kinship_heatmap(matrix, sample_labels, output_path, ...)

    Args:
        kinship_matrix: Kinship matrix (ndarray, list of lists, or file path)
        output_path_or_labels: Path to save the plot, or sample labels list
            when output_path_if_labels provides the actual output path.
        output_path_if_labels: Output path when second arg is sample labels.
        title: Optional custom title for the heatmap.
        population_labels: Optional list of population labels per sample.

    Returns:
        Dictionary with status and metadata.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create kinship heatmap")
        return {"status": "failed", "error": "matplotlib not available"}

    logger.info("Creating kinship heatmap")

    try:
        sample_labels: Optional[List[str]] = None
        actual_output_path: Optional[Path] = None

        # Handle file path input (first arg is a file path string/Path)
        if isinstance(kinship_matrix, (str, Path)) and not isinstance(kinship_matrix, np.ndarray):
            mat_path = Path(kinship_matrix)
            if mat_path.is_file():
                kinship_matrix, sample_labels = _load_kinship_from_file(mat_path)
                if output_path_or_labels is not None and not isinstance(output_path_or_labels, list):
                    actual_output_path = Path(output_path_or_labels)

        # Determine output_path and sample_labels from positional args
        if isinstance(output_path_or_labels, list):
            # Called as: kinship_heatmap(matrix, sample_labels_list, output_path, ...)
            sample_labels = output_path_or_labels
            if output_path_if_labels is not None:
                actual_output_path = Path(output_path_if_labels)
        elif isinstance(output_path_or_labels, (str, Path)) and output_path_or_labels is not None:
            actual_output_path = Path(output_path_or_labels)

        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix, dtype=np.float64)

        # Ensure float dtype to avoid "Image data of dtype object cannot be converted to float"
        if isinstance(kinship_matrix, np.ndarray) and kinship_matrix.dtype == object:
            kinship_matrix = kinship_matrix.astype(np.float64)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))

        # Plot heatmap
        im = ax.imshow(kinship_matrix, cmap="viridis", aspect="equal")

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Kinship coefficient")

        # Labels and title
        plot_title = title if title else "Kinship Matrix Heatmap"
        ax.set_title(plot_title)
        ax.set_xlabel("Sample")
        ax.set_ylabel("Sample")

        # Add population color bar if population_labels provided
        if population_labels is not None and len(population_labels) == kinship_matrix.shape[0]:
            unique_pops = sorted(set(population_labels))
            pop_color_map = {pop: i for i, pop in enumerate(unique_pops)}
            pop_colors = [pop_color_map[p] for p in population_labels]

            # Add a thin color strip on top for populations
            pop_ax = fig.add_axes([0.125, 0.92, 0.62, 0.02])
            pop_arr = np.array(pop_colors).reshape(1, -1)
            pop_ax.imshow(pop_arr, aspect="auto", cmap="tab10")
            pop_ax.set_xticks([])
            pop_ax.set_yticks([])

        plt.tight_layout()

        # Save if output path provided
        if actual_output_path:
            fig.savefig(actual_output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Saved kinship heatmap to {actual_output_path}")

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": kinship_matrix.shape[0],
            "output_path": str(actual_output_path) if actual_output_path else None,
        }

    except Exception as e:
        logger.error(f"Error creating kinship heatmap: {e}")
        return {"status": "failed", "error": str(e)}


def _load_admixture_from_file(admixture_file: str | Path) -> np.ndarray:
    """Load admixture proportions from a whitespace/tab-delimited text file.

    Each line contains the ancestry proportions for one sample.

    Args:
        admixture_file: Path to admixture proportions file.

    Returns:
        2D numpy array of shape (n_samples, n_ancestries).
    """
    admixture_file = Path(admixture_file)
    lines = admixture_file.read_text().strip().split("\n")
    rows = []
    for line in lines:
        parts = line.strip().split()
        rows.append([float(v) for v in parts])
    return np.array(rows, dtype=np.float64)


def admixture_plot(
    admixture_data: Dict[str, Any] | str | Path,
    output_path: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 6),
    k: Optional[int] = None,
) -> Dict[str, Any]:
    """Create an admixture plot showing population ancestry proportions.

    Args:
        admixture_data: Dictionary containing admixture proportions and sample info,
            OR a file path (str/Path) to an ADMIXTURE-format text file.
        output_path: Path to save the plot
        figsize: Figure size
        k: Number of ancestral populations (informational, auto-detected from data)

    Returns:
        Dictionary with plot status and metadata
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if not HAS_MATPLOTLIB:
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # If admixture_data is a file path, load it
        if isinstance(admixture_data, (str, Path)):
            admixture_matrix = _load_admixture_from_file(admixture_data)
            sample_names: List[str] = []
            population_labels: List[str] = []
            ancestry_names: List[str] = []
        else:
            # Extract data from admixture_data dict
            admixture_proportions = admixture_data.get("admixture_proportions", [])
            sample_names = admixture_data.get("sample_names", [])
            population_labels = admixture_data.get("population_labels", [])
            ancestry_names = admixture_data.get("ancestry_names", [])

            if not admixture_proportions:
                return {"status": "failed", "error": "No admixture proportions provided"}

            admixture_matrix = np.array(admixture_proportions)

        if admixture_matrix.ndim != 2:
            return {"status": "failed", "error": "Admixture proportions must be 2D array"}

        n_samples, n_ancestries = admixture_matrix.shape

        # Sort samples by population if population labels are provided
        if population_labels and len(population_labels) == n_samples:
            # Group by population
            pop_indices: Dict[str, List[int]] = {}
            for i, pop in enumerate(population_labels):
                if pop not in pop_indices:
                    pop_indices[pop] = []
                pop_indices[pop].append(i)

            # Sort within populations by ancestry proportion
            sorted_indices: List[int] = []
            for pop in sorted(pop_indices.keys()):
                indices = pop_indices[pop]
                sorted_pop_indices = sorted(indices, key=lambda x: admixture_matrix[x, 0], reverse=True)
                sorted_indices.extend(sorted_pop_indices)
        else:
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
        plot_title = f"Population Admixture Proportions (K={k})" if k else "Population Admixture Proportions"
        ax.set_title(plot_title)
        ax.set_xlim(-0.5, n_samples - 0.5)
        ax.set_ylim(0, 1)

        # Create legend
        legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i], alpha=0.8) for i in range(n_ancestries)]
        ax.legend(legend_elements, ancestry_labels, loc="center left", bbox_to_anchor=(1.02, 0.5))

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": n_samples,
            "n_ancestries": n_ancestries,
            "k": k if k else n_ancestries,
            "populations": list(set(population_labels)) if population_labels else None,
            "output_path": str(output_path) if output_path else None,
        }

    except Exception as e:
        return {"status": "failed", "error": str(e)}
