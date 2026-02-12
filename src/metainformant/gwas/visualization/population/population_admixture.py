"""Admixture and kinship population visualization for GWAS.

This module provides admixture plots, kinship heatmaps, dendrograms, and
clustermaps for visualizing population relatedness and ancestry proportions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


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


def kinship_dendrogram(
    kinship_matrix: np.ndarray | List[List[float]],
    sample_labels: Optional[List[str]] = None,
    output_file: Optional[str | Path] = None,
    method: str = "ward",
    color_by: Optional[str] = None,
    metadata: Optional[Dict[str, Dict[str, str]]] = None,
) -> Dict[str, Any]:
    """Create a hierarchical clustering dendrogram from a kinship matrix.

    Converts the kinship similarity matrix to a distance matrix (1 - kinship)
    and performs hierarchical clustering using scipy.

    Args:
        kinship_matrix: Square kinship matrix (similarity, not distance).
            Values typically range from 0 to 1 where 1 = self-kinship.
        sample_labels: Optional list of sample labels for leaf nodes.
        output_file: Optional path to save the figure.
        method: Linkage method for scipy (e.g., "ward", "complete", "average", "single").
        color_by: Metadata field name to color leaves by (e.g., "population").
        metadata: Dictionary of {sample_id: {field: value}} for coloring leaves.
            Used together with color_by to assign colors per population/group.

    Returns:
        Dictionary with keys:
            - status: "success", "failed", or "skipped"
            - output_path: Path to saved file or None
            - n_samples: Number of samples in the dendrogram
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    try:
        from scipy.cluster.hierarchy import dendrogram as scipy_dendrogram
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform

        HAS_SCIPY = True
    except ImportError:
        HAS_SCIPY = False

    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create kinship dendrogram")
        return {"status": "skipped", "output_path": None, "n_samples": 0}

    if not HAS_SCIPY:
        logger.warning("scipy not available, cannot create kinship dendrogram")
        return {"status": "skipped", "output_path": None, "n_samples": 0}

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix, dtype=np.float64)
        elif kinship_matrix.dtype != np.float64:
            kinship_matrix = kinship_matrix.astype(np.float64)

        n_samples = kinship_matrix.shape[0]

        # Edge case: fewer than 2 samples
        if n_samples < 2:
            logger.warning("Need at least 2 samples for a dendrogram")
            return {"status": "failed", "output_path": None, "n_samples": n_samples}

        # Generate default labels if not provided
        if sample_labels is None or len(sample_labels) != n_samples:
            sample_labels = [f"S{i}" for i in range(n_samples)]

        # Convert kinship (similarity) to distance: distance = 1 - kinship
        # Clip to avoid negative distances from numerical noise
        distance_matrix = np.clip(1.0 - kinship_matrix, 0.0, None)

        # Force symmetry and zero diagonal
        distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
        np.fill_diagonal(distance_matrix, 0.0)

        # Handle all-identical matrix (all zeros in distance) gracefully
        if np.allclose(distance_matrix, 0.0):
            logger.warning("All kinship values are identical; dendrogram will be flat")

        # Convert to condensed form for scipy linkage
        condensed_dist = squareform(distance_matrix, checks=False)

        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed_dist, method=method)

        # Determine leaf colors from metadata
        leaf_colors: Optional[Dict[str, str]] = None
        if color_by and metadata:
            # Build group-to-color mapping
            groups: List[str] = []
            for label in sample_labels:
                sample_meta = metadata.get(label, {})
                groups.append(sample_meta.get(color_by, "unknown"))

            unique_groups = sorted(set(groups))
            cmap = plt.cm.tab10(np.linspace(0, 1, max(len(unique_groups), 1)))
            group_to_color = {g: cmap[i] for i, g in enumerate(unique_groups)}

            # scipy dendrogram uses link_color_func; we colorize after via leaf labeling
            leaf_colors = {}
            for i, label in enumerate(sample_labels):
                grp = groups[i]
                rgba = group_to_color[grp]
                # Convert RGBA to hex for matplotlib
                hex_color = "#{:02x}{:02x}{:02x}".format(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))
                leaf_colors[label] = hex_color

        # Create figure
        fig_width = max(10, n_samples * 0.3)
        fig, ax = plt.subplots(figsize=(fig_width, 8))

        # Draw dendrogram
        dendro_result = scipy_dendrogram(
            linkage_matrix,
            labels=sample_labels,
            leaf_rotation=90,
            leaf_font_size=max(6, min(12, 200 // n_samples)),
            ax=ax,
        )

        # Color the leaf labels if metadata coloring is requested
        if leaf_colors:
            xlbls = ax.get_xticklabels()
            for lbl in xlbls:
                label_text = lbl.get_text()
                if label_text in leaf_colors:
                    lbl.set_color(leaf_colors[label_text])

            # Add a legend for population groups
            if color_by and metadata:
                from matplotlib.patches import Patch

                legend_elements = [Patch(facecolor=group_to_color[g], label=g) for g in unique_groups if g != "unknown"]
                if legend_elements:
                    ax.legend(
                        handles=legend_elements,
                        title=color_by.capitalize(),
                        bbox_to_anchor=(1.05, 1),
                        loc="upper left",
                    )

        ax.set_title("Kinship Dendrogram", fontsize=14, pad=20)
        ax.set_ylabel("Distance (1 - kinship)", fontsize=12)
        ax.set_xlabel("Samples", fontsize=12)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved kinship dendrogram to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_samples": n_samples,
        }

    except Exception as e:
        logger.error(f"Error creating kinship dendrogram: {e}")
        return {"status": "failed", "output_path": None, "n_samples": 0}


def kinship_clustermap(
    kinship_matrix: np.ndarray | List[List[float]],
    sample_labels: Optional[List[str]] = None,
    output_file: Optional[str | Path] = None,
    method: str = "ward",
    metadata: Optional[Dict[str, Dict[str, str]]] = None,
    annotate_populations: bool = True,
) -> Dict[str, Any]:
    """Create a combined heatmap with dendrograms on both axes from a kinship matrix.

    Produces a clustermap-style visualization with hierarchical clustering
    dendrograms along both axes and optional population color annotation bars.

    Args:
        kinship_matrix: Square kinship matrix (similarity values).
        sample_labels: Optional list of sample labels.
        output_file: Optional path to save the figure.
        method: Linkage method for scipy clustering.
        metadata: Dictionary of {sample_id: {population: str, ...}} for annotations.
        annotate_populations: Whether to add population color bars when metadata
            is available. Defaults to True.

    Returns:
        Dictionary with keys:
            - status: "success", "failed", or "skipped"
            - output_path: Path to saved file or None
            - n_samples: Number of samples in the clustermap
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.gridspec as gridspec
        import matplotlib.pyplot as plt

        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    try:
        from scipy.cluster.hierarchy import dendrogram as scipy_dendrogram
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform

        HAS_SCIPY = True
    except ImportError:
        HAS_SCIPY = False

    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available, cannot create kinship clustermap")
        return {"status": "skipped", "output_path": None, "n_samples": 0}

    if not HAS_SCIPY:
        logger.warning("scipy not available, cannot create kinship clustermap")
        return {"status": "skipped", "output_path": None, "n_samples": 0}

    try:
        # Convert to numpy array if needed
        if isinstance(kinship_matrix, list):
            kinship_matrix = np.array(kinship_matrix, dtype=np.float64)
        elif kinship_matrix.dtype != np.float64:
            kinship_matrix = kinship_matrix.astype(np.float64)

        n_samples = kinship_matrix.shape[0]

        # Edge case: fewer than 2 samples
        if n_samples < 2:
            logger.warning("Need at least 2 samples for a clustermap")
            return {"status": "failed", "output_path": None, "n_samples": n_samples}

        # Generate default labels if not provided
        if sample_labels is None or len(sample_labels) != n_samples:
            sample_labels = [f"S{i}" for i in range(n_samples)]

        # Convert kinship (similarity) to distance
        distance_matrix = np.clip(1.0 - kinship_matrix, 0.0, None)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2.0
        np.fill_diagonal(distance_matrix, 0.0)

        # Condensed distance for scipy
        condensed_dist = squareform(distance_matrix, checks=False)

        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed_dist, method=method)

        # Determine population annotations
        has_pop_annotations = False
        pop_groups: List[str] = []
        unique_pops: List[str] = []
        pop_color_map: Dict[str, Any] = {}
        if annotate_populations and metadata:
            for label in sample_labels:
                sample_meta = metadata.get(label, {})
                pop_groups.append(sample_meta.get("population", "unknown"))
            unique_pops = sorted(set(pop_groups))
            if len(unique_pops) > 1 or (len(unique_pops) == 1 and unique_pops[0] != "unknown"):
                has_pop_annotations = True
                cmap_colors = plt.cm.tab10(np.linspace(0, 1, max(len(unique_pops), 1)))
                pop_color_map = {pop: cmap_colors[i] for i, pop in enumerate(unique_pops)}

        # Build the layout using GridSpec
        # Layout: top dendrogram, optional color bar, heatmap + left dendrogram
        fig_size = max(10, n_samples * 0.25)
        fig = plt.figure(figsize=(fig_size + 2, fig_size))

        # Define grid ratios
        if has_pop_annotations:
            # Rows: top-dendrogram, color-bar, heatmap
            # Cols: left-dendrogram, color-bar, heatmap
            gs = gridspec.GridSpec(
                3,
                3,
                width_ratios=[0.15, 0.02, 0.83],
                height_ratios=[0.15, 0.02, 0.83],
                wspace=0.01,
                hspace=0.01,
            )
        else:
            # Rows: top-dendrogram, heatmap
            # Cols: left-dendrogram, heatmap
            gs = gridspec.GridSpec(
                2,
                2,
                width_ratios=[0.15, 0.85],
                height_ratios=[0.15, 0.85],
                wspace=0.01,
                hspace=0.01,
            )

        # --- Top dendrogram ---
        if has_pop_annotations:
            ax_top_dendro = fig.add_subplot(gs[0, 2])
        else:
            ax_top_dendro = fig.add_subplot(gs[0, 1])

        dendro_top = scipy_dendrogram(
            linkage_matrix,
            orientation="top",
            no_labels=True,
            ax=ax_top_dendro,
        )
        ax_top_dendro.set_xticks([])
        ax_top_dendro.set_yticks([])
        ax_top_dendro.spines["top"].set_visible(False)
        ax_top_dendro.spines["right"].set_visible(False)
        ax_top_dendro.spines["bottom"].set_visible(False)
        ax_top_dendro.spines["left"].set_visible(False)

        # Get the reordered indices from the dendrogram
        reorder_idx = dendro_top["leaves"]

        # --- Left dendrogram ---
        if has_pop_annotations:
            ax_left_dendro = fig.add_subplot(gs[2, 0])
        else:
            ax_left_dendro = fig.add_subplot(gs[1, 0])

        scipy_dendrogram(
            linkage_matrix,
            orientation="left",
            no_labels=True,
            ax=ax_left_dendro,
        )
        ax_left_dendro.set_xticks([])
        ax_left_dendro.set_yticks([])
        ax_left_dendro.spines["top"].set_visible(False)
        ax_left_dendro.spines["right"].set_visible(False)
        ax_left_dendro.spines["bottom"].set_visible(False)
        ax_left_dendro.spines["left"].set_visible(False)
        ax_left_dendro.invert_yaxis()

        # --- Population color bars ---
        if has_pop_annotations:
            # Top color bar (columns)
            ax_col_colors = fig.add_subplot(gs[1, 2])
            col_colors_arr = np.array([pop_color_map[pop_groups[i]] for i in reorder_idx]).reshape(1, -1, 4)
            ax_col_colors.imshow(col_colors_arr, aspect="auto", interpolation="nearest")
            ax_col_colors.set_xticks([])
            ax_col_colors.set_yticks([])

            # Left color bar (rows)
            ax_row_colors = fig.add_subplot(gs[2, 1])
            row_colors_arr = np.array([pop_color_map[pop_groups[i]] for i in reorder_idx]).reshape(-1, 1, 4)
            ax_row_colors.imshow(row_colors_arr, aspect="auto", interpolation="nearest")
            ax_row_colors.set_xticks([])
            ax_row_colors.set_yticks([])

        # --- Main heatmap ---
        if has_pop_annotations:
            ax_heatmap = fig.add_subplot(gs[2, 2])
        else:
            ax_heatmap = fig.add_subplot(gs[1, 1])

        # Reorder kinship matrix according to dendrogram
        reordered_matrix = kinship_matrix[np.ix_(reorder_idx, reorder_idx)]

        im = ax_heatmap.imshow(reordered_matrix, cmap="viridis", aspect="equal", interpolation="nearest")

        # Add sample labels if not too many
        reordered_labels = [sample_labels[i] for i in reorder_idx]
        if n_samples <= 50:
            font_size = max(5, min(10, 150 // n_samples))
            ax_heatmap.set_xticks(range(n_samples))
            ax_heatmap.set_xticklabels(reordered_labels, rotation=90, fontsize=font_size)
            ax_heatmap.set_yticks(range(n_samples))
            ax_heatmap.set_yticklabels(reordered_labels, fontsize=font_size)
        else:
            ax_heatmap.set_xticks([])
            ax_heatmap.set_yticks([])

        # Colorbar
        cbar = fig.colorbar(im, ax=ax_heatmap, fraction=0.046, pad=0.04)
        cbar.set_label("Kinship coefficient", fontsize=10)

        # Add population legend if applicable
        if has_pop_annotations:
            from matplotlib.patches import Patch

            legend_elements = [Patch(facecolor=pop_color_map[p], label=p) for p in unique_pops]
            fig.legend(
                handles=legend_elements,
                title="Population",
                loc="upper right",
                bbox_to_anchor=(0.98, 0.98),
                fontsize=8,
            )

        fig.suptitle("Kinship Clustermap", fontsize=14, y=1.02)

        output_path_str: Optional[str] = None
        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            output_path_str = str(output_file)
            logger.info(f"Saved kinship clustermap to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_samples": n_samples,
        }

    except Exception as e:
        logger.error(f"Error creating kinship clustermap: {e}")
        return {"status": "failed", "output_path": None, "n_samples": 0}
