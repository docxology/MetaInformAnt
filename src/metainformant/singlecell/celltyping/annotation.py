"""Cell type annotation and classification using marker gene sets.

Provides marker-based cell type annotation using multiple methods (overlap
scoring, correlation, average expression), label transfer from annotated
reference datasets via kNN in PCA space, novel cell type detection for cells
that do not confidently match any known type, and cell type composition
analysis across samples or experimental groups.

Methods:
    - Overlap scoring: Enrichment based on intersection of expressed genes
      with marker sets, normalized by set size and background.
    - Correlation: Pearson correlation of cell expression profiles against
      mean reference profiles for each cell type.
    - Scoring: Mean expression of marker genes minus mean expression of a
      randomly selected background set.
    - Label transfer: kNN classification in a shared PCA embedding space
      between reference and query datasets.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


def annotate_by_markers(
    expression_matrix: Any,
    marker_genes: dict[str, list[str]],
    gene_names: list[str] | None = None,
    method: str = "overlap",
    threshold: float = 0.0,
) -> dict:
    """Annotate cell types using known marker gene sets.

    Assigns each cell (row) in the expression matrix a cell type label based
    on how strongly its expression profile matches the provided marker gene
    sets. Three scoring methods are available: overlap enrichment, correlation
    with mean marker expression, and average marker expression scoring.

    Args:
        expression_matrix: Expression matrix (cells x genes). Accepts lists
            of lists or numpy arrays.
        marker_genes: Dictionary mapping cell type names to lists of marker
            gene names. Example: {"T_cell": ["CD3D", "CD3E"], "B_cell": ["CD19"]}.
        gene_names: List of gene names corresponding to columns in the
            expression matrix. Required if expression_matrix is not annotated.
        method: Scoring method. One of "overlap" (enrichment score based on
            expressed gene overlap), "correlation" (Pearson correlation with
            mean marker profile), or "scoring" (mean marker expression minus
            background).
        threshold: Minimum score required for a confident annotation. Cells
            scoring below this for all types are labeled "unassigned".

    Returns:
        Dictionary with keys:
            - cell_labels: list[str] of assigned cell type for each cell.
            - confidence_scores: list[float] of best score per cell.
            - ambiguous_cells: list[int] of indices where top two scores
              are within 10% of each other.
            - marker_stats: dict mapping cell type to per-type statistics
              (mean_score, n_cells_assigned, marker_coverage).

    Raises:
        ValueError: If method is not recognized, or dimensions are inconsistent.
    """
    valid_methods = ("overlap", "correlation", "scoring")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    # Convert to list-of-lists for uniform handling
    matrix = _to_list_matrix(expression_matrix)
    n_cells = len(matrix)
    if n_cells == 0:
        raise ValueError("Expression matrix must contain at least one cell (row)")

    n_genes = len(matrix[0])
    if gene_names is not None and len(gene_names) != n_genes:
        raise ValueError(
            f"gene_names length ({len(gene_names)}) does not match " f"expression matrix columns ({n_genes})"
        )

    logger.info(f"Annotating {n_cells} cells with {len(marker_genes)} cell types " f"using method='{method}'")

    # Build gene name to index mapping
    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]
    gene_to_idx: dict[str, int] = {g: i for i, g in enumerate(gene_names)}

    # Resolve marker indices for each cell type
    marker_indices: dict[str, list[int]] = {}
    for ct, markers in marker_genes.items():
        indices = [gene_to_idx[g] for g in markers if g in gene_to_idx]
        if indices:
            marker_indices[ct] = indices
        else:
            logger.warning(f"No marker genes found in expression matrix for cell type '{ct}'")

    if not marker_indices:
        raise ValueError("None of the provided marker genes were found in the expression matrix")

    # Score every cell against every cell type
    cell_type_names = list(marker_indices.keys())
    scores_per_cell: list[list[float]] = []

    for cell_idx in range(n_cells):
        row = matrix[cell_idx]
        cell_scores: list[float] = []
        for ct in cell_type_names:
            if method == "overlap":
                s = _overlap_score(row, marker_indices[ct], n_genes)
            elif method == "correlation":
                s = _correlation_score(row, marker_indices[ct])
            else:  # scoring
                s = _expression_score(row, marker_indices[ct])
            cell_scores.append(s)
        scores_per_cell.append(cell_scores)

    # Assign labels based on best score
    cell_labels: list[str] = []
    confidence_scores: list[float] = []
    ambiguous_cells: list[int] = []

    for cell_idx, cell_scores in enumerate(scores_per_cell):
        best_idx = _argmax(cell_scores)
        best_score = cell_scores[best_idx]

        if best_score < threshold:
            cell_labels.append("unassigned")
            confidence_scores.append(best_score)
        else:
            cell_labels.append(cell_type_names[best_idx])
            confidence_scores.append(best_score)

        # Check for ambiguity: top two scores within 10%
        if len(cell_scores) >= 2:
            sorted_scores = sorted(cell_scores, reverse=True)
            if sorted_scores[0] > 0 and sorted_scores[1] / sorted_scores[0] > 0.9:
                ambiguous_cells.append(cell_idx)

    # Compute per-type statistics
    marker_stats: dict[str, dict[str, Any]] = {}
    label_counts = Counter(cell_labels)
    for ct in cell_type_names:
        ct_idx = cell_type_names.index(ct)
        ct_scores = [scores_per_cell[i][ct_idx] for i in range(n_cells)]
        mean_score = sum(ct_scores) / len(ct_scores) if ct_scores else 0.0
        n_markers_found = len(marker_indices[ct])
        n_markers_total = len(marker_genes[ct])
        marker_stats[ct] = {
            "mean_score": mean_score,
            "n_cells_assigned": label_counts.get(ct, 0),
            "marker_coverage": n_markers_found / n_markers_total if n_markers_total > 0 else 0.0,
            "n_markers_found": n_markers_found,
            "n_markers_total": n_markers_total,
        }

    logger.info(
        f"Annotation complete: {n_cells - label_counts.get('unassigned', 0)}/{n_cells} "
        f"cells assigned, {len(ambiguous_cells)} ambiguous"
    )

    return {
        "cell_labels": cell_labels,
        "confidence_scores": confidence_scores,
        "ambiguous_cells": ambiguous_cells,
        "marker_stats": marker_stats,
    }


def score_cell_type(
    expression: list[float],
    gene_names: list[str],
    marker_set: list[str],
    n_background: int = 50,
    seed: int | None = None,
) -> float:
    """Compute cell type score for a single cell using marker gene expression.

    Calculates the mean expression of marker genes minus the mean expression
    of a randomly selected background gene set of the same size. This
    corrects for overall expression level and highlights marker-specific
    enrichment.

    Args:
        expression: Expression values for a single cell (one value per gene).
        gene_names: Gene names corresponding to expression values.
        marker_set: List of marker gene names for the cell type of interest.
        n_background: Number of background genes to sample for normalization.
        seed: Random seed for reproducible background selection.

    Returns:
        Float score representing marker enrichment. Positive values indicate
        the cell expresses the marker set above background levels.

    Raises:
        ValueError: If expression and gene_names have different lengths.
    """
    if len(expression) != len(gene_names):
        raise ValueError(f"expression length ({len(expression)}) must match " f"gene_names length ({len(gene_names)})")

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Get marker gene expression values
    marker_values: list[float] = []
    for gene in marker_set:
        if gene in gene_to_idx:
            marker_values.append(expression[gene_to_idx[gene]])

    if not marker_values:
        return 0.0

    marker_mean = sum(marker_values) / len(marker_values)

    # Sample background genes (excluding markers)
    marker_set_lower = {g.lower() for g in marker_set}
    non_marker_indices = [i for i, g in enumerate(gene_names) if g.lower() not in marker_set_lower]

    if not non_marker_indices:
        return marker_mean

    rng = random.Random(seed)
    n_bg = min(n_background, len(non_marker_indices))
    bg_indices = rng.sample(non_marker_indices, n_bg)
    bg_values = [expression[i] for i in bg_indices]
    bg_mean = sum(bg_values) / len(bg_values) if bg_values else 0.0

    return marker_mean - bg_mean


def transfer_labels(
    reference_data: Any,
    reference_labels: list[str],
    query_data: Any,
    n_neighbors: int = 30,
    n_components: int = 50,
) -> dict:
    """Transfer cell type labels from a reference dataset to a query dataset.

    Projects both reference and query data into a shared PCA space, then
    uses a k-nearest-neighbor classifier to assign reference labels to
    query cells. Prediction confidence is derived from the fraction of
    neighbors agreeing on the assigned label.

    Args:
        reference_data: Reference expression matrix (cells x genes).
        reference_labels: Cell type labels for each reference cell.
        query_data: Query expression matrix (cells x genes). Must have
            the same number of genes (columns) as reference_data.
        n_neighbors: Number of nearest neighbors for kNN classification.
        n_components: Number of PCA components for dimensionality reduction
            before neighbor search.

    Returns:
        Dictionary with keys:
            - predicted_labels: list[str] of predicted cell type per query cell.
            - prediction_scores: list[float] of confidence (fraction of
              neighbors with the predicted label) per query cell.
            - mapping_quality: dict with overall transfer statistics
              (mean_confidence, n_high_confidence, n_low_confidence).

    Raises:
        ValueError: If reference and query have different gene counts, or
            reference_labels length does not match reference rows.
    """
    ref_matrix = _to_list_matrix(reference_data)
    query_matrix = _to_list_matrix(query_data)

    n_ref = len(ref_matrix)
    n_query = len(query_matrix)

    if n_ref == 0 or n_query == 0:
        raise ValueError("Both reference and query must contain at least one cell")

    n_genes_ref = len(ref_matrix[0])
    n_genes_query = len(query_matrix[0])
    if n_genes_ref != n_genes_query:
        raise ValueError(
            f"Reference ({n_genes_ref} genes) and query ({n_genes_query} genes) " f"must have the same number of genes"
        )

    if len(reference_labels) != n_ref:
        raise ValueError(f"reference_labels length ({len(reference_labels)}) must match " f"reference rows ({n_ref})")

    logger.info(
        f"Transferring labels from {n_ref} reference cells to {n_query} query cells "
        f"(k={n_neighbors}, n_components={n_components})"
    )

    # Simple PCA via covariance eigendecomposition (pure Python fallback)
    combined = ref_matrix + query_matrix
    n_total = n_ref + n_query
    n_genes = n_genes_ref
    actual_components = min(n_components, n_genes, n_total)

    # Center the data
    means = [0.0] * n_genes
    for row in combined:
        for j in range(n_genes):
            means[j] += row[j]
    means = [m / n_total for m in means]

    centered = [[row[j] - means[j] for j in range(n_genes)] for row in combined]

    # Use numpy PCA if available, otherwise use simplified projection
    if HAS_NUMPY:
        X = np.array(centered, dtype=np.float64)
        # Truncated SVD for PCA
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        projected = (U[:, :actual_components] * S[:actual_components]).tolist()
    else:
        # Fallback: use first actual_components features as pseudo-projection
        projected = [row[:actual_components] for row in centered]

    ref_projected = projected[:n_ref]
    query_projected = projected[n_ref:]

    # kNN classification
    effective_k = min(n_neighbors, n_ref)
    predicted_labels: list[str] = []
    prediction_scores: list[float] = []

    for qi in range(n_query):
        # Compute distances to all reference cells
        dists: list[tuple[float, int]] = []
        for ri in range(n_ref):
            d = _euclidean_distance(query_projected[qi], ref_projected[ri])
            dists.append((d, ri))
        dists.sort(key=lambda x: x[0])

        # Count labels among k nearest neighbors
        neighbor_labels = [reference_labels[dists[i][1]] for i in range(effective_k)]
        label_counts = Counter(neighbor_labels)
        best_label, best_count = label_counts.most_common(1)[0]

        predicted_labels.append(best_label)
        prediction_scores.append(best_count / effective_k)

    # Mapping quality summary
    high_conf = sum(1 for s in prediction_scores if s >= 0.7)
    low_conf = sum(1 for s in prediction_scores if s < 0.3)
    mean_conf = sum(prediction_scores) / len(prediction_scores) if prediction_scores else 0.0

    mapping_quality = {
        "mean_confidence": mean_conf,
        "n_high_confidence": high_conf,
        "n_low_confidence": low_conf,
        "n_query_cells": n_query,
        "n_reference_cells": n_ref,
        "n_neighbors": effective_k,
        "n_components": actual_components,
    }

    logger.info(
        f"Label transfer complete: mean confidence={mean_conf:.3f}, "
        f"{high_conf} high-confidence, {low_conf} low-confidence"
    )

    return {
        "predicted_labels": predicted_labels,
        "prediction_scores": prediction_scores,
        "mapping_quality": mapping_quality,
    }


def find_novel_types(
    expression_matrix: Any,
    known_labels: list[str],
    marker_genes: dict[str, list[str]] | None = None,
    gene_names: list[str] | None = None,
    threshold: float = 0.5,
) -> dict:
    """Identify cells that do not match any known cell type.

    Examines the confidence scores from marker-based annotation (or the
    provided known_labels with their distribution) and flags cells whose
    best cell type score falls below the given threshold, suggesting they
    may represent novel or transitional cell types.

    Args:
        expression_matrix: Expression matrix (cells x genes).
        known_labels: Current cell type labels for each cell (e.g., from
            annotate_by_markers). Cells labeled "unassigned" or with low
            confidence are candidates.
        marker_genes: Marker gene sets used for annotation. If provided,
            cells are re-scored to identify low-confidence assignments.
        gene_names: Gene names for columns in expression_matrix.
        threshold: Confidence threshold below which cells are considered
            potentially novel. Default 0.5.

    Returns:
        Dictionary with keys:
            - novel_cell_indices: list[int] of cell indices flagged as
              potentially novel.
            - novel_cell_scores: list[float] of best-match scores for
              novel cells (all below threshold).
            - n_novel: int total count of novel cells.
            - fraction_novel: float fraction of all cells that are novel.
            - cluster_summary: dict mapping provisional cluster IDs to
              lists of novel cell indices (simple distance-based grouping).
    """
    matrix = _to_list_matrix(expression_matrix)
    n_cells = len(matrix)

    if len(known_labels) != n_cells:
        raise ValueError(f"known_labels length ({len(known_labels)}) must match " f"expression matrix rows ({n_cells})")

    logger.info(f"Searching for novel cell types among {n_cells} cells (threshold={threshold})")

    novel_indices: list[int] = []
    novel_scores: list[float] = []

    if marker_genes is not None and gene_names is not None:
        # Re-score cells to get confidence values
        result = annotate_by_markers(expression_matrix, marker_genes, gene_names=gene_names, method="scoring")
        for i, score in enumerate(result["confidence_scores"]):
            if score < threshold:
                novel_indices.append(i)
                novel_scores.append(score)
    else:
        # Use known_labels: cells marked "unassigned" or "unknown"
        unassigned_tags = {"unassigned", "unknown", "na", "none", ""}
        for i, label in enumerate(known_labels):
            if label.lower().strip() in unassigned_tags:
                novel_indices.append(i)
                novel_scores.append(0.0)

    # Simple distance-based grouping of novel cells
    cluster_summary: dict[int, list[int]] = {}
    if novel_indices and len(novel_indices) > 1:
        novel_rows = [matrix[i] for i in novel_indices]
        clusters = _simple_cluster(novel_rows, max_clusters=min(10, len(novel_indices)))
        for idx, cluster_id in enumerate(clusters):
            if cluster_id not in cluster_summary:
                cluster_summary[cluster_id] = []
            cluster_summary[cluster_id].append(novel_indices[idx])
    elif novel_indices:
        cluster_summary[0] = novel_indices

    fraction = len(novel_indices) / n_cells if n_cells > 0 else 0.0

    logger.info(f"Found {len(novel_indices)} potentially novel cells " f"({fraction:.1%} of total)")

    return {
        "novel_cell_indices": novel_indices,
        "novel_cell_scores": novel_scores,
        "n_novel": len(novel_indices),
        "fraction_novel": fraction,
        "cluster_summary": cluster_summary,
    }


def cell_type_composition(
    labels: list[str],
    groups: list[str] | None = None,
    confidence_level: float = 0.95,
) -> dict:
    """Compute cell type proportions per sample or group.

    Calculates the fraction of each cell type within each group (or
    overall if no groups are provided), along with confidence intervals
    computed using the Wilson score interval for binomial proportions.

    Args:
        labels: Cell type label for each cell.
        groups: Optional group/sample label for each cell. If provided,
            proportions are computed separately per group. Must have the
            same length as labels.
        confidence_level: Confidence level for proportion intervals
            (default 0.95).

    Returns:
        Dictionary with keys:
            - overall: dict mapping cell type to proportion across all cells.
            - per_group: dict mapping group name to dict of cell type
              proportions (only present if groups is provided).
            - confidence_intervals: dict mapping (group, cell_type) tuples
              as "group::cell_type" strings to (lower, upper) tuples.
            - n_cells: int total cell count.
            - n_types: int number of unique cell types.
            - n_groups: int number of unique groups (1 if no groups).

    Raises:
        ValueError: If labels is empty or groups length mismatches.
    """
    if not labels:
        raise ValueError("labels must not be empty")

    if groups is not None and len(groups) != len(labels):
        raise ValueError(f"groups length ({len(groups)}) must match labels length ({len(labels)})")

    n_cells = len(labels)
    all_types = sorted(set(labels))
    overall_counts = Counter(labels)
    overall_proportions = {ct: count / n_cells for ct, count in overall_counts.items()}

    logger.info(f"Computing composition for {n_cells} cells, " f"{len(all_types)} cell types")

    result: dict[str, Any] = {
        "overall": overall_proportions,
        "n_cells": n_cells,
        "n_types": len(all_types),
    }

    # Z-score for confidence interval
    z = _z_score_for_confidence(confidence_level)

    confidence_intervals: dict[str, tuple[float, float]] = {}

    if groups is not None:
        unique_groups = sorted(set(groups))
        per_group: dict[str, dict[str, float]] = {}

        for grp in unique_groups:
            grp_indices = [i for i, g in enumerate(groups) if g == grp]
            grp_labels = [labels[i] for i in grp_indices]
            grp_n = len(grp_labels)
            grp_counts = Counter(grp_labels)

            grp_props: dict[str, float] = {}
            for ct in all_types:
                count = grp_counts.get(ct, 0)
                prop = count / grp_n if grp_n > 0 else 0.0
                grp_props[ct] = prop

                # Wilson score interval
                lower, upper = _wilson_interval(count, grp_n, z)
                key = f"{grp}::{ct}"
                confidence_intervals[key] = (lower, upper)

            per_group[grp] = grp_props

        result["per_group"] = per_group
        result["n_groups"] = len(unique_groups)
    else:
        result["n_groups"] = 1
        for ct in all_types:
            count = overall_counts[ct]
            lower, upper = _wilson_interval(count, n_cells, z)
            key = f"all::{ct}"
            confidence_intervals[key] = (lower, upper)

    result["confidence_intervals"] = confidence_intervals

    logger.info(f"Composition analysis complete for {result['n_groups']} group(s)")
    return result


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _to_list_matrix(data: Any) -> list[list[float]]:
    """Convert various matrix types to list of lists of floats.

    Args:
        data: Expression matrix. Supports list[list[float]], numpy arrays,
            scipy sparse matrices, or objects with a .toarray() method.

    Returns:
        List of lists of floats.
    """
    if isinstance(data, list):
        return data

    # numpy array
    if HAS_NUMPY and isinstance(data, np.ndarray):
        return data.tolist()

    # scipy sparse or anndata-like
    if hasattr(data, "toarray"):
        return data.toarray().tolist()

    # pandas DataFrame
    if hasattr(data, "values"):
        return data.values.tolist()

    raise TypeError(f"Unsupported expression matrix type: {type(data)}")


def _overlap_score(
    row: list[float],
    marker_indices: list[int],
    n_genes: int,
) -> float:
    """Compute overlap enrichment score.

    Measures the fraction of marker genes expressed (value > 0) in the
    cell, normalized by the overall expression rate.

    Args:
        row: Expression values for a single cell.
        marker_indices: Column indices of marker genes.
        n_genes: Total number of genes.

    Returns:
        Enrichment score (higher = stronger marker expression).
    """
    n_expressed_total = sum(1 for v in row if v > 0)
    if n_expressed_total == 0:
        return 0.0

    n_markers_expressed = sum(1 for idx in marker_indices if row[idx] > 0)
    n_markers = len(marker_indices)

    if n_markers == 0:
        return 0.0

    # Enrichment: observed / expected
    marker_rate = n_markers_expressed / n_markers
    background_rate = n_expressed_total / n_genes if n_genes > 0 else 0.0

    if background_rate == 0:
        return marker_rate

    return marker_rate / background_rate


def _correlation_score(
    row: list[float],
    marker_indices: list[int],
) -> float:
    """Compute Pearson correlation between cell expression and marker profile.

    Creates a binary marker profile (1 at marker positions, 0 elsewhere)
    and correlates it with the cell's expression values at those positions.

    Args:
        row: Expression values for a single cell.
        marker_indices: Column indices of marker genes.

    Returns:
        Pearson correlation coefficient (range -1 to 1).
    """
    if not marker_indices:
        return 0.0

    marker_values = [row[i] for i in marker_indices]
    n = len(marker_values)

    if n < 2:
        return marker_values[0] if marker_values else 0.0

    mean_val = sum(marker_values) / n
    # Use all-ones as the "ideal" marker profile, so correlation is
    # simply the normalized deviation from the mean
    mean_ideal = 1.0
    ideal = [1.0] * n

    numerator = sum((marker_values[i] - mean_val) * (ideal[i] - mean_ideal) for i in range(n))
    denom_a = math.sqrt(sum((v - mean_val) ** 2 for v in marker_values))
    denom_b = math.sqrt(sum((v - mean_ideal) ** 2 for v in ideal))

    if denom_a == 0 or denom_b == 0:
        return 0.0

    return numerator / (denom_a * denom_b)


def _expression_score(
    row: list[float],
    marker_indices: list[int],
) -> float:
    """Compute mean expression score for marker genes.

    Returns the mean expression of marker genes in the cell. Higher values
    indicate stronger cell type identity.

    Args:
        row: Expression values for a single cell.
        marker_indices: Column indices of marker genes.

    Returns:
        Mean expression of marker genes.
    """
    if not marker_indices:
        return 0.0

    values = [row[i] for i in marker_indices]
    return sum(values) / len(values)


def _argmax(values: list[float]) -> int:
    """Return index of maximum value in a list.

    Args:
        values: List of numeric values.

    Returns:
        Index of the maximum value. Returns 0 for empty lists.
    """
    if not values:
        return 0
    max_idx = 0
    max_val = values[0]
    for i in range(1, len(values)):
        if values[i] > max_val:
            max_val = values[i]
            max_idx = i
    return max_idx


def _euclidean_distance(a: list[float], b: list[float]) -> float:
    """Compute Euclidean distance between two vectors.

    Args:
        a: First vector.
        b: Second vector (same length as a).

    Returns:
        Euclidean distance.
    """
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


def _simple_cluster(
    rows: list[list[float]],
    max_clusters: int = 10,
) -> list[int]:
    """Simple k-means-style clustering for grouping novel cells.

    Uses a basic iterative assignment algorithm. Not intended as a
    production clustering method; used only for grouping small sets
    of novel cells for summary reporting.

    Args:
        rows: Data rows to cluster.
        max_clusters: Maximum number of clusters.

    Returns:
        List of cluster assignments (0-indexed).
    """
    n = len(rows)
    k = min(max_clusters, n)

    if k <= 1:
        return [0] * n

    # Initialize centroids with first k points
    centroids = [list(rows[i]) for i in range(k)]
    assignments = [0] * n
    n_dims = len(rows[0])

    for _iteration in range(20):
        # Assign each point to nearest centroid
        changed = False
        for i in range(n):
            best_c = 0
            best_dist = _euclidean_distance(rows[i], centroids[0])
            for c in range(1, k):
                d = _euclidean_distance(rows[i], centroids[c])
                if d < best_dist:
                    best_dist = d
                    best_c = c
            if assignments[i] != best_c:
                assignments[i] = best_c
                changed = True

        if not changed:
            break

        # Recompute centroids
        for c in range(k):
            members = [rows[i] for i in range(n) if assignments[i] == c]
            if members:
                centroids[c] = [sum(m[j] for m in members) / len(members) for j in range(n_dims)]

    return assignments


def _wilson_interval(
    successes: int,
    total: int,
    z: float,
) -> tuple[float, float]:
    """Compute Wilson score confidence interval for a binomial proportion.

    Args:
        successes: Number of successes.
        total: Total number of trials.
        z: Z-score for desired confidence level.

    Returns:
        Tuple of (lower_bound, upper_bound) for the proportion.
    """
    if total == 0:
        return (0.0, 0.0)

    p_hat = successes / total
    denom = 1 + z * z / total
    center = (p_hat + z * z / (2 * total)) / denom
    spread = z * math.sqrt((p_hat * (1 - p_hat) + z * z / (4 * total)) / total) / denom

    lower = max(0.0, center - spread)
    upper = min(1.0, center + spread)
    return (lower, upper)


def _z_score_for_confidence(confidence: float) -> float:
    """Return approximate z-score for a given confidence level.

    Args:
        confidence: Confidence level (e.g., 0.95).

    Returns:
        Corresponding z-score.
    """
    # Common z-scores for typical confidence levels
    z_map = {
        0.90: 1.645,
        0.95: 1.960,
        0.99: 2.576,
    }
    if confidence in z_map:
        return z_map[confidence]

    # Approximate using inverse error function approximation
    # For confidence c, z = sqrt(2) * erfinv(c)
    # Using a rational approximation
    alpha = 1.0 - confidence
    if alpha <= 0 or alpha >= 1:
        return 1.96  # fallback

    # Abramowitz and Stegun approximation for probit
    t = math.sqrt(-2.0 * math.log(alpha / 2.0))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    z = t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t)
    return z


__all__ = [
    "annotate_by_markers",
    "cell_type_composition",
    "find_novel_types",
    "score_cell_type",
    "transfer_labels",
]
