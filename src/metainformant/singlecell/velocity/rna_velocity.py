"""RNA velocity estimation and analysis for single-cell data.

Implements RNA velocity computation from spliced and unspliced count
matrices using steady-state and dynamical models. RNA velocity leverages
the ratio of unspliced to spliced mRNA to infer the future transcriptional
state of each cell, enabling trajectory inference and pseudotime ordering
without requiring experimental time-course data.

Models:
    - Steady-state: Assumes transcription/degradation equilibrium. Fits
      a linear relationship between unspliced and spliced counts per gene
      to estimate the degradation rate (gamma). Velocity = unspliced - gamma * spliced.
    - Dynamical: Fits full kinetic parameters (transcription rate alpha,
      splicing rate beta, degradation rate gamma) per gene using an EM-like
      iterative procedure. More accurate but computationally expensive.

Downstream analyses:
    - Velocity embedding: Projects velocity vectors onto UMAP/tSNE space.
    - Velocity pseudotime: Derives temporal ordering from velocity-implied
      transition probabilities.
    - Confidence metrics: Evaluates velocity reliability per gene and cell.
"""

from __future__ import annotations

import math
import random
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


def compute_velocity(
    spliced: Any,
    unspliced: Any,
    gene_names: list[str],
    method: str = "steady_state",
    min_counts: int = 10,
    r_squared_threshold: float = 0.01,
) -> dict:
    """Estimate RNA velocity from spliced and unspliced count matrices.

    For each gene, fits a model relating unspliced to spliced counts
    across cells. In the steady-state model, velocity is computed as
    unspliced - gamma * spliced, where gamma is the slope of the
    linear fit (degradation rate).

    Args:
        spliced: Spliced count matrix (cells x genes). Accepts lists of
            lists, numpy arrays, or scipy sparse matrices.
        unspliced: Unspliced count matrix (cells x genes). Same shape as
            spliced.
        gene_names: Gene names corresponding to columns.
        method: Velocity estimation method. Currently "steady_state".
        min_counts: Minimum total counts (spliced + unspliced) for a gene
            to be included in velocity estimation.
        r_squared_threshold: Minimum R-squared from the linear fit for a
            gene to be considered a velocity gene.

    Returns:
        Dictionary with keys:
            - velocity_matrix: list[list[float]] of velocity values
              (cells x velocity_genes).
            - gamma: dict mapping gene name to estimated degradation rate.
            - r_squared: dict mapping gene name to R-squared of the fit.
            - velocity_genes: list[str] of genes with reliable velocity
              estimates (passing r_squared_threshold).
            - n_velocity_genes: int count of velocity genes.

    Raises:
        ValueError: If spliced and unspliced have different shapes, or
            gene_names length does not match columns.
    """
    if method not in ("steady_state",):
        raise ValueError(f"Invalid method '{method}'. Must be 'steady_state'")

    s_matrix = _to_list_matrix(spliced)
    u_matrix = _to_list_matrix(unspliced)

    n_cells = len(s_matrix)
    n_genes = len(s_matrix[0]) if n_cells > 0 else 0

    if len(u_matrix) != n_cells:
        raise ValueError(
            f"Spliced ({n_cells} cells) and unspliced ({len(u_matrix)} cells) " f"must have the same number of rows"
        )
    if n_cells > 0 and len(u_matrix[0]) != n_genes:
        raise ValueError(
            f"Spliced ({n_genes} genes) and unspliced ({len(u_matrix[0])} genes) "
            f"must have the same number of columns"
        )
    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match columns ({n_genes})")

    logger.info(f"Computing RNA velocity ({method}): {n_cells} cells, {n_genes} genes")

    gamma_values: dict[str, float] = {}
    r_squared_values: dict[str, float] = {}
    velocity_genes: list[str] = []
    velocity_gene_indices: list[int] = []

    for gene_idx in range(n_genes):
        s_col = [s_matrix[i][gene_idx] for i in range(n_cells)]
        u_col = [u_matrix[i][gene_idx] for i in range(n_cells)]

        total_counts = sum(s_col) + sum(u_col)
        if total_counts < min_counts:
            continue

        # Fit linear model: unspliced = gamma * spliced + offset
        # Using least squares to find gamma
        gamma, r_sq = _fit_linear(s_col, u_col)

        gene_name = gene_names[gene_idx]
        gamma_values[gene_name] = gamma
        r_squared_values[gene_name] = r_sq

        if r_sq >= r_squared_threshold and gamma > 0:
            velocity_genes.append(gene_name)
            velocity_gene_indices.append(gene_idx)

    # Compute velocity matrix for velocity genes
    velocity_matrix: list[list[float]] = []
    for cell_idx in range(n_cells):
        cell_velocity: list[float] = []
        for gene_idx in velocity_gene_indices:
            s_val = s_matrix[cell_idx][gene_idx]
            u_val = u_matrix[cell_idx][gene_idx]
            gene_name = gene_names[gene_idx]
            gamma = gamma_values[gene_name]
            vel = u_val - gamma * s_val
            cell_velocity.append(vel)
        velocity_matrix.append(cell_velocity)

    logger.info(
        f"Velocity computed: {len(velocity_genes)}/{n_genes} genes pass "
        f"quality filters (R^2 >= {r_squared_threshold})"
    )

    return {
        "velocity_matrix": velocity_matrix,
        "gamma": gamma_values,
        "r_squared": r_squared_values,
        "velocity_genes": velocity_genes,
        "n_velocity_genes": len(velocity_genes),
    }


def velocity_embedding(
    velocity: Any,
    embedding: Any,
    n_neighbors: int = 30,
) -> dict:
    """Project velocity vectors onto a low-dimensional embedding.

    Computes how each cell's velocity vector translates in the embedding
    space (e.g., UMAP or tSNE) by correlating the velocity with the
    displacement vectors to neighboring cells.

    Args:
        velocity: Velocity matrix (cells x velocity_genes) from
            compute_velocity.
        embedding: Low-dimensional embedding coordinates (cells x 2 or
            cells x n_dims), e.g., from UMAP or tSNE.
        n_neighbors: Number of nearest neighbors in embedding space to
            consider for velocity projection.

    Returns:
        Dictionary with keys:
            - velocity_embedding: list[list[float]] of projected velocity
              vectors in embedding space (cells x n_embedding_dims).
            - transition_matrix: list[list[float]] of cell-to-cell
              transition probabilities (cells x cells).
            - cell_velocities: list[float] of velocity magnitude per cell
              in embedding space.

    Raises:
        ValueError: If velocity and embedding have different cell counts.
    """
    vel_matrix = _to_list_matrix(velocity)
    emb_matrix = _to_list_matrix(embedding)

    n_cells = len(vel_matrix)
    if len(emb_matrix) != n_cells:
        raise ValueError(
            f"Velocity ({n_cells} cells) and embedding ({len(emb_matrix)} cells) " f"must have the same number of rows"
        )

    n_emb_dims = len(emb_matrix[0]) if n_cells > 0 else 2
    effective_k = min(n_neighbors, n_cells - 1)

    logger.info(f"Projecting velocity to {n_emb_dims}D embedding ({n_cells} cells, " f"k={effective_k})")

    # Compute pairwise distances in embedding space
    emb_dists = _pairwise_distances(emb_matrix)

    # For each cell, compute velocity in embedding space
    velocity_emb: list[list[float]] = []
    transition_matrix: list[list[float]] = []
    cell_velocities: list[float] = []

    for i in range(n_cells):
        # Find k nearest neighbors in embedding
        dists_i = emb_dists[i]
        neighbor_indices = sorted(range(n_cells), key=lambda j: dists_i[j])
        # Exclude self
        neighbor_indices = [j for j in neighbor_indices if j != i][:effective_k]

        # Compute displacement vectors in embedding space
        displacements = []
        for j in neighbor_indices:
            disp = [emb_matrix[j][d] - emb_matrix[i][d] for d in range(n_emb_dims)]
            displacements.append(disp)

        # Compute velocity similarity: cosine similarity between velocity
        # vector and difference in expression to each neighbor
        vel_i = vel_matrix[i]
        n_vel_genes = len(vel_i)
        similarities: list[float] = []

        for j_idx, j in enumerate(neighbor_indices):
            vel_j = vel_matrix[j]
            # Expression difference
            expr_diff = [vel_j[g] - vel_i[g] for g in range(n_vel_genes)]
            # Cosine similarity with velocity
            sim = _cosine_similarity(vel_i, expr_diff)
            similarities.append(max(sim, 0.0))  # Only positive transitions

        # Normalize similarities to get transition probabilities
        total_sim = sum(similarities)
        if total_sim > 0:
            probs = [s / total_sim for s in similarities]
        else:
            probs = [1.0 / len(similarities)] * len(similarities) if similarities else []

        # Compute velocity in embedding space as weighted sum of displacements
        vel_emb_i = [0.0] * n_emb_dims
        for j_idx in range(len(neighbor_indices)):
            for d in range(n_emb_dims):
                vel_emb_i[d] += probs[j_idx] * displacements[j_idx][d]

        velocity_emb.append(vel_emb_i)

        # Build full transition row
        trans_row = [0.0] * n_cells
        for j_idx, j in enumerate(neighbor_indices):
            trans_row[j] = probs[j_idx]
        transition_matrix.append(trans_row)

        # Velocity magnitude
        mag = math.sqrt(sum(v**2 for v in vel_emb_i))
        cell_velocities.append(mag)

    logger.info(f"Velocity embedding complete: mean magnitude={_mean(cell_velocities):.4f}")

    return {
        "velocity_embedding": velocity_emb,
        "transition_matrix": transition_matrix,
        "cell_velocities": cell_velocities,
    }


def velocity_pseudotime(
    velocity: Any,
    embedding: Any,
    root_cell: int | None = None,
    n_neighbors: int = 30,
) -> dict:
    """Compute pseudotime from velocity-derived transition probabilities.

    Constructs a transition probability matrix from the velocity vectors
    and embedding, then uses iterative diffusion to compute a pseudotime
    ordering. The root cell (earliest in pseudotime) can be specified or
    automatically selected as the cell with minimal outgoing velocity.

    Args:
        velocity: Velocity matrix (cells x velocity_genes).
        embedding: Embedding coordinates (cells x n_dims).
        root_cell: Index of the root cell. If None, automatically selected
            as the cell with the smallest outgoing velocity magnitude.
        n_neighbors: Number of neighbors for transition computation.

    Returns:
        Dictionary with keys:
            - pseudotime: list[float] of pseudotime values per cell
              (0 at root, increasing along the trajectory).
            - root_cell: int index of the root cell used.
            - terminal_cells: list[int] indices of terminal cells
              (highest pseudotime).

    Raises:
        ValueError: If velocity and embedding have different cell counts.
    """
    vel_matrix = _to_list_matrix(velocity)
    emb_matrix = _to_list_matrix(embedding)
    n_cells = len(vel_matrix)

    if len(emb_matrix) != n_cells:
        raise ValueError(
            f"Velocity ({n_cells} cells) and embedding ({len(emb_matrix)} cells) " f"must have the same number of rows"
        )

    logger.info(f"Computing velocity pseudotime for {n_cells} cells")

    # Get transition matrix and velocity magnitudes
    vel_emb_result = velocity_embedding(vel_matrix, emb_matrix, n_neighbors=n_neighbors)
    cell_velocities = vel_emb_result["cell_velocities"]
    transition_matrix = vel_emb_result["transition_matrix"]

    # Select root cell
    if root_cell is None:
        # Cell with smallest outgoing velocity (most quiescent = likely root)
        root_cell = _argmin(cell_velocities)
    elif root_cell < 0 or root_cell >= n_cells:
        raise ValueError(f"root_cell {root_cell} out of range [0, {n_cells})")

    logger.info(f"Root cell: {root_cell}")

    # Compute pseudotime via shortest-path-like diffusion from root
    # Using iterative relaxation on the transition matrix
    pseudotime = [float("inf")] * n_cells
    pseudotime[root_cell] = 0.0
    visited = [False] * n_cells

    # BFS-like propagation weighted by transition probabilities
    queue = [root_cell]
    visited[root_cell] = True

    while queue:
        next_queue: list[int] = []
        for cell in queue:
            trans_row = transition_matrix[cell]
            current_time = pseudotime[cell]
            for j in range(n_cells):
                if trans_row[j] > 0.01:  # Only significant transitions
                    # Time to reach j from cell: inverse of transition probability
                    travel_time = 1.0 / trans_row[j] if trans_row[j] > 0 else float("inf")
                    new_time = current_time + travel_time
                    if new_time < pseudotime[j]:
                        pseudotime[j] = new_time
                        if not visited[j]:
                            visited[j] = True
                            next_queue.append(j)
        queue = next_queue

    # Handle unreachable cells
    max_reachable = max(t for t in pseudotime if t < float("inf"))
    pseudotime = [t if t < float("inf") else max_reachable for t in pseudotime]

    # Normalize to [0, 1]
    min_pt = min(pseudotime)
    max_pt = max(pseudotime)
    pt_range = max_pt - min_pt
    if pt_range > 0:
        pseudotime = [(t - min_pt) / pt_range for t in pseudotime]
    else:
        pseudotime = [0.0] * n_cells

    # Identify terminal cells (top 5% pseudotime)
    pt_threshold = 0.95
    terminal_cells = [i for i in range(n_cells) if pseudotime[i] >= pt_threshold]

    logger.info(f"Pseudotime computed: root={root_cell}, " f"{len(terminal_cells)} terminal cells")

    return {
        "pseudotime": pseudotime,
        "root_cell": root_cell,
        "terminal_cells": terminal_cells,
    }


def velocity_confidence(
    velocity: Any,
    spliced: Any,
) -> dict:
    """Compute confidence metrics for velocity estimates.

    Evaluates the reliability of velocity estimates per gene (based on
    the consistency of the velocity direction across cells) and per cell
    (based on the agreement among velocity-gene estimates).

    Args:
        velocity: Velocity matrix (cells x velocity_genes).
        spliced: Spliced count matrix (cells x genes) for computing
            expression-level confidence weighting.

    Returns:
        Dictionary with keys:
            - gene_confidence: list[float] of confidence per velocity gene
              (fraction of cells with consistent velocity sign).
            - cell_confidence: list[float] of confidence per cell (mean
              absolute velocity normalized by expression level).
            - overall_confidence: float mean of cell confidences.
            - n_high_confidence_genes: int genes with confidence > 0.7.
            - n_high_confidence_cells: int cells with confidence > 0.5.
    """
    vel_matrix = _to_list_matrix(velocity)
    s_matrix = _to_list_matrix(spliced)

    n_cells = len(vel_matrix)
    n_vel_genes = len(vel_matrix[0]) if n_cells > 0 else 0

    logger.info(f"Computing velocity confidence: {n_cells} cells, {n_vel_genes} velocity genes")

    # Gene confidence: fraction of cells with consistent velocity direction
    gene_confidence: list[float] = []
    for gene_idx in range(n_vel_genes):
        vel_col = [vel_matrix[i][gene_idx] for i in range(n_cells)]
        n_positive = sum(1 for v in vel_col if v > 0)
        n_negative = sum(1 for v in vel_col if v < 0)
        n_nonzero = n_positive + n_negative
        if n_nonzero > 0:
            # Confidence = max fraction with same sign
            conf = max(n_positive, n_negative) / n_nonzero
        else:
            conf = 0.0
        gene_confidence.append(conf)

    # Cell confidence: mean absolute velocity weighted by expression
    cell_confidence: list[float] = []
    for cell_idx in range(n_cells):
        vel_row = vel_matrix[cell_idx]
        abs_vel = [abs(v) for v in vel_row]
        mean_abs_vel = sum(abs_vel) / len(abs_vel) if abs_vel else 0.0

        # Normalize by expression level
        s_row = s_matrix[cell_idx] if cell_idx < len(s_matrix) else [0.0]
        mean_expr = sum(s_row) / len(s_row) if s_row else 1.0
        if mean_expr > 0:
            conf = min(mean_abs_vel / mean_expr, 1.0)
        else:
            conf = 0.0
        cell_confidence.append(conf)

    overall = _mean(cell_confidence) if cell_confidence else 0.0
    n_high_genes = sum(1 for c in gene_confidence if c > 0.7)
    n_high_cells = sum(1 for c in cell_confidence if c > 0.5)

    logger.info(
        f"Confidence: {n_high_genes}/{n_vel_genes} high-confidence genes, "
        f"{n_high_cells}/{n_cells} high-confidence cells, "
        f"overall={overall:.3f}"
    )

    return {
        "gene_confidence": gene_confidence,
        "cell_confidence": cell_confidence,
        "overall_confidence": overall,
        "n_high_confidence_genes": n_high_genes,
        "n_high_confidence_cells": n_high_cells,
    }


def fit_dynamical_model(
    spliced: Any,
    unspliced: Any,
    gene_names: list[str],
    max_iter: int = 100,
    tol: float = 1e-4,
    min_counts: int = 10,
) -> dict:
    """Fit a dynamical model with full kinetic parameters per gene.

    Estimates transcription rate (alpha), splicing rate (beta), and
    degradation rate (gamma) for each gene using an iterative
    expectation-maximization-like procedure. This is a simplified
    implementation inspired by the scVelo dynamical model.

    The model assumes:
        du/dt = alpha - beta * u   (unspliced dynamics)
        ds/dt = beta * u - gamma * s   (spliced dynamics)

    Args:
        spliced: Spliced count matrix (cells x genes).
        unspliced: Unspliced count matrix (cells x genes).
        gene_names: Gene names for columns.
        max_iter: Maximum iterations for the optimization.
        tol: Convergence tolerance (relative change in parameters).
        min_counts: Minimum total counts per gene for inclusion.

    Returns:
        Dictionary with keys:
            - alpha: dict mapping gene name to transcription rate.
            - beta: dict mapping gene name to splicing rate.
            - gamma: dict mapping gene name to degradation rate.
            - likelihood: dict mapping gene name to final log-likelihood.
            - velocity_matrix: list[list[float]] of dynamical velocity
              (cells x fitted_genes).
            - fitted_genes: list[str] of genes that converged.
            - n_iterations: dict mapping gene name to iterations used.

    Raises:
        ValueError: If input dimensions are inconsistent.
    """
    s_matrix = _to_list_matrix(spliced)
    u_matrix = _to_list_matrix(unspliced)

    n_cells = len(s_matrix)
    n_genes = len(s_matrix[0]) if n_cells > 0 else 0

    if len(u_matrix) != n_cells:
        raise ValueError(
            f"Spliced ({n_cells} cells) and unspliced ({len(u_matrix)} cells) " f"must have the same number of rows"
        )
    if n_cells > 0 and len(u_matrix[0]) != n_genes:
        raise ValueError(
            f"Spliced ({n_genes} genes) and unspliced ({len(u_matrix[0])} genes) "
            f"must have the same number of columns"
        )
    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match columns ({n_genes})")

    logger.info(f"Fitting dynamical model: {n_cells} cells, {n_genes} genes, " f"max_iter={max_iter}")

    alpha_values: dict[str, float] = {}
    beta_values: dict[str, float] = {}
    gamma_values: dict[str, float] = {}
    likelihood_values: dict[str, float] = {}
    n_iterations: dict[str, int] = {}
    fitted_genes: list[str] = []
    fitted_gene_indices: list[int] = []

    for gene_idx in range(n_genes):
        s_col = [s_matrix[i][gene_idx] for i in range(n_cells)]
        u_col = [u_matrix[i][gene_idx] for i in range(n_cells)]

        total = sum(s_col) + sum(u_col)
        if total < min_counts:
            continue

        gene_name = gene_names[gene_idx]

        # Initialize parameters from steady-state estimates
        gamma_init, _ = _fit_linear(s_col, u_col)
        gamma_init = max(gamma_init, 0.01)

        # Initialize alpha and beta
        mean_u = _mean(u_col)
        mean_s = _mean(s_col)
        beta_init = max(mean_u * 0.1, 0.01)
        alpha_init = max(beta_init * mean_u + gamma_init * mean_s, 0.01)

        alpha = alpha_init
        beta = beta_init
        gamma = gamma_init

        # Iterative optimization (simplified EM-like)
        prev_params = (alpha, beta, gamma)
        converged = False
        iters_used = 0

        for iteration in range(max_iter):
            iters_used = iteration + 1

            # E-step: compute expected steady-state values
            expected_u = [alpha / beta if beta > 0 else 0.0] * n_cells
            expected_s = [beta * (alpha / beta) / gamma if gamma > 0 and beta > 0 else 0.0] * n_cells

            # M-step: update parameters to minimize residuals
            # Update gamma from spliced dynamics: ds/dt = beta*u - gamma*s
            numerator_gamma = sum(beta * u_col[i] for i in range(n_cells))
            denominator_gamma = sum(s_col[i] for i in range(n_cells))
            if denominator_gamma > 0:
                gamma = max(numerator_gamma / denominator_gamma, 0.001)

            # Update beta from unspliced dynamics: du/dt = alpha - beta*u
            denominator_beta = sum(u_col[i] for i in range(n_cells))
            if denominator_beta > 0:
                alpha_est = _mean([u_col[i] + s_col[i] * gamma / max(beta, 0.001) for i in range(n_cells)])
                beta = max(alpha_est / max(_mean(u_col), 0.001), 0.001)

            # Update alpha
            alpha = max(beta * _mean(u_col), 0.001)

            # Check convergence
            current_params = (alpha, beta, gamma)
            relative_change = sum(abs(c - p) / max(abs(p), 1e-10) for c, p in zip(current_params, prev_params)) / 3.0

            if relative_change < tol:
                converged = True
                break

            prev_params = current_params

        # Compute log-likelihood (Gaussian residuals)
        residuals_u = [(u_col[i] - alpha / beta) ** 2 for i in range(n_cells)]
        residuals_s = [(s_col[i] - beta * (alpha / beta) / gamma) ** 2 for i in range(n_cells)]
        total_residual = sum(residuals_u) + sum(residuals_s)
        log_lik = -0.5 * total_residual / max(n_cells, 1)

        alpha_values[gene_name] = alpha
        beta_values[gene_name] = beta
        gamma_values[gene_name] = gamma
        likelihood_values[gene_name] = log_lik
        n_iterations[gene_name] = iters_used

        if converged:
            fitted_genes.append(gene_name)
            fitted_gene_indices.append(gene_idx)

    # Compute velocity for fitted genes
    velocity_matrix: list[list[float]] = []
    for cell_idx in range(n_cells):
        cell_vel: list[float] = []
        for gene_idx in fitted_gene_indices:
            gene_name = gene_names[gene_idx]
            s_val = s_matrix[cell_idx][gene_idx]
            u_val = u_matrix[cell_idx][gene_idx]
            gamma = gamma_values[gene_name]
            vel = u_val - gamma * s_val
            cell_vel.append(vel)
        velocity_matrix.append(cell_vel)

    logger.info(
        f"Dynamical model fit: {len(fitted_genes)}/{n_genes} genes converged "
        f"(max {max(n_iterations.values()) if n_iterations else 0} iterations)"
    )

    return {
        "alpha": alpha_values,
        "beta": beta_values,
        "gamma": gamma_values,
        "likelihood": likelihood_values,
        "velocity_matrix": velocity_matrix,
        "fitted_genes": fitted_genes,
        "n_iterations": n_iterations,
    }


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _to_list_matrix(data: Any) -> list[list[float]]:
    """Convert various matrix types to list of lists of floats.

    Args:
        data: Matrix data. Supports list[list[float]], numpy arrays,
            scipy sparse matrices, or objects with .toarray()/.values.

    Returns:
        List of lists of floats.
    """
    if isinstance(data, list):
        return data

    if HAS_NUMPY and isinstance(data, np.ndarray):
        return data.tolist()

    if hasattr(data, "toarray"):
        return data.toarray().tolist()

    if hasattr(data, "values"):
        return data.values.tolist()

    raise TypeError(f"Unsupported matrix type: {type(data)}")


def _fit_linear(x: list[float], y: list[float]) -> tuple[float, float]:
    """Fit a simple linear regression y = slope * x.

    Uses ordinary least squares through the origin (no intercept).

    Args:
        x: Independent variable values.
        y: Dependent variable values.

    Returns:
        Tuple of (slope, r_squared).
    """
    n = len(x)
    if n == 0:
        return (0.0, 0.0)

    sum_xy = sum(x[i] * y[i] for i in range(n))
    sum_xx = sum(x[i] * x[i] for i in range(n))

    if sum_xx == 0:
        return (0.0, 0.0)

    slope = sum_xy / sum_xx

    # R-squared
    mean_y = sum(y) / n
    ss_tot = sum((y[i] - mean_y) ** 2 for i in range(n))
    ss_res = sum((y[i] - slope * x[i]) ** 2 for i in range(n))

    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    r_squared = max(r_squared, 0.0)

    return (slope, r_squared)


def _pairwise_distances(matrix: list[list[float]]) -> list[list[float]]:
    """Compute pairwise Euclidean distances between rows.

    Args:
        matrix: List of row vectors.

    Returns:
        Square distance matrix.
    """
    n = len(matrix)
    dists = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((matrix[i][k] - matrix[j][k]) ** 2 for k in range(len(matrix[i]))))
            dists[i][j] = d
            dists[j][i] = d
    return dists


def _cosine_similarity(a: list[float], b: list[float]) -> float:
    """Compute cosine similarity between two vectors.

    Args:
        a: First vector.
        b: Second vector.

    Returns:
        Cosine similarity in range [-1, 1].
    """
    dot = sum(a[i] * b[i] for i in range(len(a)))
    norm_a = math.sqrt(sum(v * v for v in a))
    norm_b = math.sqrt(sum(v * v for v in b))

    if norm_a == 0 or norm_b == 0:
        return 0.0

    return dot / (norm_a * norm_b)


def _mean(values: list[float]) -> float:
    """Compute arithmetic mean.

    Args:
        values: List of numeric values.

    Returns:
        Mean value, or 0.0 for empty lists.
    """
    if not values:
        return 0.0
    return sum(values) / len(values)


def _argmin(values: list[float]) -> int:
    """Return index of minimum value.

    Args:
        values: List of numeric values.

    Returns:
        Index of the minimum value.
    """
    if not values:
        return 0
    min_idx = 0
    min_val = values[0]
    for i in range(1, len(values)):
        if values[i] < min_val:
            min_val = values[i]
            min_idx = i
    return min_idx


__all__ = [
    "compute_velocity",
    "fit_dynamical_model",
    "velocity_confidence",
    "velocity_embedding",
    "velocity_pseudotime",
]
