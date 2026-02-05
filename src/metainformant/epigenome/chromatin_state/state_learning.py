"""Chromatin state learning and analysis.

This module implements a simplified ChromHMM-style approach to chromatin
state discovery from histone modification data.  It uses multivariate
Gaussian mixture models fit via the Expectation-Maximization algorithm,
Viterbi-like state assignment, biological state interpretation based on
known histone mark patterns, enrichment testing against genomic
annotations, genome segmentation, and cross-condition comparison.
"""

from __future__ import annotations

import math
import random
import statistics
from collections import defaultdict
from typing import Any, Optional

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional numpy for vectorized operations
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

# Optional scipy for statistical tests
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _to_2d_list(signal_matrix: list[list[float]] | Any) -> list[list[float]]:
    """Normalise signal_matrix into a plain list-of-lists.

    Accepts numpy arrays, nested lists, or any array-like, and returns
    a consistent ``list[list[float]]`` representation.

    Args:
        signal_matrix: Input matrix (rows = bins, cols = marks).

    Returns:
        2D list of floats.

    Raises:
        ValueError: If the input is empty or has inconsistent dimensions.
    """
    if HAS_NUMPY and isinstance(signal_matrix, np.ndarray):
        if signal_matrix.ndim != 2:
            raise ValueError(f"Expected 2D array, got {signal_matrix.ndim}D")
        return signal_matrix.tolist()

    if not signal_matrix:
        raise ValueError("signal_matrix must not be empty")

    result: list[list[float]] = []
    expected_cols: int | None = None
    for row in signal_matrix:
        converted = [float(v) for v in row]
        if expected_cols is None:
            expected_cols = len(converted)
        elif len(converted) != expected_cols:
            raise ValueError(f"Inconsistent row lengths: expected {expected_cols}, " f"got {len(converted)}")
        result.append(converted)

    if expected_cols is None or expected_cols == 0:
        raise ValueError("signal_matrix must have at least one column")
    return result


def _multivariate_gaussian_logpdf(
    x: list[float],
    mean: list[float],
    variance: list[float],
) -> float:
    """Compute log-PDF of a diagonal multivariate Gaussian.

    Uses a diagonal covariance model (independent features) for
    computational efficiency.

    Args:
        x: Observation vector.
        mean: Mean vector.
        variance: Diagonal variance vector (one per dimension).

    Returns:
        Log probability density.
    """
    d = len(x)
    log_pdf = -0.5 * d * math.log(2.0 * math.pi)
    for j in range(d):
        var_j = max(variance[j], 1e-10)
        log_pdf -= 0.5 * math.log(var_j)
        log_pdf -= 0.5 * ((x[j] - mean[j]) ** 2) / var_j
    return log_pdf


def _logsumexp(values: list[float]) -> float:
    """Numerically stable log-sum-exp.

    Args:
        values: List of log-space values.

    Returns:
        log(sum(exp(values))).
    """
    if not values:
        return float("-inf")
    max_val = max(values)
    if max_val == float("-inf"):
        return float("-inf")
    return max_val + math.log(sum(math.exp(v - max_val) for v in values))


def _kmeans_init(
    data: list[list[float]],
    n_clusters: int,
    max_iter: int = 20,
) -> list[list[float]]:
    """Simple k-means++ initialisation to get starting centroids.

    Args:
        data: 2D list of observations.
        n_clusters: Number of clusters.
        max_iter: Maximum k-means iterations.

    Returns:
        List of centroid vectors.
    """
    n = len(data)
    d = len(data[0])

    if n_clusters >= n:
        return data[:n_clusters]

    # k-means++ initialisation
    rng = random.Random(42)
    centroids: list[list[float]] = [data[rng.randint(0, n - 1)]]

    for _ in range(1, n_clusters):
        # Compute distances to nearest centroid
        dists = []
        for point in data:
            min_dist = float("inf")
            for c in centroids:
                dist = sum((point[j] - c[j]) ** 2 for j in range(d))
                min_dist = min(min_dist, dist)
            dists.append(min_dist)

        total_dist = sum(dists)
        if total_dist <= 0:
            centroids.append(data[rng.randint(0, n - 1)])
            continue

        # Weighted random selection
        threshold = rng.random() * total_dist
        cumulative = 0.0
        for idx, dist in enumerate(dists):
            cumulative += dist
            if cumulative >= threshold:
                centroids.append(data[idx])
                break

    # Run a few k-means iterations to refine
    assignments = [0] * n
    for _ in range(max_iter):
        # Assign step
        for i in range(n):
            best_k = 0
            best_dist = float("inf")
            for k in range(len(centroids)):
                dist = sum((data[i][j] - centroids[k][j]) ** 2 for j in range(d))
                if dist < best_dist:
                    best_dist = dist
                    best_k = k
            assignments[i] = best_k

        # Update step
        new_centroids: list[list[float]] = []
        for k in range(len(centroids)):
            members = [data[i] for i in range(n) if assignments[i] == k]
            if members:
                centroid = [sum(m[j] for m in members) / len(members) for j in range(d)]
                new_centroids.append(centroid)
            else:
                new_centroids.append(centroids[k])
        centroids = new_centroids

    return centroids


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def learn_chromatin_states(
    signal_matrix: list[list[float]] | Any,
    n_states: int = 15,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> dict:
    """Learn chromatin states using a multivariate Gaussian mixture model.

    Implements a simplified ChromHMM-style approach using the EM
    algorithm with diagonal covariance.  Rows of *signal_matrix*
    correspond to genomic bins and columns to histone marks.

    Args:
        signal_matrix: Matrix of shape (n_bins, n_marks).  Accepts
            numpy arrays or nested lists.
        n_states: Number of chromatin states to learn.
        max_iter: Maximum EM iterations.
        tol: Convergence tolerance on log-likelihood change.

    Returns:
        Dictionary with keys:

        - ``states``: Number of states.
        - ``emission_params``: Dict with ``means`` (list[list[float]])
          and ``variances`` (list[list[float]]) per state.
        - ``transition_matrix``: State-to-state transition probabilities
          (list[list[float]]).
        - ``assignments``: State assignment for each bin (list[int]).
        - ``log_likelihood``: Final log-likelihood.
        - ``n_iterations``: Number of EM iterations performed.
        - ``converged``: Whether the algorithm converged.

    Raises:
        ValueError: If signal_matrix is empty, has fewer rows than
            states, or has inconsistent dimensions.
    """
    data = _to_2d_list(signal_matrix)
    n_bins = len(data)
    n_marks = len(data[0])

    if n_bins < n_states:
        raise ValueError(f"Need at least n_states ({n_states}) bins, got {n_bins}")

    logger.info(
        "Learning %d chromatin states from %d bins x %d marks " "(max_iter=%d, tol=%.1e)",
        n_states,
        n_bins,
        n_marks,
        max_iter,
        tol,
    )

    # Initialise parameters via k-means
    centroids = _kmeans_init(data, n_states)

    # Emission parameters: per-state mean and diagonal variance
    means: list[list[float]] = [list(c) for c in centroids]
    variances: list[list[float]] = []
    for k in range(n_states):
        # Initial variance: global variance per mark
        col_vars = []
        for j in range(n_marks):
            col_vals = [data[i][j] for i in range(n_bins)]
            if len(col_vals) > 1:
                col_var = statistics.variance(col_vals)
            else:
                col_var = 1.0
            col_vars.append(max(col_var, 1e-6))
        variances.append(col_vars)

    # Mixing weights (uniform init)
    weights = [1.0 / n_states] * n_states

    # Transition matrix (uniform init with self-loop bias)
    transition: list[list[float]] = []
    for k in range(n_states):
        row = [0.1 / (n_states - 1) if j != k else 0.9 for j in range(n_states)]
        row_sum = sum(row)
        transition.append([v / row_sum for v in row])

    # EM algorithm
    prev_ll = float("-inf")
    converged = False
    iteration = 0

    # Responsibilities matrix: n_bins x n_states
    responsibilities: list[list[float]] = [[0.0] * n_states for _ in range(n_bins)]

    for iteration in range(1, max_iter + 1):
        # --- E-step ---
        log_likelihood = 0.0

        for i in range(n_bins):
            log_probs = []
            for k in range(n_states):
                log_prior = math.log(max(weights[k], 1e-300))
                log_emission = _multivariate_gaussian_logpdf(data[i], means[k], variances[k])
                log_probs.append(log_prior + log_emission)

            log_norm = _logsumexp(log_probs)
            log_likelihood += log_norm

            for k in range(n_states):
                responsibilities[i][k] = math.exp(log_probs[k] - log_norm)

        # Check convergence
        ll_change = abs(log_likelihood - prev_ll)
        if iteration > 1 and ll_change < tol:
            converged = True
            logger.info(
                "EM converged at iteration %d (delta=%.2e)",
                iteration,
                ll_change,
            )
            break

        if iteration % 20 == 0:
            logger.debug(
                "EM iteration %d: log_likelihood=%.4f, delta=%.4e",
                iteration,
                log_likelihood,
                ll_change,
            )

        prev_ll = log_likelihood

        # --- M-step ---
        for k in range(n_states):
            n_k = sum(responsibilities[i][k] for i in range(n_bins))
            n_k = max(n_k, 1e-10)

            # Update weight
            weights[k] = n_k / n_bins

            # Update mean
            for j in range(n_marks):
                means[k][j] = sum(responsibilities[i][k] * data[i][j] for i in range(n_bins)) / n_k

            # Update variance
            for j in range(n_marks):
                var_j = sum(responsibilities[i][k] * (data[i][j] - means[k][j]) ** 2 for i in range(n_bins)) / n_k
                variances[k][j] = max(var_j, 1e-6)

        # Update transition matrix from consecutive bin assignments
        for k in range(n_states):
            row_sum_trans = 0.0
            for j in range(n_states):
                count = sum(responsibilities[i][k] * responsibilities[i + 1][j] for i in range(n_bins - 1))
                transition[k][j] = count
                row_sum_trans += count
            if row_sum_trans > 0:
                transition[k] = [v / row_sum_trans for v in transition[k]]

    # Final state assignments (MAP)
    assignments = []
    for i in range(n_bins):
        best_state = 0
        best_resp = responsibilities[i][0]
        for k in range(1, n_states):
            if responsibilities[i][k] > best_resp:
                best_resp = responsibilities[i][k]
                best_state = k
        assignments.append(best_state)

    result = {
        "states": n_states,
        "emission_params": {
            "means": means,
            "variances": variances,
        },
        "transition_matrix": transition,
        "assignments": assignments,
        "log_likelihood": round(log_likelihood, 4),
        "n_iterations": iteration,
        "converged": converged,
        "weights": [round(w, 6) for w in weights],
    }

    # Count per-state assignments
    state_counts: dict[int, int] = defaultdict(int)
    for s in assignments:
        state_counts[s] += 1
    result["state_counts"] = dict(state_counts)

    logger.info(
        "Learned %d states in %d iterations (converged=%s, LL=%.4f)",
        n_states,
        iteration,
        converged,
        log_likelihood,
    )
    return result


def assign_states(
    signal_matrix: list[list[float]] | Any,
    model: dict,
) -> list[int]:
    """Assign chromatin states to genomic bins using a learned model.

    Performs Viterbi-like decoding using the emission parameters and
    transition matrix from the model to find the most likely state
    sequence.

    Args:
        signal_matrix: Matrix of shape (n_bins, n_marks).
        model: Model dictionary as returned by :func:`learn_chromatin_states`.

    Returns:
        List of state assignments (one per bin).

    Raises:
        ValueError: If signal_matrix is empty or has wrong number of marks.
    """
    data = _to_2d_list(signal_matrix)
    n_bins = len(data)

    emission = model["emission_params"]
    means = emission["means"]
    variances = emission["variances"]
    n_states = model["states"]
    transition = model["transition_matrix"]
    weights = model.get("weights", [1.0 / n_states] * n_states)

    n_marks = len(means[0])
    if len(data[0]) != n_marks:
        raise ValueError(f"Data has {len(data[0])} marks but model expects {n_marks}")

    logger.info("Assigning states to %d bins using %d-state model", n_bins, n_states)

    # Viterbi algorithm
    # V[i][k] = log probability of most likely path ending in state k at bin i
    V: list[list[float]] = [[0.0] * n_states for _ in range(n_bins)]
    backpointer: list[list[int]] = [[0] * n_states for _ in range(n_bins)]

    # Initialise
    for k in range(n_states):
        log_pi = math.log(max(weights[k], 1e-300))
        log_emit = _multivariate_gaussian_logpdf(data[0], means[k], variances[k])
        V[0][k] = log_pi + log_emit

    # Forward pass
    for i in range(1, n_bins):
        for k in range(n_states):
            log_emit = _multivariate_gaussian_logpdf(data[i], means[k], variances[k])
            best_prev = float("-inf")
            best_prev_state = 0

            for j in range(n_states):
                log_trans = math.log(max(transition[j][k], 1e-300))
                score = V[i - 1][j] + log_trans
                if score > best_prev:
                    best_prev = score
                    best_prev_state = j

            V[i][k] = best_prev + log_emit
            backpointer[i][k] = best_prev_state

    # Backtrace
    assignments = [0] * n_bins
    # Find best final state
    best_final = 0
    best_score = V[n_bins - 1][0]
    for k in range(1, n_states):
        if V[n_bins - 1][k] > best_score:
            best_score = V[n_bins - 1][k]
            best_final = k
    assignments[n_bins - 1] = best_final

    for i in range(n_bins - 2, -1, -1):
        assignments[i] = backpointer[i + 1][assignments[i + 1]]

    state_counts: dict[int, int] = defaultdict(int)
    for s in assignments:
        state_counts[s] += 1

    logger.info(
        "Assigned states to %d bins (%d unique states used)",
        n_bins,
        len(state_counts),
    )
    return assignments


def interpret_states(
    emission_params: dict,
    mark_names: list[str],
) -> list[dict]:
    """Automatically label chromatin states based on emission patterns.

    Uses known associations between histone marks and functional
    categories:

    - H3K4me3 high -> Active Promoter
    - H3K4me1 high + H3K27ac high -> Active Enhancer
    - H3K4me1 high + H3K27ac low -> Poised Enhancer
    - H3K27me3 high -> Polycomb Repressed
    - H3K9me3 high -> Heterochromatin
    - H3K36me3 high -> Transcribed / Gene Body
    - H3K27ac high (without H3K4me1/H3K4me3) -> Active Regulatory
    - Low signal across all marks -> Quiescent

    Args:
        emission_params: Dictionary with ``means`` key containing
            per-state mean vectors.
        mark_names: List of histone mark names matching the column
            order (e.g. ``["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3",
            "H3K9me3", "H3K36me3"]``).

    Returns:
        List of dictionaries (one per state) with keys: ``state_id``,
        ``label``, ``category``, ``dominant_marks``, ``mean_signal``,
        ``description``.
    """
    means = emission_params["means"]
    n_states = len(means)
    n_marks = len(mark_names)

    logger.info(
        "Interpreting %d chromatin states with %d marks: %s",
        n_states,
        n_marks,
        mark_names,
    )

    # Normalise mark names to lowercase for pattern matching
    mark_lower = [m.lower() for m in mark_names]

    # Build mark index for known marks
    def _find_mark(pattern: str) -> int | None:
        """Find index of mark matching pattern (case-insensitive substring)."""
        for idx, name in enumerate(mark_lower):
            if pattern.lower() in name:
                return idx
        return None

    idx_k4me3 = _find_mark("k4me3")
    idx_k4me1 = _find_mark("k4me1")
    idx_k27ac = _find_mark("k27ac")
    idx_k27me3 = _find_mark("k27me3")
    idx_k9me3 = _find_mark("k9me3")
    idx_k36me3 = _find_mark("k36me3")

    # Compute per-mark thresholds (mean + 0.5*std across states)
    mark_thresholds: list[float] = []
    for j in range(n_marks):
        vals = [means[k][j] for k in range(n_states)]
        m = sum(vals) / len(vals)
        if len(vals) > 1:
            std = (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
        else:
            std = 0.0
        mark_thresholds.append(m + 0.5 * std)

    interpretations: list[dict] = []

    for k in range(n_states):
        state_means = means[k]

        # Determine which marks are "high" for this state
        high_marks: list[str] = []
        for j in range(n_marks):
            if state_means[j] >= mark_thresholds[j]:
                high_marks.append(mark_names[j])

        # Classify state based on mark patterns
        label = "Unknown"
        category = "other"
        description = "Uncharacterised chromatin state"

        has_k4me3 = idx_k4me3 is not None and state_means[idx_k4me3] >= mark_thresholds[idx_k4me3]
        has_k4me1 = idx_k4me1 is not None and state_means[idx_k4me1] >= mark_thresholds[idx_k4me1]
        has_k27ac = idx_k27ac is not None and state_means[idx_k27ac] >= mark_thresholds[idx_k27ac]
        has_k27me3 = idx_k27me3 is not None and state_means[idx_k27me3] >= mark_thresholds[idx_k27me3]
        has_k9me3 = idx_k9me3 is not None and state_means[idx_k9me3] >= mark_thresholds[idx_k9me3]
        has_k36me3 = idx_k36me3 is not None and state_means[idx_k36me3] >= mark_thresholds[idx_k36me3]

        if has_k4me3 and has_k27me3:
            label = "Bivalent Promoter"
            category = "bivalent"
            description = (
                "Bivalent chromatin with both activating (H3K4me3) and "
                "repressive (H3K27me3) marks; common in developmental genes"
            )
        elif has_k4me3 and has_k27ac:
            label = "Active Promoter"
            category = "promoter"
            description = "Active promoter marked by H3K4me3 and H3K27ac; " "associated with actively transcribed genes"
        elif has_k4me3:
            label = "Promoter"
            category = "promoter"
            description = "Promoter region marked by H3K4me3"
        elif has_k4me1 and has_k27ac:
            label = "Active Enhancer"
            category = "enhancer"
            description = (
                "Active enhancer marked by H3K4me1 and H3K27ac; " "drives gene expression in the current cell state"
            )
        elif has_k4me1 and has_k27me3:
            label = "Poised Enhancer"
            category = "enhancer"
            description = "Poised enhancer marked by H3K4me1 and H3K27me3; " "primed for activation upon signal"
        elif has_k4me1:
            label = "Weak Enhancer"
            category = "enhancer"
            description = "Weak or primed enhancer marked by H3K4me1 alone"
        elif has_k27me3:
            label = "Polycomb Repressed"
            category = "repressed"
            description = "Polycomb-repressed region marked by H3K27me3; " "silenced by Polycomb group proteins"
        elif has_k9me3:
            label = "Heterochromatin"
            category = "repressed"
            description = "Constitutive heterochromatin marked by H3K9me3; " "permanently silenced regions"
        elif has_k36me3:
            label = "Transcribed"
            category = "transcribed"
            description = "Actively transcribed gene body marked by H3K36me3"
        elif has_k27ac:
            label = "Active Regulatory"
            category = "regulatory"
            description = "Active regulatory region marked by H3K27ac"
        elif not high_marks:
            label = "Quiescent"
            category = "quiescent"
            description = "Quiescent or low-signal region with no strong " "histone modifications"

        mean_signal = round(sum(state_means) / n_marks, 4) if n_marks > 0 else 0.0

        interpretations.append(
            {
                "state_id": k,
                "label": label,
                "category": category,
                "dominant_marks": high_marks,
                "mean_signal": mean_signal,
                "description": description,
                "emission_means": [round(v, 4) for v in state_means],
            }
        )

    logger.info(
        "Interpreted %d states: %s",
        n_states,
        [s["label"] for s in interpretations],
    )
    return interpretations


def compute_state_enrichment(
    state_assignments: list[int],
    annotations: dict,
    n_states: int,
) -> dict:
    """Test enrichment of chromatin states in genomic annotations.

    For each state and each annotation category, computes the observed
    vs expected overlap and a Fisher's exact test (or chi-squared
    approximation as fallback) for enrichment.

    Args:
        state_assignments: Per-bin state assignments.
        annotations: Dictionary mapping annotation names to sets or
            lists of bin indices that belong to that annotation
            (e.g. ``{"TSS": [0, 1, 5, 10], "gene_body": [2, 3, 4]}``).
        n_states: Total number of chromatin states.

    Returns:
        Dictionary with keys:

        - ``enrichment_matrix``: Dict mapping
          ``(state_id, annotation)`` to enrichment fold and p-value.
        - ``summary``: Per-state list of most enriched annotation.
    """
    n_bins = len(state_assignments)
    if n_bins == 0:
        return {"enrichment_matrix": {}, "summary": []}

    logger.info(
        "Computing state enrichment: %d bins, %d states, %d annotations",
        n_bins,
        n_states,
        len(annotations),
    )

    # Build per-state bin sets
    state_bins: dict[int, set[int]] = defaultdict(set)
    for i, s in enumerate(state_assignments):
        state_bins[s].add(i)

    # Build annotation bin sets
    anno_bins: dict[str, set[int]] = {}
    for name, indices in annotations.items():
        anno_bins[name] = set(indices)

    enrichment_matrix: dict[str, dict] = {}
    summary: list[dict] = []

    for state_id in range(n_states):
        s_bins = state_bins.get(state_id, set())
        n_state = len(s_bins)

        best_enrichment: dict | None = None

        for anno_name, a_bins in anno_bins.items():
            n_anno = len(a_bins)

            # 2x2 contingency table
            overlap = len(s_bins & a_bins)
            state_only = n_state - overlap
            anno_only = n_anno - overlap
            neither = n_bins - n_state - n_anno + overlap

            # Expected overlap under independence
            expected = (n_state * n_anno) / n_bins if n_bins > 0 else 0
            fold = overlap / expected if expected > 0 else 0.0

            # Fisher's exact test or chi-squared approximation
            pval = _fisher_or_chi2(overlap, state_only, anno_only, neither)

            key = f"state_{state_id}__{anno_name}"
            enrichment_matrix[key] = {
                "state_id": state_id,
                "annotation": anno_name,
                "observed_overlap": overlap,
                "expected_overlap": round(expected, 4),
                "fold_enrichment": round(fold, 4),
                "p_value": pval,
                "n_state_bins": n_state,
                "n_annotation_bins": n_anno,
            }

            if best_enrichment is None or fold > best_enrichment["fold_enrichment"]:
                best_enrichment = enrichment_matrix[key]

        summary.append(
            {
                "state_id": state_id,
                "n_bins": n_state,
                "most_enriched_annotation": (best_enrichment["annotation"] if best_enrichment else None),
                "best_fold_enrichment": (best_enrichment["fold_enrichment"] if best_enrichment else 0.0),
                "best_p_value": (best_enrichment["p_value"] if best_enrichment else 1.0),
            }
        )

    result = {
        "enrichment_matrix": enrichment_matrix,
        "summary": summary,
        "n_bins": n_bins,
        "n_states": n_states,
        "n_annotations": len(annotations),
    }

    logger.info("Computed enrichment for %d state-annotation pairs", len(enrichment_matrix))
    return result


def _fisher_or_chi2(a: int, b: int, c: int, d: int) -> float:
    """Compute p-value from 2x2 contingency table.

    Uses scipy's Fisher exact test when available, otherwise falls
    back to chi-squared approximation.

    Args:
        a: Overlap count.
        b: State-only count.
        c: Annotation-only count.
        d: Neither count.

    Returns:
        One-sided p-value for enrichment.
    """
    if HAS_SCIPY:
        table = [[a, b], [c, d]]
        _, pval = scipy_stats.fisher_exact(table, alternative="greater")
        return float(pval)

    # Chi-squared approximation fallback
    n = a + b + c + d
    if n == 0:
        return 1.0

    expected_a = (a + b) * (a + c) / n
    if expected_a <= 0:
        return 1.0

    chi2 = (a - expected_a) ** 2 / expected_a
    # Approximate p-value from chi2 with 1 df using complementary error function
    pval = 0.5 * math.erfc(math.sqrt(chi2 / 2.0))
    return max(0.0, min(1.0, pval))


def segment_genome(
    state_assignments: list[int],
    bin_size: int = 200,
    min_segment: int = 1,
) -> list[dict]:
    """Convert bin-level state assignments to contiguous genomic segments.

    Scans state assignments and groups consecutive bins with the same
    state into segments with genomic coordinates.

    Args:
        state_assignments: Per-bin state assignments.
        bin_size: Size of each genomic bin in base pairs.
        min_segment: Minimum number of consecutive bins required to
            form a segment.

    Returns:
        List of segment dictionaries with keys: ``state``, ``start_bin``,
        ``end_bin``, ``start_bp``, ``end_bp``, ``length_bins``,
        ``length_bp``.

    Raises:
        ValueError: If state_assignments is empty.
    """
    if not state_assignments:
        raise ValueError("state_assignments must not be empty")

    logger.info(
        "Segmenting %d bins (bin_size=%d bp, min_segment=%d)",
        len(state_assignments),
        bin_size,
        min_segment,
    )

    segments: list[dict] = []
    current_state = state_assignments[0]
    seg_start = 0

    for i in range(1, len(state_assignments)):
        if state_assignments[i] != current_state:
            length_bins = i - seg_start
            if length_bins >= min_segment:
                segments.append(
                    {
                        "state": current_state,
                        "start_bin": seg_start,
                        "end_bin": i,
                        "start_bp": seg_start * bin_size,
                        "end_bp": i * bin_size,
                        "length_bins": length_bins,
                        "length_bp": length_bins * bin_size,
                    }
                )
            current_state = state_assignments[i]
            seg_start = i

    # Final segment
    length_bins = len(state_assignments) - seg_start
    if length_bins >= min_segment:
        segments.append(
            {
                "state": current_state,
                "start_bin": seg_start,
                "end_bin": len(state_assignments),
                "start_bp": seg_start * bin_size,
                "end_bp": len(state_assignments) * bin_size,
                "length_bins": length_bins,
                "length_bp": length_bins * bin_size,
            }
        )

    logger.info(
        "Generated %d segments covering %d bp",
        len(segments),
        sum(s["length_bp"] for s in segments),
    )
    return segments


def compare_chromatin_states(
    states_a: dict,
    states_b: dict,
) -> dict:
    """Compare chromatin state maps between two conditions or cell types.

    Analyses state switching between conditions and identifies
    differentially segmented regions.

    Args:
        states_a: Model dictionary from condition A (as returned by
            :func:`learn_chromatin_states`).
        states_b: Model dictionary from condition B.

    Returns:
        Dictionary with keys:

        - ``n_bins``: Number of bins compared.
        - ``concordance``: Fraction of bins with same state.
        - ``state_switches``: Counts of state transitions between conditions.
        - ``switch_matrix``: Full transition count matrix (A -> B).
        - ``differential_segments``: List of regions that changed state.
        - ``summary``: High-level statistics.

    Raises:
        ValueError: If assignment lists have different lengths.
    """
    assign_a = states_a["assignments"]
    assign_b = states_b["assignments"]

    if len(assign_a) != len(assign_b):
        raise ValueError(f"Assignment lengths must match: A={len(assign_a)}, B={len(assign_b)}")

    n_bins = len(assign_a)
    n_states_a = states_a["states"]
    n_states_b = states_b["states"]
    max_states = max(n_states_a, n_states_b)

    logger.info(
        "Comparing chromatin states: %d bins, A=%d states, B=%d states",
        n_bins,
        n_states_a,
        n_states_b,
    )

    # Count concordant bins
    concordant = sum(1 for a, b in zip(assign_a, assign_b) if a == b)
    concordance = concordant / n_bins if n_bins > 0 else 0.0

    # Build switch matrix
    switch_matrix: list[list[int]] = [[0] * max_states for _ in range(max_states)]
    for a, b in zip(assign_a, assign_b):
        if a < max_states and b < max_states:
            switch_matrix[a][b] += 1

    # Summarise state switches
    state_switches: dict[str, int] = {}
    for a_state in range(max_states):
        for b_state in range(max_states):
            if a_state != b_state and switch_matrix[a_state][b_state] > 0:
                key = f"{a_state}->{b_state}"
                state_switches[key] = switch_matrix[a_state][b_state]

    # Sort switches by frequency
    state_switches = dict(sorted(state_switches.items(), key=lambda x: x[1], reverse=True))

    # Identify differential segments (contiguous regions that changed)
    differential_segments: list[dict] = []
    in_diff = False
    diff_start = 0

    for i in range(n_bins):
        if assign_a[i] != assign_b[i]:
            if not in_diff:
                diff_start = i
                in_diff = True
        else:
            if in_diff:
                differential_segments.append(
                    {
                        "start_bin": diff_start,
                        "end_bin": i,
                        "length_bins": i - diff_start,
                        "states_a": list(set(assign_a[diff_start:i])),
                        "states_b": list(set(assign_b[diff_start:i])),
                    }
                )
                in_diff = False

    if in_diff:
        differential_segments.append(
            {
                "start_bin": diff_start,
                "end_bin": n_bins,
                "length_bins": n_bins - diff_start,
                "states_a": list(set(assign_a[diff_start:n_bins])),
                "states_b": list(set(assign_b[diff_start:n_bins])),
            }
        )

    result = {
        "n_bins": n_bins,
        "concordance": round(concordance, 6),
        "discordance": round(1.0 - concordance, 6),
        "concordant_bins": concordant,
        "discordant_bins": n_bins - concordant,
        "state_switches": state_switches,
        "switch_matrix": switch_matrix,
        "differential_segments": differential_segments,
        "n_differential_segments": len(differential_segments),
        "summary": {
            "n_states_a": n_states_a,
            "n_states_b": n_states_b,
            "concordance_pct": round(concordance * 100, 2),
            "n_switch_types": len(state_switches),
            "top_switch": (list(state_switches.keys())[0] if state_switches else None),
            "top_switch_count": (list(state_switches.values())[0] if state_switches else 0),
        },
    }

    logger.info(
        "Comparison: %.1f%% concordance, %d switch types, " "%d differential segments",
        concordance * 100,
        len(state_switches),
        len(differential_segments),
    )
    return result


def compute_state_transition_rates(
    assignments: list[int],
    n_states: int,
) -> dict:
    """Compute and analyse state transition probability matrix.

    Counts transitions between consecutive bins and normalises to
    obtain transition probabilities.  Also computes summary statistics
    such as self-transition rates (state stability) and entropy of
    each state's transition distribution.

    Args:
        assignments: Per-bin state assignments.
        n_states: Total number of chromatin states.

    Returns:
        Dictionary with keys:

        - ``transition_counts``: Raw count matrix (list[list[int]]).
        - ``transition_probabilities``: Row-normalised probability
          matrix (list[list[float]]).
        - ``self_transition_rates``: Per-state probability of
          staying in the same state.
        - ``transition_entropy``: Per-state entropy of the transition
          distribution (higher = more variable transitions).
        - ``mean_segment_length``: Per-state average segment length
          derived from self-transition rate.
        - ``n_transitions``: Total number of transitions.

    Raises:
        ValueError: If assignments is empty or n_states < 1.
    """
    if not assignments:
        raise ValueError("assignments must not be empty")
    if n_states < 1:
        raise ValueError(f"n_states must be >= 1, got {n_states}")

    logger.info(
        "Computing transition rates: %d bins, %d states",
        len(assignments),
        n_states,
    )

    # Count transitions
    counts: list[list[int]] = [[0] * n_states for _ in range(n_states)]
    for i in range(len(assignments) - 1):
        s_from = assignments[i]
        s_to = assignments[i + 1]
        if 0 <= s_from < n_states and 0 <= s_to < n_states:
            counts[s_from][s_to] += 1

    # Normalise to probabilities
    probabilities: list[list[float]] = []
    for k in range(n_states):
        row_sum = sum(counts[k])
        if row_sum > 0:
            probabilities.append([c / row_sum for c in counts[k]])
        else:
            # Uniform if no transitions observed from this state
            probabilities.append([1.0 / n_states] * n_states)

    # Self-transition rates
    self_rates: list[float] = []
    for k in range(n_states):
        self_rates.append(round(probabilities[k][k], 6))

    # Transition entropy per state
    entropies: list[float] = []
    for k in range(n_states):
        h = 0.0
        for p in probabilities[k]:
            if p > 0:
                h -= p * math.log2(p)
        entropies.append(round(h, 6))

    # Mean segment length: E[L] = 1 / (1 - p_self)
    mean_lengths: list[float] = []
    for k in range(n_states):
        p_self = probabilities[k][k]
        if p_self < 1.0:
            mean_lengths.append(round(1.0 / (1.0 - p_self), 4))
        else:
            mean_lengths.append(float("inf"))

    n_transitions = len(assignments) - 1

    result = {
        "transition_counts": counts,
        "transition_probabilities": [[round(p, 6) for p in row] for row in probabilities],
        "self_transition_rates": self_rates,
        "transition_entropy": entropies,
        "mean_segment_length": mean_lengths,
        "n_transitions": n_transitions,
        "n_states": n_states,
    }

    logger.info(
        "Computed transition rates: mean self-rate=%.4f, " "mean entropy=%.4f",
        sum(self_rates) / n_states if n_states > 0 else 0,
        sum(entropies) / n_states if n_states > 0 else 0,
    )
    return result
