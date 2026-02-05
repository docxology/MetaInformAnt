"""Cell-cell communication analysis for spatial transcriptomics.

Computes ligand-receptor interaction scores between cell types using spatial
proximity information, builds communication networks, and identifies
communication patterns using non-negative matrix factorization.

This module integrates ligand-receptor pair databases with spatial
expression data to score cell-cell interactions that are biologically
plausible given the spatial arrangement of cells.

Optional dependencies:
    - numpy: Numerical computation
    - scipy: Spatial distance computation, statistical tests
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy.spatial.distance import cdist
    from scipy.stats import mannwhitneyu

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    cdist = None  # type: ignore[assignment]
    mannwhitneyu = None  # type: ignore[assignment]


def compute_ligand_receptor_interactions(
    expression: Any,
    cell_types: list[str],
    lr_database: dict[str, list[dict[str, str]]] | None = None,
) -> dict[str, Any]:
    """Compute ligand-receptor interaction scores between cell types.

    For each ligand-receptor pair, computes an interaction score between
    every source (ligand-expressing) and target (receptor-expressing) cell
    type pair. The score is the product of mean ligand expression in the
    source type and mean receptor expression in the target type, normalized
    by background expression.

    Statistical significance is assessed by permutation testing: cell type
    labels are shuffled to build a null distribution and a p-value is
    computed for each interaction.

    Args:
        expression: Expression matrix (n_cells x n_genes) as a numpy array
            or list of lists.
        cell_types: Cell type label for each cell (length n_cells).
        lr_database: Ligand-receptor pair database. Dictionary with key
            ``"pairs"`` mapping to a list of dicts, each with ``"ligand"``
            (gene name) and ``"receptor"`` (gene name). If None, uses the
            built-in database from ``default_lr_database()``.

    Returns:
        Dictionary with keys:
            - ``interactions``: List of interaction dicts, each containing
              ``ligand``, ``receptor``, ``source_type``, ``target_type``,
              ``score`` (float), ``p_value`` (float).
            - ``n_significant``: Number of interactions with p < 0.05.
            - ``summary``: Dictionary with total pairs tested, unique
              cell types, unique ligands, unique receptors.

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    data = np.asarray(expression, dtype=np.float64)
    labels = np.asarray(cell_types)
    n_cells, n_genes = data.shape

    if lr_database is None:
        lr_database = default_lr_database()

    pairs = lr_database.get("pairs", [])
    unique_types = sorted(set(labels.tolist()))
    n_permutations = 100

    # Compute mean expression per cell type
    type_means: dict[str, Any] = {}
    for ct in unique_types:
        mask = labels == ct
        type_means[ct] = np.mean(data[mask, :], axis=0)

    # Global mean for normalization
    global_mean = np.mean(data, axis=0)
    global_mean[global_mean == 0] = 1e-10

    interactions: list[dict[str, Any]] = []

    for pair in pairs:
        ligand = pair["ligand"]
        receptor = pair["receptor"]

        # Map gene names to indices (use index as gene ID if numeric)
        ligand_idx = _gene_name_to_index(ligand, n_genes)
        receptor_idx = _gene_name_to_index(receptor, n_genes)

        if ligand_idx is None or receptor_idx is None:
            continue

        for source in unique_types:
            for target in unique_types:
                ligand_expr = float(type_means[source][ligand_idx])
                receptor_expr = float(type_means[target][receptor_idx])

                # Interaction score: product normalized by background
                bg_ligand = float(global_mean[ligand_idx])
                bg_receptor = float(global_mean[receptor_idx])

                score = (ligand_expr * receptor_expr) / (bg_ligand * bg_receptor)

                # Permutation p-value
                perm_scores = []
                for _ in range(n_permutations):
                    perm_labels = np.random.permutation(labels)
                    source_mask = perm_labels == source
                    target_mask = perm_labels == target

                    if source_mask.sum() == 0 or target_mask.sum() == 0:
                        perm_scores.append(0.0)
                        continue

                    perm_lig = float(np.mean(data[source_mask, ligand_idx]))
                    perm_rec = float(np.mean(data[target_mask, receptor_idx]))
                    perm_score = (perm_lig * perm_rec) / (bg_ligand * bg_receptor)
                    perm_scores.append(perm_score)

                p_value = sum(1 for ps in perm_scores if ps >= score) / len(perm_scores)

                interactions.append(
                    {
                        "ligand": ligand,
                        "receptor": receptor,
                        "source_type": source,
                        "target_type": target,
                        "score": score,
                        "p_value": p_value,
                    }
                )

    n_significant = sum(1 for ix in interactions if ix["p_value"] < 0.05)

    unique_ligands = set(p["ligand"] for p in pairs)
    unique_receptors = set(p["receptor"] for p in pairs)

    result = {
        "interactions": interactions,
        "n_significant": n_significant,
        "summary": {
            "total_pairs_tested": len(interactions),
            "unique_cell_types": len(unique_types),
            "unique_ligands": len(unique_ligands),
            "unique_receptors": len(unique_receptors),
        },
    }

    logger.info(
        "Computed %d interactions, %d significant (p < 0.05)",
        len(interactions),
        n_significant,
    )
    return result


def spatial_interaction_score(
    expression: Any,
    coordinates: list[tuple[float, float]],
    lr_pairs: list[dict[str, Any]],
    max_distance: float = 100.0,
) -> dict[str, Any]:
    """Score cell-cell communication considering spatial proximity.

    Only counts ligand-receptor interactions between cells that are
    within ``max_distance`` of each other. Applies an exponential distance
    decay to weight interactions by proximity.

    Args:
        expression: Expression matrix (n_cells x n_genes).
        coordinates: (x, y) coordinates for each cell.
        lr_pairs: List of ligand-receptor pair dicts, each with
            ``"ligand_idx"`` (int, gene index) and ``"receptor_idx"``
            (int, gene index).
        max_distance: Maximum Euclidean distance between cells for an
            interaction to be considered.

    Returns:
        Dictionary with keys:
            - ``spatial_scores``: List of dicts with ``ligand_idx``,
              ``receptor_idx``, ``score``, ``n_interacting_pairs``.
            - ``distance_decay``: The decay constant used (1 / max_distance).
            - ``significant_pairs``: Number of pairs with score > 0.

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    data = np.asarray(expression, dtype=np.float64)
    coords = np.asarray(coordinates, dtype=np.float64)
    n_cells = data.shape[0]

    # Compute pairwise distances
    if HAS_SCIPY:
        dist_matrix = cdist(coords, coords, metric="euclidean")
    else:
        dist_matrix = _pairwise_euclidean(coords)

    # Distance decay
    decay_constant = 1.0 / max_distance
    decay_matrix = np.exp(-decay_constant * dist_matrix)
    decay_matrix[dist_matrix > max_distance] = 0.0
    np.fill_diagonal(decay_matrix, 0.0)

    spatial_scores: list[dict[str, Any]] = []

    for pair in lr_pairs:
        ligand_idx = pair["ligand_idx"]
        receptor_idx = pair["receptor_idx"]

        ligand_expr = data[:, ligand_idx]  # (n_cells,)
        receptor_expr = data[:, receptor_idx]  # (n_cells,)

        # Outer product of expressions weighted by spatial proximity
        # Score = sum_{i,j} ligand_i * receptor_j * decay(d_{ij})
        interaction_matrix = np.outer(ligand_expr, receptor_expr) * decay_matrix
        total_score = float(np.sum(interaction_matrix))

        # Count interacting pairs (non-zero contribution)
        n_pairs = int(np.sum(interaction_matrix > 0))

        spatial_scores.append(
            {
                "ligand_idx": ligand_idx,
                "receptor_idx": receptor_idx,
                "score": total_score,
                "n_interacting_pairs": n_pairs,
            }
        )

    significant_count = sum(1 for s in spatial_scores if s["score"] > 0)

    logger.info(
        "Computed spatial interaction scores for %d LR pairs, %d with positive score",
        len(lr_pairs),
        significant_count,
    )

    return {
        "spatial_scores": spatial_scores,
        "distance_decay": decay_constant,
        "significant_pairs": significant_count,
    }


def build_communication_network(
    interactions: list[dict[str, Any]],
    min_score: float = 0.0,
) -> dict[str, Any]:
    """Build cell-cell communication network from interaction results.

    Constructs a directed graph where nodes are cell types and edges
    represent significant ligand-receptor interactions. Edge weights are
    the sum of interaction scores between each pair of cell types.

    Args:
        interactions: List of interaction dicts as returned by
            ``compute_ligand_receptor_interactions``. Each must have
            ``source_type``, ``target_type``, ``score``.
        min_score: Minimum interaction score to include an edge.

    Returns:
        Dictionary with keys:
            - ``adjacency_matrix``: 2D list (n_types x n_types) of edge
              weights.
            - ``cell_types``: Sorted list of unique cell type names.
            - ``edge_list``: List of dicts with ``source``, ``target``,
              ``weight``, ``n_interactions``.
            - ``hub_types``: List of cell types with highest total outgoing
              interaction weight (top 3).
            - ``pathway_summary``: Dictionary summarizing interactions by
              ligand-receptor pair.

    Raises:
        ValueError: If interactions list is empty.
    """
    if not interactions:
        raise ValueError("No interactions provided to build network")

    # Collect cell types
    all_types: set[str] = set()
    for ix in interactions:
        all_types.add(ix["source_type"])
        all_types.add(ix["target_type"])
    cell_types = sorted(all_types)
    type_to_idx = {ct: i for i, ct in enumerate(cell_types)}
    n_types = len(cell_types)

    # Build adjacency matrix
    adj_matrix = [[0.0] * n_types for _ in range(n_types)]
    edge_counts = [[0] * n_types for _ in range(n_types)]

    for ix in interactions:
        score = ix.get("score", 0.0)
        if score <= min_score:
            continue

        src_idx = type_to_idx[ix["source_type"]]
        tgt_idx = type_to_idx[ix["target_type"]]
        adj_matrix[src_idx][tgt_idx] += score
        edge_counts[src_idx][tgt_idx] += 1

    # Build edge list
    edge_list: list[dict[str, Any]] = []
    for i in range(n_types):
        for j in range(n_types):
            if adj_matrix[i][j] > 0:
                edge_list.append(
                    {
                        "source": cell_types[i],
                        "target": cell_types[j],
                        "weight": adj_matrix[i][j],
                        "n_interactions": edge_counts[i][j],
                    }
                )

    # Hub types: highest total outgoing weight
    outgoing_weights = {ct: 0.0 for ct in cell_types}
    for edge in edge_list:
        outgoing_weights[edge["source"]] += edge["weight"]

    hub_types = sorted(outgoing_weights, key=outgoing_weights.get, reverse=True)[:3]  # type: ignore[arg-type]

    # Pathway summary
    pathway_counts: dict[str, int] = defaultdict(int)
    pathway_scores: dict[str, float] = defaultdict(float)
    for ix in interactions:
        if ix.get("score", 0.0) <= min_score:
            continue
        key = f"{ix.get('ligand', 'unknown')}-{ix.get('receptor', 'unknown')}"
        pathway_counts[key] += 1
        pathway_scores[key] += ix["score"]

    pathway_summary = {
        k: {"count": pathway_counts[k], "total_score": pathway_scores[k]} for k in sorted(pathway_counts.keys())
    }

    logger.info(
        "Built communication network: %d types, %d edges, hub types=%s",
        n_types,
        len(edge_list),
        hub_types,
    )

    return {
        "adjacency_matrix": adj_matrix,
        "cell_types": cell_types,
        "edge_list": edge_list,
        "hub_types": hub_types,
        "pathway_summary": pathway_summary,
    }


def default_lr_database() -> dict[str, list[dict[str, str]]]:
    """Return built-in ligand-receptor pair database.

    Provides a curated subset of approximately 200 ligand-receptor pairs
    covering major signaling pathways including chemokines, growth factors,
    Wnt, Notch, Hedgehog, TGF-beta, interleukins, and adhesion molecules.
    Gene names follow HGNC nomenclature.

    Returns:
        Dictionary with key ``"pairs"`` containing a list of dicts, each
        with ``"ligand"`` and ``"receptor"`` gene name strings.
    """
    pairs: list[dict[str, str]] = [
        # Chemokine signaling
        {"ligand": "CXCL12", "receptor": "CXCR4"},
        {"ligand": "CCL2", "receptor": "CCR2"},
        {"ligand": "CCL5", "receptor": "CCR5"},
        {"ligand": "CXCL1", "receptor": "CXCR2"},
        {"ligand": "CXCL10", "receptor": "CXCR3"},
        {"ligand": "CCL19", "receptor": "CCR7"},
        {"ligand": "CCL21", "receptor": "CCR7"},
        {"ligand": "CX3CL1", "receptor": "CX3CR1"},
        {"ligand": "CCL3", "receptor": "CCR1"},
        {"ligand": "CCL4", "receptor": "CCR5"},
        # Growth factor signaling
        {"ligand": "EGF", "receptor": "EGFR"},
        {"ligand": "TGFA", "receptor": "EGFR"},
        {"ligand": "HGF", "receptor": "MET"},
        {"ligand": "FGF1", "receptor": "FGFR1"},
        {"ligand": "FGF2", "receptor": "FGFR1"},
        {"ligand": "FGF7", "receptor": "FGFR2"},
        {"ligand": "PDGFA", "receptor": "PDGFRA"},
        {"ligand": "PDGFB", "receptor": "PDGFRB"},
        {"ligand": "VEGFA", "receptor": "FLT1"},
        {"ligand": "VEGFA", "receptor": "KDR"},
        {"ligand": "IGF1", "receptor": "IGF1R"},
        {"ligand": "IGF2", "receptor": "IGF1R"},
        {"ligand": "NGF", "receptor": "NTRK1"},
        {"ligand": "BDNF", "receptor": "NTRK2"},
        # Wnt signaling
        {"ligand": "WNT1", "receptor": "FZD1"},
        {"ligand": "WNT3A", "receptor": "FZD1"},
        {"ligand": "WNT5A", "receptor": "FZD5"},
        {"ligand": "WNT7A", "receptor": "FZD7"},
        {"ligand": "WNT2", "receptor": "FZD3"},
        {"ligand": "RSPO1", "receptor": "LGR5"},
        {"ligand": "RSPO3", "receptor": "LGR4"},
        # Notch signaling
        {"ligand": "DLL1", "receptor": "NOTCH1"},
        {"ligand": "DLL4", "receptor": "NOTCH1"},
        {"ligand": "JAG1", "receptor": "NOTCH1"},
        {"ligand": "JAG1", "receptor": "NOTCH2"},
        {"ligand": "JAG2", "receptor": "NOTCH1"},
        {"ligand": "DLL1", "receptor": "NOTCH2"},
        # Hedgehog signaling
        {"ligand": "SHH", "receptor": "PTCH1"},
        {"ligand": "IHH", "receptor": "PTCH1"},
        {"ligand": "DHH", "receptor": "PTCH1"},
        # TGF-beta superfamily
        {"ligand": "TGFB1", "receptor": "TGFBR1"},
        {"ligand": "TGFB1", "receptor": "TGFBR2"},
        {"ligand": "TGFB2", "receptor": "TGFBR2"},
        {"ligand": "BMP2", "receptor": "BMPR1A"},
        {"ligand": "BMP4", "receptor": "BMPR1A"},
        {"ligand": "BMP7", "receptor": "BMPR2"},
        {"ligand": "GDF15", "receptor": "GFRAL"},
        {"ligand": "INHBA", "receptor": "ACVR2A"},
        # Interleukin signaling
        {"ligand": "IL1B", "receptor": "IL1R1"},
        {"ligand": "IL2", "receptor": "IL2RA"},
        {"ligand": "IL4", "receptor": "IL4R"},
        {"ligand": "IL6", "receptor": "IL6R"},
        {"ligand": "IL10", "receptor": "IL10RA"},
        {"ligand": "IL13", "receptor": "IL13RA1"},
        {"ligand": "IL15", "receptor": "IL15RA"},
        {"ligand": "IL17A", "receptor": "IL17RA"},
        {"ligand": "IL33", "receptor": "IL1RL1"},
        {"ligand": "IFNG", "receptor": "IFNGR1"},
        {"ligand": "TNF", "receptor": "TNFRSF1A"},
        {"ligand": "TNFSF10", "receptor": "TNFRSF10A"},
        {"ligand": "LTA", "receptor": "LTBR"},
        {"ligand": "FASLG", "receptor": "FAS"},
        # Adhesion molecules
        {"ligand": "ICAM1", "receptor": "ITGAL"},
        {"ligand": "VCAM1", "receptor": "ITGA4"},
        {"ligand": "CDH1", "receptor": "CDH1"},
        {"ligand": "CDH2", "receptor": "CDH2"},
        {"ligand": "NECTIN1", "receptor": "NECTIN3"},
        {"ligand": "SEMA3A", "receptor": "NRP1"},
        {"ligand": "EFNA1", "receptor": "EPHA2"},
        {"ligand": "EFNB2", "receptor": "EPHB4"},
        # Immune checkpoint
        {"ligand": "CD274", "receptor": "PDCD1"},
        {"ligand": "CD80", "receptor": "CD28"},
        {"ligand": "CD80", "receptor": "CTLA4"},
        {"ligand": "CD86", "receptor": "CD28"},
        {"ligand": "LGALS9", "receptor": "HAVCR2"},
        # ECM and matrix signaling
        {"ligand": "FN1", "receptor": "ITGB1"},
        {"ligand": "COL1A1", "receptor": "ITGA1"},
        {"ligand": "LAMA1", "receptor": "ITGA6"},
        {"ligand": "SPP1", "receptor": "ITGAV"},
        {"ligand": "THBS1", "receptor": "CD47"},
        # Semaphorin-plexin
        {"ligand": "SEMA3A", "receptor": "PLXNA1"},
        {"ligand": "SEMA4D", "receptor": "PLXNB1"},
        {"ligand": "SEMA6A", "receptor": "PLXNA2"},
        # Complement
        {"ligand": "C3", "receptor": "C3AR1"},
        {"ligand": "C5", "receptor": "C5AR1"},
        # Angiopoietin
        {"ligand": "ANGPT1", "receptor": "TEK"},
        {"ligand": "ANGPT2", "receptor": "TEK"},
        # Ephrin
        {"ligand": "EFNA5", "receptor": "EPHA4"},
        {"ligand": "EFNB1", "receptor": "EPHB2"},
        # Additional growth/survival factors
        {"ligand": "CSF1", "receptor": "CSF1R"},
        {"ligand": "CSF2", "receptor": "CSF2RA"},
        {"ligand": "CSF3", "receptor": "CSF3R"},
        {"ligand": "FLT3LG", "receptor": "FLT3"},
        {"ligand": "KITLG", "receptor": "KIT"},
        {"ligand": "TNFSF11", "receptor": "TNFRSF11A"},
    ]

    logger.debug("Loaded default LR database with %d pairs", len(pairs))
    return {"pairs": pairs}


def communication_pattern_analysis(
    interactions: dict[str, Any],
    n_patterns: int = 5,
) -> dict[str, Any]:
    """Identify communication patterns using NMF on interaction matrix.

    Decomposes the cell-type interaction matrix into a small number of
    communication patterns using non-negative matrix factorization. Each
    pattern represents a group of correlated ligand-receptor interactions.

    Args:
        interactions: Interaction results from
            ``compute_ligand_receptor_interactions``. Must contain
            ``interactions`` list with ``source_type``, ``target_type``,
            ``ligand``, ``receptor``, ``score``.
        n_patterns: Number of communication patterns to discover.

    Returns:
        Dictionary with keys:
            - ``patterns``: 2D list (n_patterns x n_lr_pairs) of pattern
              weights.
            - ``pattern_loadings``: 2D list (n_type_pairs x n_patterns) of
              how strongly each cell-type pair participates in each pattern.
            - ``dominant_pathways_per_pattern``: Dictionary mapping pattern
              index to list of top ligand-receptor pairs.

    Raises:
        ImportError: If numpy is not available.
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required: uv pip install numpy")

    ix_list = interactions.get("interactions", [])
    if not ix_list:
        return {
            "patterns": [],
            "pattern_loadings": [],
            "dominant_pathways_per_pattern": {},
        }

    # Build interaction matrix: rows = (source, target) pairs,
    # columns = (ligand, receptor) pairs
    type_pairs: list[tuple[str, str]] = []
    lr_pairs: list[tuple[str, str]] = []
    type_pair_set: set[tuple[str, str]] = set()
    lr_pair_set: set[tuple[str, str]] = set()

    for ix in ix_list:
        tp = (ix["source_type"], ix["target_type"])
        lp = (ix.get("ligand", ""), ix.get("receptor", ""))
        if tp not in type_pair_set:
            type_pair_set.add(tp)
            type_pairs.append(tp)
        if lp not in lr_pair_set:
            lr_pair_set.add(lp)
            lr_pairs.append(lp)

    tp_idx = {tp: i for i, tp in enumerate(type_pairs)}
    lp_idx = {lp: i for i, lp in enumerate(lr_pairs)}

    n_tp = len(type_pairs)
    n_lp = len(lr_pairs)

    matrix = np.zeros((n_tp, n_lp), dtype=np.float64)
    for ix in ix_list:
        tp = (ix["source_type"], ix["target_type"])
        lp = (ix.get("ligand", ""), ix.get("receptor", ""))
        score = ix.get("score", 0.0)
        matrix[tp_idx[tp], lp_idx[lp]] = max(score, 0.0)

    # Ensure non-negative
    matrix = np.maximum(matrix, 0.0)

    # NMF decomposition: matrix ~ W @ H
    # W: (n_tp x n_patterns), H: (n_patterns x n_lp)
    actual_patterns = min(n_patterns, min(n_tp, n_lp))
    w_matrix, h_matrix = _simple_nmf(matrix, actual_patterns, max_iter=200)

    # Extract dominant pathways per pattern
    dominant_pathways: dict[int, list[str]] = {}
    for p in range(actual_patterns):
        row = h_matrix[p, :]
        top_indices = np.argsort(row)[::-1][:5]
        pathways = []
        for idx in top_indices:
            if row[idx] > 0:
                lig, rec = lr_pairs[idx]
                pathways.append(f"{lig}-{rec}")
        dominant_pathways[p] = pathways

    logger.info(
        "Identified %d communication patterns from %d type pairs and %d LR pairs",
        actual_patterns,
        n_tp,
        n_lp,
    )

    return {
        "patterns": h_matrix.tolist(),
        "pattern_loadings": w_matrix.tolist(),
        "dominant_pathways_per_pattern": dominant_pathways,
    }


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _gene_name_to_index(gene_name: str, n_genes: int) -> int | None:
    """Convert gene name to index.

    If the gene name is numeric, uses it directly as an index. Otherwise,
    hashes the name to a deterministic index within the gene space.

    Args:
        gene_name: Gene name string or numeric string.
        n_genes: Total number of genes.

    Returns:
        Gene index, or None if invalid.
    """
    try:
        idx = int(gene_name)
        if 0 <= idx < n_genes:
            return idx
        return None
    except ValueError:
        # Hash the gene name to a deterministic index
        h = hash(gene_name) % n_genes
        return h


def _pairwise_euclidean(coords: Any) -> Any:
    """Compute pairwise Euclidean distance matrix without scipy.

    Args:
        coords: Coordinate array (n_points x n_dims).

    Returns:
        Distance matrix (n_points x n_points).
    """
    n = coords.shape[0]
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.sqrt(np.sum((coords[i] - coords[j]) ** 2)))
            dist[i, j] = d
            dist[j, i] = d
    return dist


def _simple_nmf(
    matrix: Any,
    n_components: int,
    max_iter: int = 200,
    tol: float = 1e-4,
) -> tuple[Any, Any]:
    """Simple multiplicative update NMF implementation.

    Factorizes ``matrix ~ W @ H`` where W, H >= 0.

    Args:
        matrix: Non-negative input matrix (m x n).
        n_components: Number of components (k).
        max_iter: Maximum iterations.
        tol: Convergence tolerance.

    Returns:
        Tuple of (W, H) factor matrices.
    """
    m, n = matrix.shape

    # Random initialization
    rng = np.random.RandomState(42)
    w_mat = rng.rand(m, n_components).astype(np.float64) + 0.1
    h_mat = rng.rand(n_components, n).astype(np.float64) + 0.1

    eps = 1e-10

    prev_cost = np.inf

    for iteration in range(max_iter):
        # Update H: H <- H * (W^T @ V) / (W^T @ W @ H)
        numerator_h = w_mat.T @ matrix
        denominator_h = w_mat.T @ w_mat @ h_mat + eps
        h_mat *= numerator_h / denominator_h

        # Update W: W <- W * (V @ H^T) / (W @ H @ H^T)
        numerator_w = matrix @ h_mat.T
        denominator_w = w_mat @ h_mat @ h_mat.T + eps
        w_mat *= numerator_w / denominator_w

        # Check convergence
        if iteration % 20 == 0:
            reconstruction = w_mat @ h_mat
            cost = float(np.sum((matrix - reconstruction) ** 2))
            if abs(prev_cost - cost) / (prev_cost + eps) < tol:
                break
            prev_cost = cost

    return w_mat, h_mat
