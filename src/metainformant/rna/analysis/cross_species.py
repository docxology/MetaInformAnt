"""Cross-species gene expression comparison module.

Provides tools for comparing gene expression patterns across species
using ortholog mappings. Supports multi-species comparative studies
including expression conservation scoring, divergence analysis,
phylogenetic expression profiling, and cross-species PCA.

All implementations use numpy, scipy, and pandas. No mocking, no
placeholder data.

Main Functions:
    Ortholog Mapping:
        - build_ortholog_map: Build gene ortholog mapping from DataFrame
        - map_expression_to_orthologs: Map expression matrix to ortholog space

    Expression Comparison:
        - compare_expression_across_species: Multi-species expression comparison
        - compute_expression_conservation: Gene-level conservation scoring
        - identify_divergent_genes: Find divergently expressed genes

    Divergence & Phylogeny:
        - compute_expression_divergence_matrix: Pairwise species divergence
        - phylogenetic_expression_profile: Expression along phylogenetic tree

    Dimensionality Reduction:
        - cross_species_pca: PCA across species in shared ortholog space

Example:
    >>> from metainformant.rna.analysis import cross_species
    >>> import pandas as pd
    >>> # Build ortholog mapping
    >>> orthologs = pd.DataFrame({
    ...     "human_gene": ["TP53", "BRCA1", "EGFR"],
    ...     "mouse_gene": ["Trp53", "Brca1", "Egfr"],
    ... })
    >>> orth_map = cross_species.build_ortholog_map(orthologs, "human_gene", "mouse_gene")
    >>> # Map expression to ortholog space
    >>> human_expr = pd.DataFrame(...)
    >>> mapped = cross_species.map_expression_to_orthologs(human_expr, orth_map)
"""

from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


# =============================================================================
# Type Definitions
# =============================================================================

AggregationMethod = Literal["mean", "max", "sum"]
CorrelationMethod = Literal["spearman", "pearson", "euclidean"]


# =============================================================================
# Ortholog Mapping
# =============================================================================


def build_ortholog_map(
    orthologs_df: pd.DataFrame,
    source_col: str,
    target_col: str,
) -> Dict[str, List[str]]:
    """Build a gene ortholog mapping from a DataFrame of ortholog relationships.

    Parses an ortholog relationship table and constructs a dictionary mapping
    each source gene to its corresponding target gene(s). Handles one-to-one,
    one-to-many, and many-to-many ortholog relationships.

    Args:
        orthologs_df: DataFrame containing ortholog relationships. Must contain
            at least the two columns specified by source_col and target_col.
            Each row represents one ortholog pair. Duplicate pairs are deduplicated.
        source_col: Name of the column containing source species gene IDs.
        target_col: Name of the column containing target species gene IDs.

    Returns:
        Dictionary mapping each source gene ID (str) to a list of target gene
        IDs (List[str]). One-to-one orthologs produce single-element lists.
        One-to-many orthologs produce multi-element lists.

    Raises:
        ValueError: If source_col or target_col is not found in orthologs_df,
            or if orthologs_df is empty.

    Example:
        >>> orthologs = pd.DataFrame({
        ...     "human": ["TP53", "BRCA1", "EGFR", "EGFR"],
        ...     "mouse": ["Trp53", "Brca1", "Egfr", "Erbb1"],
        ... })
        >>> orth_map = build_ortholog_map(orthologs, "human", "mouse")
        >>> orth_map["TP53"]
        ['Trp53']
        >>> orth_map["EGFR"]
        ['Egfr', 'Erbb1']
    """
    if orthologs_df.empty:
        raise ValueError("orthologs_df cannot be empty")

    if source_col not in orthologs_df.columns:
        raise ValueError(
            f"Source column '{source_col}' not found in orthologs_df. "
            f"Available columns: {orthologs_df.columns.tolist()}"
        )

    if target_col not in orthologs_df.columns:
        raise ValueError(
            f"Target column '{target_col}' not found in orthologs_df. "
            f"Available columns: {orthologs_df.columns.tolist()}"
        )

    # Drop rows with missing values in the key columns
    valid = orthologs_df[[source_col, target_col]].dropna()

    if valid.empty:
        logger.warning("No valid ortholog pairs after removing NaN values")
        return {}

    # Deduplicate pairs
    valid = valid.drop_duplicates(subset=[source_col, target_col])

    # Build mapping: source gene -> list of unique target genes
    ortholog_map: Dict[str, List[str]] = {}
    for _, row in valid.iterrows():
        src = str(row[source_col])
        tgt = str(row[target_col])
        if src not in ortholog_map:
            ortholog_map[src] = []
        if tgt not in ortholog_map[src]:
            ortholog_map[src].append(tgt)

    n_one_to_one = sum(1 for v in ortholog_map.values() if len(v) == 1)
    n_one_to_many = sum(1 for v in ortholog_map.values() if len(v) > 1)
    logger.info(
        f"Built ortholog map: {len(ortholog_map)} source genes "
        f"({n_one_to_one} one-to-one, {n_one_to_many} one-to-many)"
    )

    return ortholog_map


def map_expression_to_orthologs(
    expression_df: pd.DataFrame,
    ortholog_map: Dict[str, List[str]],
    aggregation: AggregationMethod = "mean",
) -> pd.DataFrame:
    """Map an expression matrix to ortholog space.

    Transforms a gene expression matrix from one species' gene space into
    ortholog gene IDs. For one-to-many ortholog relationships, expression
    values are aggregated across the multiple source genes mapping to the
    same target using the specified method.

    Args:
        expression_df: Expression matrix with genes as rows and samples as
            columns. Row index should contain gene IDs matching the keys
            of ortholog_map.
        ortholog_map: Dictionary mapping source gene IDs to lists of target
            gene IDs (from build_ortholog_map).
        aggregation: Method for aggregating expression when multiple source
            genes map to the same target gene:
            - "mean": Average expression across source genes
            - "max": Maximum expression across source genes
            - "sum": Sum of expression across source genes

    Returns:
        Expression matrix with target ortholog gene IDs as rows and original
        samples as columns. Genes not found in ortholog_map are excluded.
        Target genes with multiple source mappings are aggregated.

    Raises:
        ValueError: If expression_df is empty or aggregation method is unknown.

    Example:
        >>> expr = pd.DataFrame({
        ...     "sample1": [10.0, 20.0, 30.0],
        ...     "sample2": [15.0, 25.0, 35.0],
        ... }, index=["GeneA", "GeneB", "GeneC"])
        >>> orth_map = {"GeneA": ["ortho1"], "GeneB": ["ortho1"], "GeneC": ["ortho2"]}
        >>> mapped = map_expression_to_orthologs(expr, orth_map, aggregation="mean")
        >>> # ortho1 is average of GeneA and GeneB; ortho2 is GeneC
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    if aggregation not in ("mean", "max", "sum"):
        raise ValueError(f"Unknown aggregation method: {aggregation}. Valid: mean, max, sum")

    # Build reverse map: target_gene -> list of source genes present in expression_df
    target_to_sources: Dict[str, List[str]] = {}
    genes_in_expr = set(expression_df.index)

    for source_gene, target_genes in ortholog_map.items():
        if source_gene not in genes_in_expr:
            continue
        for target_gene in target_genes:
            if target_gene not in target_to_sources:
                target_to_sources[target_gene] = []
            if source_gene not in target_to_sources[target_gene]:
                target_to_sources[target_gene].append(source_gene)

    if not target_to_sources:
        logger.warning("No genes in expression_df matched ortholog_map keys")
        return pd.DataFrame(columns=expression_df.columns)

    # Aggregate expression for each target gene
    result_rows: Dict[str, pd.Series] = {}
    for target_gene, source_genes in target_to_sources.items():
        source_expr = expression_df.loc[source_genes]

        if len(source_genes) == 1:
            result_rows[target_gene] = source_expr.iloc[0]
        else:
            if aggregation == "mean":
                result_rows[target_gene] = source_expr.mean(axis=0)
            elif aggregation == "max":
                result_rows[target_gene] = source_expr.max(axis=0)
            elif aggregation == "sum":
                result_rows[target_gene] = source_expr.sum(axis=0)

    mapped_df = pd.DataFrame(result_rows).T
    mapped_df.index.name = "gene"

    n_mapped = len(mapped_df)
    n_original = len(expression_df)
    logger.info(f"Mapped {n_original} source genes to {n_mapped} ortholog genes " f"(aggregation={aggregation})")

    return mapped_df


# =============================================================================
# Expression Conservation
# =============================================================================


def compute_expression_conservation(
    expr_a: pd.DataFrame,
    expr_b: pd.DataFrame,
    method: CorrelationMethod = "spearman",
) -> pd.DataFrame:
    """Compute gene-level expression conservation between two species.

    For each gene present in both expression matrices, computes a conservation
    score measuring how similar the expression pattern is across species.
    Uses correlation-based or distance-based measures across shared samples
    or conditions.

    Args:
        expr_a: Expression matrix for species A with genes as rows and
            samples/conditions as columns. Should be normalized.
        expr_b: Expression matrix for species B with genes as rows and
            samples/conditions as columns. Gene IDs (index) must overlap
            with expr_a (i.e., mapped to shared ortholog space). Column
            count must match expr_a for correlation methods.
        method: Conservation metric:
            - "spearman": Spearman rank correlation across samples
            - "pearson": Pearson correlation across samples
            - "euclidean": Negative normalized Euclidean distance (higher = more conserved)

    Returns:
        DataFrame with columns:
            - gene_id: Shared ortholog gene identifier
            - correlation: Conservation score (higher = more conserved).
              For correlation methods, range is [-1, 1]. For euclidean,
              range is [0, 1] where 1 means identical.
            - p_value: Statistical significance (for correlation methods).
              For euclidean, this is derived from a permutation-based z-score.
            - conserved: Boolean flag (True if correlation > 0.5 and p < 0.05)

    Raises:
        ValueError: If no shared genes exist between expr_a and expr_b,
            or if method is unknown.

    Example:
        >>> human_expr = pd.DataFrame({"cond1": [10, 20], "cond2": [15, 25]},
        ...                           index=["ortho1", "ortho2"])
        >>> mouse_expr = pd.DataFrame({"cond1": [12, 18], "cond2": [14, 22]},
        ...                           index=["ortho1", "ortho2"])
        >>> conservation = compute_expression_conservation(human_expr, mouse_expr)
    """
    if method not in ("spearman", "pearson", "euclidean"):
        raise ValueError(f"Unknown method: {method}. Valid: spearman, pearson, euclidean")

    # Find shared genes
    shared_genes = sorted(set(expr_a.index) & set(expr_b.index))
    if not shared_genes:
        raise ValueError(
            "No shared genes between expr_a and expr_b. " "Ensure both are mapped to the same ortholog space."
        )

    logger.info(f"Computing expression conservation for {len(shared_genes)} shared genes " f"using {method} method")

    results = []

    for gene in shared_genes:
        vals_a = expr_a.loc[gene].values.astype(float)
        vals_b = expr_b.loc[gene].values.astype(float)

        # Ensure matching dimensions by using the minimum length
        min_len = min(len(vals_a), len(vals_b))
        vals_a = vals_a[:min_len]
        vals_b = vals_b[:min_len]

        if min_len < 2:
            # Cannot compute correlation with fewer than 2 data points
            results.append(
                {
                    "gene_id": gene,
                    "correlation": np.nan,
                    "p_value": np.nan,
                    "conserved": False,
                }
            )
            continue

        if method == "spearman":
            corr, pval = stats.spearmanr(vals_a, vals_b)
        elif method == "pearson":
            corr, pval = stats.pearsonr(vals_a, vals_b)
        elif method == "euclidean":
            # Normalize both vectors to unit scale before distance
            range_a = vals_a.max() - vals_a.min()
            range_b = vals_b.max() - vals_b.min()
            norm_a = (vals_a - vals_a.min()) / range_a if range_a > 0 else np.zeros_like(vals_a)
            norm_b = (vals_b - vals_b.min()) / range_b if range_b > 0 else np.zeros_like(vals_b)

            eucl_dist = np.sqrt(np.mean((norm_a - norm_b) ** 2))
            # Convert distance to similarity score in [0, 1]
            # Maximum possible RMSE for [0,1]-normalized vectors is 1.0
            corr = float(1.0 - eucl_dist)

            # Approximate p-value using permutation z-score
            n_perm = 200
            rng = np.random.default_rng(42)
            perm_dists = np.empty(n_perm)
            for p_idx in range(n_perm):
                perm_b = rng.permutation(norm_b)
                perm_dists[p_idx] = np.sqrt(np.mean((norm_a - perm_b) ** 2))
            perm_similarities = 1.0 - perm_dists
            perm_mean = perm_similarities.mean()
            perm_std = perm_similarities.std()
            if perm_std > 0:
                z_score = (corr - perm_mean) / perm_std
                pval = float(stats.norm.sf(z_score))
            else:
                pval = 1.0 if corr <= perm_mean else 0.0

        # Handle NaN from correlation functions
        if np.isnan(corr):
            corr = 0.0
        if np.isnan(pval):
            pval = 1.0

        conserved = bool(corr > 0.5 and pval < 0.05)

        results.append(
            {
                "gene_id": gene,
                "correlation": float(corr),
                "p_value": float(pval),
                "conserved": conserved,
            }
        )

    result_df = pd.DataFrame(results)

    n_conserved = result_df["conserved"].sum()
    logger.info(f"Expression conservation: {n_conserved}/{len(result_df)} genes conserved " f"(corr > 0.5, p < 0.05)")

    return result_df


def identify_divergent_genes(
    conservation_df: pd.DataFrame,
    threshold: float = 0.3,
) -> pd.DataFrame:
    """Identify genes with divergent expression patterns across species.

    Filters the output of compute_expression_conservation to find genes
    whose expression conservation falls below a specified threshold,
    indicating evolutionary divergence in regulation.

    Args:
        conservation_df: DataFrame from compute_expression_conservation with
            columns: gene_id, correlation, p_value, conserved.
        threshold: Maximum correlation score for a gene to be considered
            divergent. Genes with correlation <= threshold are returned.
            Default is 0.3, indicating weak or no conservation.

    Returns:
        DataFrame with same columns as input, filtered to genes with
        correlation <= threshold. Sorted by correlation ascending (most
        divergent first). Returns empty DataFrame if no genes meet criteria.

    Raises:
        ValueError: If conservation_df is missing required columns.

    Example:
        >>> conservation = compute_expression_conservation(human_expr, mouse_expr)
        >>> divergent = identify_divergent_genes(conservation, threshold=0.3)
        >>> print(f"Found {len(divergent)} divergently expressed genes")
    """
    required_cols = {"gene_id", "correlation"}
    missing = required_cols - set(conservation_df.columns)
    if missing:
        raise ValueError(f"conservation_df missing required columns: {missing}")

    if conservation_df.empty:
        logger.info("Empty conservation_df, returning empty DataFrame")
        return conservation_df.copy()

    # Filter genes below threshold, handling NaN correlations as divergent
    mask = conservation_df["correlation"].fillna(-1.0) <= threshold
    divergent = conservation_df.loc[mask].copy()

    # Sort by correlation ascending (most divergent first)
    divergent = divergent.sort_values("correlation", ascending=True).reset_index(drop=True)

    logger.info(
        f"Identified {len(divergent)} divergent genes "
        f"(correlation <= {threshold}) out of {len(conservation_df)} total"
    )

    return divergent


# =============================================================================
# Multi-Species Comparison
# =============================================================================


def compare_expression_across_species(
    species_expressions: Dict[str, pd.DataFrame],
    ortholog_maps: Dict[Tuple[str, str], Dict[str, List[str]]],
) -> pd.DataFrame:
    """Compare expression patterns across multiple species using shared orthologs.

    For each pair of species, maps expression to a shared ortholog space and
    computes per-gene conservation scores. Aggregates across all pairwise
    comparisons to produce a single conservation metric per gene.

    Args:
        species_expressions: Dictionary mapping species name (str) to expression
            DataFrame (genes x samples). Each DataFrame should be normalized.
        ortholog_maps: Dictionary mapping species pairs (source, target) to
            ortholog dictionaries (from build_ortholog_map). Keys are tuples
            of (species_a, species_b) where species_a genes map to species_b.

    Returns:
        DataFrame with columns:
            - gene_id: Ortholog gene identifier (from the first species listed)
            - mean_conservation: Mean conservation score across all species pairs
            - min_conservation: Minimum conservation score across pairs
            - max_conservation: Maximum conservation score across pairs
            - n_species_compared: Number of pairwise comparisons for this gene
            - conserved_in_all: True if gene is conserved in all pairwise comparisons

    Raises:
        ValueError: If fewer than 2 species provided or no ortholog maps match.

    Example:
        >>> species_expr = {
        ...     "human": human_df,
        ...     "mouse": mouse_df,
        ...     "zebrafish": zebrafish_df,
        ... }
        >>> orth_maps = {
        ...     ("human", "mouse"): h2m_map,
        ...     ("human", "zebrafish"): h2z_map,
        ... }
        >>> comparison = compare_expression_across_species(species_expr, orth_maps)
    """
    species_names = list(species_expressions.keys())
    if len(species_names) < 2:
        raise ValueError(f"Need at least 2 species for comparison, got {len(species_names)}")

    logger.info(f"Comparing expression across {len(species_names)} species: {species_names}")

    # Collect pairwise conservation results
    gene_scores: Dict[str, List[float]] = {}

    for (sp_a, sp_b), orth_map in ortholog_maps.items():
        if sp_a not in species_expressions or sp_b not in species_expressions:
            logger.warning(f"Skipping pair ({sp_a}, {sp_b}): species not in species_expressions")
            continue

        expr_a = species_expressions[sp_a]
        expr_b = species_expressions[sp_b]

        # Map species A expression to ortholog space
        mapped_a = map_expression_to_orthologs(expr_a, orth_map, aggregation="mean")
        if mapped_a.empty:
            logger.warning(f"No genes mapped for ({sp_a}, {sp_b}), skipping")
            continue

        # Find shared genes in ortholog space
        shared_genes = sorted(set(mapped_a.index) & set(expr_b.index))
        if not shared_genes:
            logger.warning(f"No shared ortholog genes between {sp_a} and {sp_b}")
            continue

        # Compute conservation for shared genes
        conservation = compute_expression_conservation(
            mapped_a.loc[shared_genes],
            expr_b.loc[shared_genes],
            method="spearman",
        )

        for _, row in conservation.iterrows():
            gene = row["gene_id"]
            corr = row["correlation"]
            if not np.isnan(corr):
                if gene not in gene_scores:
                    gene_scores[gene] = []
                gene_scores[gene].append(corr)

    if not gene_scores:
        logger.warning("No pairwise comparisons produced results")
        return pd.DataFrame(
            columns=[
                "gene_id",
                "mean_conservation",
                "min_conservation",
                "max_conservation",
                "n_species_compared",
                "conserved_in_all",
            ]
        )

    # Aggregate across pairwise comparisons
    results = []
    for gene, scores in gene_scores.items():
        scores_arr = np.array(scores)
        results.append(
            {
                "gene_id": gene,
                "mean_conservation": float(np.mean(scores_arr)),
                "min_conservation": float(np.min(scores_arr)),
                "max_conservation": float(np.max(scores_arr)),
                "n_species_compared": len(scores),
                "conserved_in_all": bool(np.all(scores_arr > 0.5)),
            }
        )

    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values("mean_conservation", ascending=False).reset_index(drop=True)

    n_conserved = result_df["conserved_in_all"].sum()
    logger.info(
        f"Cross-species comparison: {n_conserved}/{len(result_df)} genes " f"conserved across all pairwise comparisons"
    )

    return result_df


# =============================================================================
# Divergence Matrix
# =============================================================================


def compute_expression_divergence_matrix(
    species_expressions: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Compute pairwise expression divergence between all species.

    For each pair of species, computes an overall expression divergence
    score based on 1 minus the median Spearman correlation across shared
    genes. Returns a symmetric distance matrix suitable for hierarchical
    clustering or tree construction.

    Args:
        species_expressions: Dictionary mapping species name (str) to
            expression DataFrame (genes x samples). All DataFrames must
            share a common gene index (i.e., already in ortholog space).

    Returns:
        Symmetric DataFrame with species as both rows and columns,
        containing divergence scores in [0, 2] (1 - correlation).
        Diagonal is 0.0 (species compared to itself). Higher values
        indicate greater expression divergence.

    Raises:
        ValueError: If fewer than 2 species provided.

    Example:
        >>> # All expression matrices should be in shared ortholog space
        >>> species_expr = {"human": human_orth, "mouse": mouse_orth, "fly": fly_orth}
        >>> div_matrix = compute_expression_divergence_matrix(species_expr)
        >>> div_matrix.loc["human", "mouse"]  # divergence between human and mouse
    """
    species_names = sorted(species_expressions.keys())
    n_species = len(species_names)

    if n_species < 2:
        raise ValueError(f"Need at least 2 species, got {n_species}")

    logger.info(f"Computing expression divergence matrix for {n_species} species")

    divergence = np.zeros((n_species, n_species))

    for i in range(n_species):
        for j in range(i + 1, n_species):
            sp_a = species_names[i]
            sp_b = species_names[j]

            expr_a = species_expressions[sp_a]
            expr_b = species_expressions[sp_b]

            # Find shared genes
            shared_genes = sorted(set(expr_a.index) & set(expr_b.index))

            if len(shared_genes) < 2:
                # No shared genes: maximum divergence
                logger.warning(f"Fewer than 2 shared genes between {sp_a} and {sp_b}")
                divergence[i, j] = 2.0
                divergence[j, i] = 2.0
                continue

            # Compute per-gene Spearman correlations across shared samples
            correlations = []
            min_samples = min(expr_a.shape[1], expr_b.shape[1])

            if min_samples >= 2:
                # Sample-level correlation per gene
                for gene in shared_genes:
                    vals_a = expr_a.loc[gene].values[:min_samples].astype(float)
                    vals_b = expr_b.loc[gene].values[:min_samples].astype(float)
                    corr, _ = stats.spearmanr(vals_a, vals_b)
                    if not np.isnan(corr):
                        correlations.append(corr)
            else:
                # Single sample per species: use gene-level correlation
                vals_a = expr_a.loc[shared_genes].iloc[:, 0].values.astype(float)
                vals_b = expr_b.loc[shared_genes].iloc[:, 0].values.astype(float)
                corr, _ = stats.spearmanr(vals_a, vals_b)
                if not np.isnan(corr):
                    correlations.append(corr)

            if correlations:
                median_corr = float(np.median(correlations))
                div_score = 1.0 - median_corr
            else:
                div_score = 2.0  # Maximum divergence

            divergence[i, j] = div_score
            divergence[j, i] = div_score

    result_df = pd.DataFrame(divergence, index=species_names, columns=species_names)

    logger.info(
        f"Divergence matrix computed. "
        f"Range: [{divergence[np.triu_indices(n_species, k=1)].min():.3f}, "
        f"{divergence[np.triu_indices(n_species, k=1)].max():.3f}]"
    )

    return result_df


# =============================================================================
# Phylogenetic Expression Profiling
# =============================================================================


def phylogenetic_expression_profile(
    expression_df: pd.DataFrame,
    species_tree: Dict[str, Any],
) -> pd.DataFrame:
    """Profile expression patterns along a phylogenetic tree.

    Traverses a phylogenetic tree structure and computes summary expression
    statistics at each node, enabling identification of lineage-specific
    expression changes. Leaf nodes correspond to species with expression
    data; internal nodes aggregate descendant expression.

    Args:
        expression_df: Expression matrix with genes as rows and columns
            named to match leaf node names in species_tree. Values should
            be normalized expression (e.g., mean expression per species).
        species_tree: Dictionary representing a phylogenetic tree with keys:
            - "name" (str): Node name (species name for leaves)
            - "children" (List[Dict], optional): Child nodes
            - "distance" (float, optional): Branch length to parent

    Returns:
        DataFrame with columns:
            - gene_id: Gene identifier
            - node: Tree node name
            - mean_expression: Mean expression at this node (leaf value or
              average of children)
            - variance: Expression variance across descendant leaves
            - n_leaves: Number of descendant leaves with data
            - branch_distance: Branch distance from parent
            - expression_change: Absolute change from parent node mean
              (NaN for root)

    Raises:
        ValueError: If expression_df is empty or species_tree is malformed.

    Example:
        >>> tree = {
        ...     "name": "root",
        ...     "children": [
        ...         {"name": "human", "distance": 0.1},
        ...         {"name": "mouse", "distance": 0.15},
        ...     ],
        ...     "distance": 0.0,
        ... }
        >>> expr = pd.DataFrame({"human": [10, 20], "mouse": [12, 18]},
        ...                     index=["gene1", "gene2"])
        >>> profile = phylogenetic_expression_profile(expr, tree)
    """
    if expression_df.empty:
        raise ValueError("expression_df cannot be empty")

    if "name" not in species_tree:
        raise ValueError("species_tree must contain a 'name' key")

    available_species = set(expression_df.columns)
    genes = expression_df.index.tolist()

    logger.info(f"Profiling expression for {len(genes)} genes across phylogenetic tree")

    # Collect all node-level expression statistics via recursive traversal
    results: List[Dict[str, Any]] = []

    def _traverse(
        node: Dict[str, Any],
        parent_means: Optional[Dict[str, float]],
    ) -> Dict[str, np.ndarray]:
        """Recursively traverse tree, collecting per-gene leaf expression arrays.

        Returns dict mapping gene_id to array of expression values at descendant leaves.
        """
        node_name = node["name"]
        branch_dist = node.get("distance", 0.0)
        children = node.get("children", [])

        if not children:
            # Leaf node: use expression data if available
            if node_name in available_species:
                leaf_values: Dict[str, np.ndarray] = {}
                for gene in genes:
                    val = float(expression_df.loc[gene, node_name])
                    leaf_values[gene] = np.array([val])

                    parent_change = np.nan
                    if parent_means is not None and gene in parent_means:
                        parent_change = abs(val - parent_means[gene])

                    results.append(
                        {
                            "gene_id": gene,
                            "node": node_name,
                            "mean_expression": val,
                            "variance": 0.0,
                            "n_leaves": 1,
                            "branch_distance": branch_dist,
                            "expression_change": parent_change,
                        }
                    )
                return leaf_values
            else:
                logger.warning(f"Leaf node '{node_name}' not found in expression data")
                return {}
        else:
            # Internal node: aggregate from children
            child_results: List[Dict[str, np.ndarray]] = []
            # First compute this node's means for passing to children
            # We need a two-pass approach: collect child data first
            for child in children:
                child_data = _traverse(child, None)  # placeholder parent
                child_results.append(child_data)

            # Merge child leaf values
            merged: Dict[str, np.ndarray] = {}
            for child_data in child_results:
                for gene, vals in child_data.items():
                    if gene in merged:
                        merged[gene] = np.concatenate([merged[gene], vals])
                    else:
                        merged[gene] = vals.copy()

            # Compute node-level statistics
            current_means: Dict[str, float] = {}
            for gene in genes:
                if gene in merged:
                    all_vals = merged[gene]
                    mean_val = float(np.mean(all_vals))
                    var_val = float(np.var(all_vals, ddof=0)) if len(all_vals) > 1 else 0.0
                    n_leaves = len(all_vals)
                    current_means[gene] = mean_val

                    parent_change = np.nan
                    if parent_means is not None and gene in parent_means:
                        parent_change = abs(mean_val - parent_means[gene])

                    results.append(
                        {
                            "gene_id": gene,
                            "node": node_name,
                            "mean_expression": mean_val,
                            "variance": var_val,
                            "n_leaves": n_leaves,
                            "branch_distance": branch_dist,
                            "expression_change": parent_change,
                        }
                    )

            # Now re-traverse children with correct parent means for expression_change
            # We already stored results, so update expression_change for child entries
            # For efficiency, update in-place
            for row in results:
                if row["node"] != node_name and np.isnan(row["expression_change"]):
                    gene = row["gene_id"]
                    if gene in current_means:
                        val = row["mean_expression"]
                        row["expression_change"] = abs(val - current_means[gene])

            return merged

    _traverse(species_tree, None)

    if not results:
        logger.warning("No results produced from tree traversal")
        return pd.DataFrame(
            columns=[
                "gene_id",
                "node",
                "mean_expression",
                "variance",
                "n_leaves",
                "branch_distance",
                "expression_change",
            ]
        )

    result_df = pd.DataFrame(results)

    n_nodes = result_df["node"].nunique()
    logger.info(f"Phylogenetic profile: {len(genes)} genes across {n_nodes} tree nodes")

    return result_df


# =============================================================================
# Cross-Species PCA
# =============================================================================


def cross_species_pca(
    species_expressions: Dict[str, pd.DataFrame],
    n_components: int = 2,
) -> Dict[str, Any]:
    """PCA across species using shared ortholog gene space.

    Combines expression data from multiple species into a single matrix
    (using genes present in all species) and performs principal component
    analysis. Each sample is labeled with its species of origin.

    Args:
        species_expressions: Dictionary mapping species name (str) to
            expression DataFrame (genes x samples). All DataFrames should
            share a common gene index (already in ortholog space).
        n_components: Number of principal components to compute.
            Will be reduced if exceeding max possible components.

    Returns:
        Dictionary with keys:
            - "coordinates": DataFrame of PC coordinates with columns
              PC1, PC2, ..., and an additional "species" column. Index
              is sample names (prefixed with species name if duplicates).
            - "explained_variance": Array of variance explained ratio per PC.
            - "loadings": DataFrame of gene loadings (genes x components).
            - "species_labels": List of species labels for each sample.
            - "shared_genes": List of gene IDs used in the analysis.

    Raises:
        ValueError: If fewer than 2 species provided or no shared genes.

    Example:
        >>> species_expr = {"human": human_orth_expr, "mouse": mouse_orth_expr}
        >>> pca_result = cross_species_pca(species_expr, n_components=3)
        >>> coords = pca_result["coordinates"]
        >>> # Plot PC1 vs PC2 colored by species
    """
    species_names = sorted(species_expressions.keys())
    n_species = len(species_names)

    if n_species < 2:
        raise ValueError(f"Need at least 2 species for cross-species PCA, got {n_species}")

    # Find genes shared across ALL species
    gene_sets = [set(species_expressions[sp].index) for sp in species_names]
    shared_genes = sorted(set.intersection(*gene_sets))

    if not shared_genes:
        raise ValueError(
            "No genes shared across all species. " "Ensure all expression matrices are in shared ortholog space."
        )

    logger.info(
        f"Cross-species PCA: {len(shared_genes)} shared genes, "
        f"{n_species} species, {n_components} components requested"
    )

    # Build combined expression matrix: shared genes x all samples
    combined_parts = []
    sample_labels = []
    species_labels = []

    for sp in species_names:
        expr = species_expressions[sp].loc[shared_genes]
        # Prefix sample names with species to avoid collisions
        renamed_cols = [f"{sp}:{col}" for col in expr.columns]
        expr_renamed = expr.copy()
        expr_renamed.columns = renamed_cols

        combined_parts.append(expr_renamed)
        sample_labels.extend(renamed_cols)
        species_labels.extend([sp] * len(renamed_cols))

    combined_df = pd.concat(combined_parts, axis=1)

    # Transpose: samples as rows, genes as columns
    X = combined_df.T.values.astype(float)
    n_samples, n_features = X.shape

    # Handle NaN
    if np.isnan(X).any():
        logger.warning("NaN values detected, imputing with column means")
        col_means = np.nanmean(X, axis=0)
        nan_idx = np.where(np.isnan(X))
        X[nan_idx] = np.take(col_means, nan_idx[1])

    # Constrain n_components
    max_components = min(n_samples, n_features)
    if n_components > max_components:
        logger.warning(f"Reducing n_components from {n_components} to {max_components}")
        n_components = max_components

    # Center and scale
    X_mean = X.mean(axis=0)
    X_centered = X - X_mean
    X_std = X.std(axis=0)
    X_std[X_std == 0] = 1.0
    X_scaled = X_centered / X_std

    # SVD
    try:
        U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
    except np.linalg.LinAlgError:
        logger.error("SVD did not converge for cross-species PCA")
        return {
            "coordinates": pd.DataFrame(),
            "explained_variance": np.array([]),
            "loadings": pd.DataFrame(),
            "species_labels": species_labels,
            "shared_genes": shared_genes,
        }

    # Select components
    U_k = U[:, :n_components]
    S_k = S[:n_components]
    Vt_k = Vt[:n_components, :]

    # PC scores
    transformed = U_k * S_k

    # Explained variance
    total_variance = (X_scaled**2).sum()
    explained_var = S_k**2 / (n_samples - 1)
    explained_var_ratio = explained_var / (total_variance / (n_samples - 1))

    # Gene loadings
    loadings = Vt_k.T * S_k / np.sqrt(n_samples - 1)

    # Build result DataFrames
    pc_names = [f"PC{i + 1}" for i in range(n_components)]

    coords_df = pd.DataFrame(transformed, index=sample_labels, columns=pc_names)
    coords_df["species"] = species_labels

    loadings_df = pd.DataFrame(loadings, index=shared_genes, columns=pc_names)

    logger.info(f"Cross-species PCA complete. " f"Variance explained: {[f'{v:.3f}' for v in explained_var_ratio]}")

    return {
        "coordinates": coords_df,
        "explained_variance": explained_var_ratio,
        "loadings": loadings_df,
        "species_labels": species_labels,
        "shared_genes": shared_genes,
    }
