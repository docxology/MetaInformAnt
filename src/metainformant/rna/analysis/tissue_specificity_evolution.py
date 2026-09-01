"""Tissue-specificity evolution across a species phylogeny (descriptive-only).

Adapts the Mantica et al. 2024 workflow ("Evolution of tissue-specific
expression of ancestral genes across vertebrates and insects",
doi:10.1038/s41559-024-02398-5; code: github.com/fedemantica/bilaterian_GE)
to the MetaInformAnt hymenoptera panel:

1. Tau tissue-specificity computation (Yanai et al. 2005) with an optional
   minimum-mean-expression filter following Xu & Colgan 2025
   (doi:10.1093/molbev/msaf063: TPM normalization upstream, lowest-10%-mean-
   expression filter before tau).
2. Binary tissue-specificity state assignment per gene (orthogroup) per
   species per tissue using a tau cutoff (Mantica et al. default 0.75).
3. Descriptive Fitch-style parsimony gain/loss counting of the binary
   tissue-specificity state on a Newick species tree. This is a direct
   description of the minimum number of transitions consistent with the
   observed tip states; it is NOT a probabilistic reconstruction and infers
   no rates.
4. Duplication-coupling summary: orthology-class (single-copy vs multi-copy
   within the panel) versus tissue-specificity-gain status, following the
   Mantica et al. duplication-specialization coupling summary.

BOUNDARY: all outputs are descriptive statistics. No inferential
cross-species testing, no p-value claims, and no ancestral-probability
estimates are produced by this module. Any Wilcoxon/chi-square style
comparison (as in Xu & Colgan 2025) belongs to a post-evidence-freeze
gated analysis, not here.

Example:
    >>> import pandas as pd
    >>> from metainformant.rna.analysis import tissue_specificity_evolution as tse
    >>> expr = pd.DataFrame(
    ...     {"brain": [10.0, 1.0], "muscle": [1.0, 10.0]},
    ...     index=["og1", "og2"],
    ... )
    >>> tau = tse.compute_tau(expr)
    >>> tau.loc["og1", "tau"] > 0.9
    True
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from metainformant.core.utils import logging
from metainformant.dna.phylogeny.tree_analysis import from_newick

logger = logging.get_logger(__name__)

Tree = Dict[str, Any]

#: Tissue-specificity tau cutoff from Mantica et al. 2024 (Supplementary Fig. 2).
DEFAULT_TAU_CUTOFF = 0.75


# =============================================================================
# Tau computation
# =============================================================================


def compute_tau(
    expression_df: pd.DataFrame,
    min_mean_expression: Optional[float] = None,
) -> pd.DataFrame:
    """Compute tau tissue-specificity per gene across tissues (Yanai 2005).

    tau = sum_t (1 - x_t / max_t x_t) / (n_tissues - 1)

    Convention note (Xu & Colgan 2025): tau is computed on TPM-normalized
    values; because tau is scale-invariant per gene, log2 transformation
    does not change the result. Callers should pass TPM-scale expression and
    apply ``min_mean_expression`` (e.g. the 10th percentile of gene mean
    expression) to drop the lowest-expressed genes before tau, as in that
    paper's filtering step.

    Args:
        expression_df: Genes (rows) x tissues (columns) expression matrix.
        min_mean_expression: Optional minimum row-mean expression to retain a
            gene. ``None`` keeps all genes.

    Returns:
        DataFrame with a single ``tau`` column indexed by gene.

    Raises:
        ValueError: If fewer than two tissues are provided.
    """
    if expression_df.shape[1] < 2:
        raise ValueError("tau requires at least two tissues (columns)")

    matrix = expression_df.to_numpy(dtype=float)
    if np.any(matrix < 0):
        raise ValueError("expression values must be non-negative")

    retained = np.ones(matrix.shape[0], dtype=bool)
    if min_mean_expression is not None:
        retained = matrix.mean(axis=1) >= float(min_mean_expression)
        matrix = matrix[retained]

    row_max = matrix.max(axis=1)
    zero_max = row_max <= 0
    safe_max = np.where(zero_max, 1.0, row_max)
    proportions = matrix / safe_max[:, None]
    tau_values = (1.0 - proportions).sum(axis=1) / (matrix.shape[1] - 1)
    tau_values = np.where(zero_max, np.nan, tau_values)

    return pd.DataFrame({"tau": tau_values}, index=expression_df.index[retained])


def assign_tissue_specificity_states(
    expression_df: pd.DataFrame,
    tau_df: pd.DataFrame,
    tau_cutoff: float = DEFAULT_TAU_CUTOFF,
) -> pd.DataFrame:
    """Assign per-gene, per-species binary tissue-specificity states.

    A gene is tissue-specific (state = associated tissue name) when its tau
    meets the cutoff; otherwise the state is ``None``. The associated tissue
    is the argmax tissue of the expression row (Mantica et al.: "associated
    tissue" = highest-expression tissue for tau-passing genes).

    Args:
        expression_df: Genes x tissues expression (same genes as ``tau_df``).
        tau_df: Result of :func:`compute_tau` (or equivalent tau table).
        tau_cutoff: Tau value at or above which a gene is tissue-specific.

    Returns:
        DataFrame with columns ``tau`` and ``tissue_specific_state`` indexed
        by gene. ``tissue_specific_state`` is a tissue name or ``None``.
    """
    common = [gene for gene in tau_df.index if gene in expression_df.index]
    expression = expression_df.loc[common]
    taus = tau_df.loc[common, "tau"].to_numpy(dtype=float)

    states: List[Optional[str]] = []
    for idx in range(len(common)):
        tau_value = taus[idx]
        if np.isnan(tau_value) or tau_value < tau_cutoff:
            states.append(None)
        else:
            row = expression.iloc[idx]
            states.append(str(row.idxmax()))

    return pd.DataFrame(
        {"tau": taus, "tissue_specific_state": states},
        index=pd.Index(common, name=expression_df.index.name),
        dtype=object,
    )


# =============================================================================
# Parsimony gain/loss counting on a species tree
# =============================================================================


def _iter_tree_edges(tree: Tree, node: str) -> Iterable[Tuple[str, str]]:
    """Yield (parent, child) edges of the Tree dictionary structure."""
    children = tree.get(node)
    if not isinstance(children, dict):
        return
    for child in children:
        yield node, child
        yield from _iter_tree_edges(tree, child)


def _leaf_names(tree: Tree, node: str) -> List[str]:
    """Collect leaf (taxon) names under ``node``."""
    children = tree.get(node)
    if not isinstance(children, dict) or not children:
        return [node]
    leaves: List[str] = []
    for child in children:
        leaves.extend(_leaf_names(tree, child))
    return leaves


def _find_root(tree: Tree) -> str:
    """Find the root node (the node that is nobody's child)."""
    all_nodes = set(tree.keys())
    child_nodes: set[str] = set()
    for node_data in tree.values():
        if isinstance(node_data, dict):
            child_nodes.update(node_data.keys())
    return (all_nodes - child_nodes).pop()


def _fitch_set(
    tree: Tree,
    node: str,
    tip_states: Mapping[str, int],
) -> Tuple[set[int], int]:
    """Fitch (1971) postorder set algorithm for one binary character.

    Returns (candidate set, parsimony score for the subtree). Descriptive
    only: this counts the minimum number of state changes on the observed
    tree, with no rate model or probabilistic inference.
    """
    children = tree.get(node)
    if not isinstance(children, dict) or not children:
        return {tip_states[node]}, 0

    sets: List[set[int]] = []
    score = 0
    for child in children:
        child_set, child_score = _fitch_set(tree, child, tip_states)
        sets.append(child_set)
        score += child_score

    intersection = set.intersection(*sets)
    if intersection:
        return intersection, score
    union = set.union(*sets)
    return union, score + 1


def count_parsimony_transitions(
    newick: str,
    tip_states: Mapping[str, int],
) -> Dict[str, Any]:
    """Count minimal state changes (gains + losses combined) on a Newick tree.

    Descriptive Fitch parsimony on one binary character (e.g. 0 = not
    tissue-specific in this tissue, 1 = tissue-specific in this tissue).
    Species absent from the tree raise ``ValueError``; tree tips without a
    state raise ``ValueError`` (partial observation is handled by callers).

    Args:
        newick: Newick string of the species tree (must end with ``;``).
        tip_states: Taxon -> binary state mapping covering every tree tip.

    Returns:
        Dict with ``parsimony_score`` (total minimal transitions),
        ``root_state_candidates`` (the Fitch root candidate set), and
        ``scored_tips`` (number of tips scored).
    """
    tree = from_newick(newick)
    root = _find_root(tree)
    tree_leaves = set(_leaf_names(tree, root))

    missing = sorted(set(tip_states) - tree_leaves)
    if missing:
        raise ValueError(f"Tip states given for taxa absent from tree: {missing}")
    unscored = sorted(tree_leaves - set(tip_states))
    if unscored:
        raise ValueError(f"Tree tips missing tip states: {unscored}")

    scored = {taxon: int(state) for taxon, state in tip_states.items()}
    root_set, score = _fitch_set(tree, root, scored)

    return {
        "parsimony_score": score,
        "root_state_candidates": sorted(root_set),
        "scored_tips": len(scored),
    }


def count_tissue_specificity_gains_losses(
    newick: str,
    states_df: pd.DataFrame,
    tau_cutoff: float = DEFAULT_TAU_CUTOFF,
) -> pd.DataFrame:
    """Per-orthogroup, per-tissue descriptive gain/loss parsimony count.

    Builds a binary tip-state vector for one tissue (1 when the orthogroup's
    tau >= ``tau_cutoff`` in that species) and counts minimal transitions on
    the tree. Orthogroups observed in fewer than every tree tip return
    ``None`` scores (not inferred); orthogroups observed in fewer than two
    species also return ``None``.

    Args:
        newick: Newick species tree.
        states_df: Concatenated per-species state tables (rows indexed by
            orthogroup) with columns ``tau`` and ``species``.
        tau_cutoff: Tau cutoff for the binary state.

    Returns:
        DataFrame indexed by orthogroup with columns ``n_observed``,
        ``parsimony_score``, ``n_specific_species``, and ``descriptive``
        (constant ``True`` marker).
    """
    tree = from_newick(newick)
    n_tips = len(set(_leaf_names(tree, _find_root(tree))))

    records: List[Dict[str, Any]] = []
    for orthogroup, group in states_df.groupby(level=0):
        observed = group[group["tau"].notna()]
        n_specific = int((observed["tau"] >= tau_cutoff).sum())
        if len(observed) < 2 or len(observed) < n_tips:
            records.append(
                {
                    "orthogroup": orthogroup,
                    "n_observed": len(observed),
                    "parsimony_score": None,
                    "n_specific_species": n_specific,
                    "descriptive": True,
                }
            )
            continue
        tip_states = {str(row["species"]): int(row["tau"] >= tau_cutoff) for _, row in observed.iterrows()}
        counts = count_parsimony_transitions(newick, tip_states)
        records.append(
            {
                "orthogroup": orthogroup,
                "n_observed": counts["scored_tips"],
                "parsimony_score": counts["parsimony_score"],
                "n_specific_species": n_specific,
                "descriptive": True,
            }
        )

    return pd.DataFrame(records).set_index("orthogroup")


# =============================================================================
# Duplication coupling summary
# =============================================================================


def orthology_class_coupling_table(
    gains_df: pd.DataFrame,
    orthogroup_copy_counts: Mapping[str, int],
    multicopy_threshold: int = 2,
) -> pd.DataFrame:
    """Summarize tissue-specificity-gain status by orthology class.

    Descriptive contingency of orthogroup classes (single-copy vs multi-copy
    within the panel, per Mantica et al.'s 1:1 vs multicopy classes) against
    whether the orthogroup shows at least one tissue where the minimal
    transition count implies a state change on some branch (parsimony_score
    >= 1 with at least one specific species).

    Args:
        gains_df: Result of :func:`count_tissue_specificity_gains_losses`
            (or a concatenated multi-tissue equivalent with columns
            ``parsimony_score`` and ``n_specific_species`` and orthogroup in
            the index).
        orthogroup_copy_counts: Orthogroup -> number of panel genes.
        multicopy_threshold: Copy count at or above which an orthogroup is
            multi-copy (default 2).

    Returns:
        DataFrame indexed by orthology class with counts and the fraction of
        orthogroups with at least one implied transition.
    """
    classes: Dict[str, List[bool]] = {"single_copy": [], "multicopy": []}
    for orthogroup, group in gains_df.groupby(level=0):
        copies = int(orthogroup_copy_counts.get(orthogroup, 1))
        key = "multicopy" if copies >= multicopy_threshold else "single_copy"
        has_transition = bool(((group["parsimony_score"].fillna(0) >= 1) & (group["n_specific_species"] >= 1)).any())
        classes[key].append(has_transition)

    rows = []
    for label, values in classes.items():
        n = len(values)
        rows.append(
            {
                "orthology_class": label,
                "n_orthogroups": n,
                "n_with_transition": int(sum(values)),
                "fraction_with_transition": (sum(values) / n) if n else np.nan,
                "descriptive": True,
            }
        )
    return pd.DataFrame(rows).set_index("orthology_class")


# =============================================================================
# Species-tree input helpers (provenance-documented inputs)
# =============================================================================


def validate_species_tree_coverage(newick: str, species: Sequence[str]) -> Dict[str, Any]:
    """Check that cohort species appear as tips of the Newick species tree.

    Descriptive provenance check for the design documented in the child repo
    (doc/05_cross_species/TISSUE_SPECIFICITY_EVOLUTION_DESIGN.md): the panel
    tree must cover every cohort species before any gain/loss counting.

    Returns:
        Dict with ``missing_species``, ``n_covered``, and ``n_requested``.
    """
    tree = from_newick(newick)
    leaves = set(_leaf_names(tree, _find_root(tree)))

    species_list = list(species)
    missing = sorted(set(species_list) - leaves)
    covered = sorted(set(species_list) & leaves)
    return {
        "missing_species": missing,
        "n_covered": len(covered),
        "n_requested": len(species_list),
    }


__all__ = [
    "DEFAULT_TAU_CUTOFF",
    "assign_tissue_specificity_states",
    "compute_tau",
    "count_parsimony_transitions",
    "count_tissue_specificity_gains_losses",
    "orthology_class_coupling_table",
    "validate_species_tree_coverage",
]
