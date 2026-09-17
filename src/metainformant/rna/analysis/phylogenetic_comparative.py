"""Phylogenetic comparative methods (PGLS) on contract-validated species trees.

Implements MJ-03 of the hymenoptera_amalgkit campaign methods backlog:

1. Brownian variance-covariance (VCV) construction from a rooted,
   bifurcating species tree (Newick string or the nested-dict
   ``name``/``children``/``distance`` representation, both validated
   read-only through
   ``metainformant.rna.analysis.statistics_contract.validate_species_tree_invariants``).
2. PGLS (phylogenetic generalized least squares) regression of a trait or
   expression response on predictors, with Pagel's ``lambda`` covariance
   parameter estimated by restricted maximum likelihood (REML) or fixed
   by the caller.
3. Coefficient standard errors, t/p statistics from a t distribution with
   ``n - p`` residual degrees of freedom, and model diagnostics (REML
   log-likelihood, AIC, VCV conditioning, GLS residuals, Shapiro-Wilk
   residual-normality check).
4. Optional tree-uncertainty quantification via seeded resampling (with
   replacement) across a supplied set of input trees.

Fail-closed behavior:

- Tree structural validation is delegated to the shared statistics
  contract; :class:`TreeInvariantError` and :class:`ProvenanceError`
  propagate unchanged. Rootedness is CALLER-DECLARED PROVENANCE
  (statistical_analysis_plan.md section 5.3): PGLS requires ``rooted``
  to be declared from recorded tree provenance, because topology alone
  cannot establish biological rooting.
- Every non-root branch must carry an explicit finite non-negative
  branch length. A missing or invalid length fails closed with
  :class:`TreeInvariantError` even when the structural invariants pass
  (the contract does not require branch lengths; PGLS cannot proceed
  without them). The root's own ``distance`` is meaningless for a
  rooted tree and is ignored.
- Missing species, duplicate labels, non-finite values, collinear
  predictors, and designs without residual degrees of freedom raise
  ``ValueError``; no silent imputation, subsetting, or dropping.

Newick limitations (shared with the contract's validator): unquoted
names, optional internal labels, and plain decimal branch lengths.
Quoted labels and Newick comments are not supported.

References: Grafen (1989, Am Nat); Pagel (1999, Evolution);
Freckleton, Harvey & Pagel (2002, Evolution) for the lambda-PGLS
formulation; REML per Harville (1977).

BOUNDARY: the Brownian model here is a covariance model for continuous
trait regression, not a character-history reconstruction; no ancestral
states, divergence times, or stochastic-character-map inference are
produced by this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import linalg as _scipy_linalg
from scipy import optimize as _scipy_optimize
from scipy import stats as _scipy_stats

from metainformant.core.utils import logging
from metainformant.rna.analysis.statistics_contract import (
    TreeInvariantError,
    validate_species_tree_invariants,
)

logger = logging.get_logger(__name__)

TreeLike = str | Mapping[str, Any]
Tree = dict[str, Any]

#: Reserved coefficient name of the intercept column (when ``add_intercept``).
INTERCEPT_NAME = "intercept"

__all__ = [
    "INTERCEPT_NAME",
    "PGLSResult",
    "TreeLike",
    "brownian_vcv",
    "fit_pgls",
    "fit_pgls_tree_uncertainty",
    "lambda_adjusted_vcv",
    "prepare_species_tree",
    "simulate_brownian_traits",
]


# =============================================================================
# Tree parsing, validation, and normalization
# =============================================================================


def _checked_branch_length(node_name: str, value: Any) -> float:
    """Return ``value`` as a finite non-negative float, failing closed."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TreeInvariantError(
            f"branch length for node {node_name!r} must be a real number, got {value!r}"
        )
    distance = float(value)
    if not np.isfinite(distance) or distance < 0.0:
        raise TreeInvariantError(
            f"branch length for node {node_name!r} must be finite and non-negative, got {value!r}"
        )
    return distance


def _parse_newick_branches(newick: str) -> Tree:
    """Parse a plain Newick string into nested dicts with branch lengths.

    Returns nodes shaped like the contract's dict representation
    (``name``/``children``/``distance``). Unnamed internal nodes receive
    deterministic ``Internal_N`` labels. Unlike the contract validator
    (which treats branch lengths as optional), every non-root branch must
    carry a valid length: PGLS cannot proceed without one.

    The caller must already have run
    ``validate_species_tree_invariants`` on ``newick``; this parser then
    only fails on missing/invalid branch lengths or residual syntax
    defects, raising :class:`TreeInvariantError`.
    """
    text = newick.strip()
    if not text.endswith(";"):
        raise TreeInvariantError("Newick string must terminate with ';'")
    text = text[:-1]
    position = 0
    counter = [0]

    def parse_node() -> tuple[Tree, float | None]:
        nonlocal position
        parsed_children: list[tuple[Tree, float | None]] = []
        if position < len(text) and text[position] == "(":
            position += 1
            parsed_children.append(parse_node())
            while position < len(text) and text[position] == ",":
                position += 1
                parsed_children.append(parse_node())
            if position >= len(text) or text[position] != ")":
                raise TreeInvariantError(f"unbalanced parentheses at character {position}")
            position += 1
        start = position
        while position < len(text) and text[position] not in "(),:":
            position += 1
        name = text[start:position].strip()
        distance: float | None = None
        if position < len(text) and text[position] == ":":
            position += 1
            length_start = position
            while position < len(text) and text[position] not in "(),":
                position += 1
            length_text = text[length_start:position].strip()
            try:
                distance = float(length_text)
            except ValueError:
                raise TreeInvariantError(
                    f"invalid branch length {length_text!r} in Newick string"
                ) from None
        children: list[Tree] = []
        for child, child_distance in parsed_children:
            if child_distance is None:
                raise TreeInvariantError(
                    f"branch to node {child['name']!r} is missing a branch length"
                )
            child["distance"] = _checked_branch_length(child["name"], child_distance)
            children.append(child)
        if children and name == "":
            counter[0] += 1
            name = f"Internal_{counter[0]}"
        if not children and name == "":
            raise TreeInvariantError("Newick string contains an anonymous leaf label")
        return {"name": name, "children": children, "distance": distance}, distance

    root, _root_distance = parse_node()
    if position != len(text):
        raise TreeInvariantError(f"unexpected characters after Newick root at position {position}")
    return root


def _normalize_dict_tree(tree: Mapping[str, Any], *, is_root: bool) -> Tree:
    """Deep-copy a validated nested-dict tree, checking branch lengths."""
    name = tree["name"]
    children: list[Tree] = []
    for child in tree.get("children", []):
        children.append(_normalize_dict_tree(child, is_root=False))
    node: Tree = {"name": name, "children": children}
    distance = tree.get("distance")
    if is_root:
        # The root's own branch length carries no information about the
        # covariance structure of a rooted tree; record it but do not
        # validate or use it.
        node["distance"] = float(distance) if distance is not None else 0.0
    elif distance is None:
        raise TreeInvariantError(f"branch to node {name!r} is missing a branch length")
    else:
        node["distance"] = _checked_branch_length(name, distance)
    return node


def prepare_species_tree(
    tree: TreeLike,
    *,
    rooted: bool | None = None,
    require_bifurcating_root: bool = True,
) -> Tree:
    """Validate a species tree through the contract and normalize it.

    Accepts a Newick string or a nested-dict mapping (keys ``name``,
    ``children``, ``distance``) and returns a private deep copy in the
    dict representation with every branch length checked. Validation is
    fail-closed per the module docstring.

    Raises:
        TreeInvariantError: On malformed structure, invalid branch
            lengths, or a caller-declared-unrooted tree.
        ProvenanceError: When rootedness was not explicitly declared.
    """
    validate_species_tree_invariants(
        tree,
        rooted=rooted,
        require_rooted=True,
        require_bifurcating_root=require_bifurcating_root,
    )
    if isinstance(tree, str):
        normalized = _parse_newick_branches(tree)
    elif isinstance(tree, Mapping):
        normalized = _normalize_dict_tree(tree, is_root=True)
    else:
        raise TreeInvariantError("species tree must be a Newick string or a nested-dict tree mapping")
    return normalized


def brownian_vcv(
    tree: TreeLike,
    *,
    rooted: bool | None = None,
    require_bifurcating_root: bool = True,
) -> pd.DataFrame:
    """Brownian tip covariance matrix from a validated rooted tree.

    For a Brownian process on a rooted tree, ``V[i, j]`` is the path
    length shared by tips ``i`` and ``j``: the root-to-MRCA distance.
    Diagonal entries are each tip's root-to-tip depth. Root children
    share only the root (depth 0), so they are independent, as expected.

    Returns a tip-by-tip ``pandas.DataFrame`` indexed in deterministic
    depth-first leaf order.

    Raises:
        TreeInvariantError: On invalid tree input or branch lengths.
        ProvenanceError: When rootedness was not explicitly declared.
    """
    root = prepare_species_tree(
        tree, rooted=rooted, require_bifurcating_root=require_bifurcating_root
    )
    leaves: list[str] = []

    def collect(node: Mapping[str, Any]) -> None:
        children = node["children"]
        if not children:
            leaves.append(node["name"])
            return
        for child in children:
            collect(child)

    collect(root)
    n = len(leaves)
    vcv = np.zeros((n, n), dtype=float)
    index = {name: i for i, name in enumerate(leaves)}

    def visit(node: Mapping[str, Any], depth: float) -> list[int]:
        children = node["children"]
        if not children:
            i = index[node["name"]]
            vcv[i, i] = depth
            return [i]
        groups = [visit(child, depth + child["distance"]) for child in children]
        # Postorder: pairs within one child group were already assigned at
        # a deeper node; pairs across distinct groups have their MRCA at
        # THIS node, so their covariance is this node's root distance.
        for first in range(len(groups)):
            for second in range(first + 1, len(groups)):
                for i in groups[first]:
                    for j in groups[second]:
                        vcv[i, j] = depth
                        vcv[j, i] = depth
        return [i for group in groups for i in group]

    visit(root, 0.0)
    return pd.DataFrame(vcv, index=leaves, columns=leaves)


def lambda_adjusted_vcv(vcv: pd.DataFrame, lambda_: float) -> pd.DataFrame:
    """Apply Pagel's lambda to a Brownian tip covariance matrix.

    ``V_lambda = lambda * V + (1 - lambda) * diag(V)``: shared ancestral
    variance is scaled by ``lambda`` while each tip's total variance is
    preserved, so ``lambda = 1`` recovers pure Brownian covariance and
    ``lambda = 0`` yields a star phylogeny (independent tips).
    """
    if not isinstance(vcv, pd.DataFrame):
        raise TypeError("vcv must be a pandas DataFrame from brownian_vcv")
    if lambda_ is None or not np.isfinite(lambda_) or not 0.0 <= float(lambda_) <= 1.0:
        raise ValueError(f"lambda must be a finite value in [0, 1], got {lambda_!r}")
    values = vcv.to_numpy(dtype=float)
    adjusted = float(lambda_) * values + (1.0 - float(lambda_)) * np.diag(np.diag(values))
    return pd.DataFrame(adjusted, index=vcv.index, columns=vcv.columns)


# =============================================================================
# GLS core and REML lambda estimation
# =============================================================================


def _positive_definite_cholesky(v: np.ndarray) -> np.ndarray:
    """Cholesky factor of ``v``, failing closed on a non-PSD matrix."""
    try:
        return np.linalg.cholesky(v)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "phylogenetic covariance matrix is not positive definite; check the "
            "tree branch lengths (a tip with zero total root-to-tip depth makes "
            f"the Brownian VCV singular): {exc}"
        ) from None


def _gls_fit(v: np.ndarray, x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """GLS fit of ``y ~ x`` with covariance ``v`` (Brownian model, unit scale).

    Returns GLS coefficients, coefficient covariance at unit residual
    variance (``(X' V^-1 X)^-1``), the quadratic residual form, REML
    pieces, and log-determinants. Residual variance is rescaled by the
    caller's residual degrees of freedom.
    """
    chol = _positive_definite_cholesky(v)
    logdet_v = 2.0 * float(np.log(np.diag(chol)).sum())
    factor = _scipy_linalg.cho_factor(v, lower=True, check_finite=False)
    v_inverse_x = _scipy_linalg.cho_solve(factor, x)
    v_inverse_y = _scipy_linalg.cho_solve(factor, y)
    xt_v_inverse_x = x.T @ v_inverse_x
    sign, logdet_xt_v_inverse_x = np.linalg.slogdet(xt_v_inverse_x)
    if sign <= 0:
        raise ValueError(
            "PGLS design is rank deficient: X' V^-1 X is singular; remove "
            "collinear predictors or add observations"
        )
    xt_v_inverse_x_inverse = np.linalg.inv(xt_v_inverse_x)
    beta = xt_v_inverse_x_inverse @ (x.T @ v_inverse_y)
    residuals = y - x @ beta
    residual_quadratic = float(residuals @ _scipy_linalg.cho_solve(factor, residuals))
    return {
        "beta": beta,
        "xt_v_inverse_x_inverse": xt_v_inverse_x_inverse,
        "residuals": residuals,
        "residual_quadratic": residual_quadratic,
        "logdet_v": logdet_v,
        "logdet_xt_v_inverse_x": logdet_xt_v_inverse_x,
    }


def _reml_loglik(fit: Mapping[str, Any], residual_df: int) -> float:
    """REML log-likelihood of a fitted GLS model (Harville 1977)."""
    sigma2_reml = fit["residual_quadratic"] / residual_df
    if sigma2_reml <= 0.0:
        raise ValueError(
            "estimated residual variance is not positive; the model has no "
            "residual variation to estimate"
        )
    return float(
        -0.5
        * (
            residual_df * (np.log(2.0 * np.pi * sigma2_reml) + 1.0)
            + fit["logdet_v"]
            + fit["logdet_xt_v_inverse_x"]
        )
    )


def _profile_reml_objective(
    lambda_value: float, vcv_values: np.ndarray, x: np.ndarray, y: np.ndarray, residual_df: int
) -> float:
    """Negative REML profile log-likelihood at a fixed lambda (for minimize_scalar)."""
    adjusted = lambda_value * vcv_values + (1.0 - lambda_value) * np.diag(np.diag(vcv_values))
    fit = _gls_fit(adjusted, x, y)
    return -_reml_loglik(fit, residual_df)


def _estimate_lambda_reml(
    vcv_values: np.ndarray, x: np.ndarray, y: np.ndarray, residual_df: int
) -> float:
    """Maximize the REML profile log-likelihood over lambda in [0, 1]."""
    result = _scipy_optimize.minimize_scalar(
        _profile_reml_objective,
        args=(vcv_values, x, y, residual_df),
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-6},
    )
    if not result.success:
        raise ValueError(f"lambda REML optimization failed: {result.message}")
    return float(result.x)


# =============================================================================
# Design-matrix alignment and fail-closed data validation
# =============================================================================


def _validated_response(response: pd.Series) -> tuple[np.ndarray, list[str]]:
    """Validate a species-indexed response series, returning values and labels."""
    if not isinstance(response, pd.Series):
        raise TypeError("response must be a pandas Series indexed by species names")
    if response.index.has_duplicates:
        raise ValueError("response index contains duplicate species labels")
    labels: list[str] = []
    for label in response.index:
        if not isinstance(label, str) or not label.strip():
            raise ValueError(
                f"response index must contain non-empty species name strings, got {label!r}"
            )
        labels.append(label)
    try:
        values = response.to_numpy(dtype=float)
    except (TypeError, ValueError):
        raise ValueError("response values must be numeric") from None
    non_finite = [labels[i] for i, value in enumerate(values) if not np.isfinite(value)]
    if non_finite:
        raise ValueError(f"response contains missing or non-finite values for: {non_finite}")
    return values, labels


def _validated_predictors(predictors: pd.Series | pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    """Validate predictors against the response species, returning the aligned frame."""
    if isinstance(predictors, pd.Series):
        name = predictors.name
        column_name = name if isinstance(name, str) and name.strip() else "predictor"
        frame = predictors.to_frame(name=column_name)
    elif isinstance(predictors, pd.DataFrame):
        frame = predictors
    else:
        raise TypeError("predictors must be a pandas Series or DataFrame")
    if frame.columns.duplicated().any():
        raise ValueError("predictors contain duplicate column names")
    if frame.shape[1] == 0:
        raise ValueError("at least one predictor is required")
    columns: list[str] = []
    for column in frame.columns:
        if not isinstance(column, str) or not column.strip():
            raise ValueError(f"predictor column names must be non-empty strings, got {column!r}")
        columns.append(column)
    missing = [label for label in labels if label not in frame.index]
    if missing:
        raise ValueError(f"predictors are missing species present in the response index: {missing}")
    aligned = frame.loc[labels, columns]
    try:
        values = aligned.to_numpy(dtype=float)
    except (TypeError, ValueError):
        raise ValueError("predictor values must be numeric") from None
    for column_index, column in enumerate(columns):
        offenders = [
            labels[i]
            for i in range(len(labels))
            if not np.isfinite(values[i, column_index])
        ]
        if offenders:
            raise ValueError(f"predictor {column!r} contains missing or non-finite values for: {offenders}")
    return pd.DataFrame(values, index=labels, columns=columns)


def _require_species_covered(vcv: pd.DataFrame, labels: Sequence[str]) -> None:
    """Fail closed when a tree tip set does not cover every analysis species."""
    missing = [label for label in labels if label not in set(vcv.index)]
    if missing:
        raise ValueError(f"species absent from the species tree tips: {missing}")


# =============================================================================
# PGLS fitting
# =============================================================================


@dataclass(frozen=True)
class PGLSResult:
    """Fitted PGLS model (GLS under a lambda-adjusted Brownian covariance).

    Attributes:
        coefficients: GLS estimates, one per design-matrix column.
        standard_errors: ``sqrt(diag(sigma2 * (X' V^-1 X)^-1))``.
        t_statistics: coefficients divided by standard errors.
        p_values: two-sided t p-values with ``residual_df`` degrees of
            freedom (asymptotic; comparative panels are small).
        lambda_: the covariance parameter used (estimated or fixed).
        lambda_estimated: whether ``lambda_`` came from REML profiling.
        sigma2: REML residual variance ``e' V^-1 e / (n - p)``.
        log_likelihood: REML log-likelihood (comparable across models
            sharing the same fixed design).
        aic: ``2k - 2 loglik`` with ``k = p + 1 (+1 when lambda is
            estimated)``; comparable only within the same design.
        residual_df: ``n - p``.
        n_obs: number of species analyzed.
        covariance: the lambda-adjusted tip covariance actually used.
        brownian_vcv: the raw Brownian tip covariance of the tree.
        diagnostics: condition number, GLS residuals, Shapiro-Wilk
            residual-normality check (None when ``n < 3``), species list.
    """

    coefficients: pd.Series
    standard_errors: pd.Series
    t_statistics: pd.Series
    p_values: pd.Series
    lambda_: float
    lambda_estimated: bool
    sigma2: float
    log_likelihood: float
    aic: float
    residual_df: int
    n_obs: int
    covariance: pd.DataFrame
    brownian_vcv: pd.DataFrame
    diagnostics: dict[str, Any]


def fit_pgls(
    tree: TreeLike,
    response: pd.Series,
    predictors: pd.Series | pd.DataFrame,
    *,
    rooted: bool | None = None,
    lambda_: float | None = None,
    add_intercept: bool = True,
    require_bifurcating_root: bool = True,
) -> PGLSResult:
    """Fit a PGLS regression of ``response`` on ``predictors`` for ``tree``.

    The tree is validated fail-closed (see module docstring); the
    Brownian tip covariance is built, lambda-adjusted (Pagel's lambda
    estimated by REML when ``lambda_ is None``), and the model is fit by
    GLS. Standard errors derive from the REML residual variance and the
    GLS information matrix.

    Args:
        tree: Rooted, bifurcating species tree (Newick or nested dict)
            with explicit finite non-negative branch lengths.
        response: Trait or expression response, indexed by species names
            present in the tree.
        predictors: One predictor series or a species-indexed DataFrame.
        rooted: Caller-declared rootedness provenance; must be True
            (topology alone cannot establish rooting, so None fails
            closed with ``ProvenanceError``).
        lambda_: Fixed covariance parameter in [0, 1]; None estimates it
            by REML.
        add_intercept: Prepend an intercept column (coefficient name
            ``intercept``; a predictor column with that name collides and
            fails closed).
        require_bifurcating_root: Passed to the contract validator;
            internal polytomies are permitted and yield MRCA covariance
            at the polytomy's root distance.

    Raises:
        TreeInvariantError: On malformed trees or invalid branch lengths.
        ProvenanceError: When rootedness was not explicitly declared.
        ValueError: On invalid data, species not in the tree, collinear
            predictors, or a design without residual degrees of freedom.
    """
    vcv_brownian = brownian_vcv(
        tree, rooted=rooted, require_bifurcating_root=require_bifurcating_root
    )
    y, labels = _validated_response(response)
    design = _validated_predictors(predictors, labels)
    if add_intercept:
        if INTERCEPT_NAME in design.columns:
            raise ValueError(
                f"predictor column collides with the reserved {INTERCEPT_NAME!r} name; "
                "rename the column or pass add_intercept=False"
            )
        x = np.column_stack([np.ones(len(labels)), design.to_numpy(dtype=float)])
        names = [INTERCEPT_NAME, *design.columns]
    else:
        x = design.to_numpy(dtype=float)
        names = list(design.columns)
    n, p = x.shape
    residual_df = n - p
    if residual_df < 1:
        raise ValueError(
            f"PGLS needs at least one residual degree of freedom: got n={n} species "
            f"and p={p} parameters"
        )
    _require_species_covered(vcv_brownian, labels)
    vcv_values = vcv_brownian.loc[labels, labels].to_numpy(dtype=float)

    if lambda_ is None:
        lambda_estimated = True
        lambda_value = _estimate_lambda_reml(vcv_values, x, y, residual_df)
    else:
        lambda_estimated = False
        if not np.isfinite(lambda_) or not 0.0 <= float(lambda_) <= 1.0:
            raise ValueError(f"lambda must be a finite value in [0, 1], got {lambda_!r}")
        lambda_value = float(lambda_)
    covariance_values = lambda_value * vcv_values + (1.0 - lambda_value) * np.diag(
        np.diag(vcv_values)
    )
    fit = _gls_fit(covariance_values, x, y)
    sigma2 = fit["residual_quadratic"] / residual_df
    if sigma2 <= 0.0:
        raise ValueError(
            "estimated residual variance is not positive; the model has no "
            "residual variation to estimate"
        )
    log_likelihood = _reml_loglik(fit, residual_df)
    n_parameters = p + 1 + (1 if lambda_estimated else 0)
    aic = 2.0 * n_parameters - 2.0 * log_likelihood

    coefficient_covariance = sigma2 * fit["xt_v_inverse_x_inverse"]
    standard_errors = np.sqrt(np.diag(coefficient_covariance))
    t_statistics = fit["beta"] / standard_errors
    p_values = 2.0 * _scipy_stats.t.sf(np.abs(t_statistics), residual_df)

    diagnostics: dict[str, Any] = {
        "lambda_estimated": lambda_estimated,
        "covariance_condition_number": float(np.linalg.cond(covariance_values)),
        "gls_residuals": pd.Series(fit["residuals"], index=labels),
        "n_species": n,
    }
    if n >= 3:
        shapiro = _scipy_stats.shapiro(fit["residuals"])
        diagnostics["shapiro_wilk_statistic"] = float(shapiro.statistic)
        diagnostics["shapiro_wilk_pvalue"] = float(shapiro.pvalue)
    else:
        diagnostics["shapiro_wilk_statistic"] = None
        diagnostics["shapiro_wilk_pvalue"] = None

    index = pd.Index(names)
    return PGLSResult(
        coefficients=pd.Series(fit["beta"], index=index),
        standard_errors=pd.Series(standard_errors, index=index),
        t_statistics=pd.Series(t_statistics, index=index),
        p_values=pd.Series(p_values, index=index),
        lambda_=lambda_value,
        lambda_estimated=lambda_estimated,
        sigma2=sigma2,
        log_likelihood=log_likelihood,
        aic=aic,
        residual_df=residual_df,
        n_obs=n,
        covariance=lambda_adjusted_vcv(vcv_brownian, lambda_value),
        brownian_vcv=vcv_brownian,
        diagnostics=diagnostics,
    )


# =============================================================================
# Tree-uncertainty resampling
# =============================================================================


def fit_pgls_tree_uncertainty(
    trees: Sequence[TreeLike],
    response: pd.Series,
    predictors: pd.Series | pd.DataFrame,
    *,
    n_resamples: int = 1000,
    seed: int | None = None,
    ci: tuple[float, float] = (0.025, 0.975),
    rooted: bool | None = None,
    add_intercept: bool = True,
    require_bifurcating_root: bool = True,
) -> dict[str, Any]:
    """Quantify tree uncertainty by seeded resampling across input trees.

    Every supplied tree is validated fail-closed up front (and must cover
    every analysis species), then trees are resampled WITH replacement
    ``n_resamples`` times using ``numpy.random.default_rng(seed)`` and
    PGLS is refit per resample (lambda re-estimated each time). Returns
    per-coefficient draws, mean, standard deviation, and percentile CI.

    Raises:
        TreeInvariantError / ProvenanceError: On any invalid tree.
        ValueError: On invalid data, empty tree set, non-positive
            ``n_resamples``, or an invalid CI interval.
    """
    if len(trees) == 0:
        raise ValueError("at least one tree is required for tree-uncertainty resampling")
    if not isinstance(n_resamples, int) or isinstance(n_resamples, bool) or n_resamples < 1:
        raise ValueError(f"n_resamples must be a positive integer, got {n_resamples!r}")
    ci_low, ci_high = float(ci[0]), float(ci[1])
    if not (np.isfinite(ci_low) and np.isfinite(ci_high)) or not (
        0.0 <= ci_low < ci_high <= 1.0
    ):
        raise ValueError(f"ci must be (low, high) with 0 <= low < high <= 1, got {ci!r}")

    labels = _validated_response(response)[1]
    for tree in trees:
        _require_species_covered(
            brownian_vcv(tree, rooted=rooted, require_bifurcating_root=require_bifurcating_root),
            labels,
        )

    full = fit_pgls(
        trees[0],
        response,
        predictors,
        rooted=rooted,
        add_intercept=add_intercept,
        require_bifurcating_root=require_bifurcating_root,
    )
    names = list(full.coefficients.index)
    rng = np.random.default_rng(seed)
    draws = np.empty((len(names), n_resamples), dtype=float)
    lambdas = np.empty(n_resamples, dtype=float)
    for resample in range(n_resamples):
        tree = trees[int(rng.integers(0, len(trees)))]
        result = fit_pgls(
            tree,
            response,
            predictors,
            rooted=rooted,
            add_intercept=add_intercept,
            require_bifurcating_root=require_bifurcating_root,
        )
        draws[:, resample] = result.coefficients.to_numpy()
        lambdas[resample] = result.lambda_

    draws_frame = pd.DataFrame(draws, index=names)
    summary = pd.DataFrame(
        {
            "mean": draws_frame.mean(axis=1),
            "sd": draws_frame.std(axis=1, ddof=1),
            "ci_low": draws_frame.quantile(ci_low, axis=1),
            "ci_high": draws_frame.quantile(ci_high, axis=1),
        },
        index=names,
    )
    return {
        "n_resamples": n_resamples,
        "seed": seed,
        "ci": (ci_low, ci_high),
        "coefficient_names": names,
        "draws": draws_frame,
        "summary": summary,
        "lambda_mean": float(np.mean(lambdas)),
        "lambda_sd": float(np.std(lambdas, ddof=1)),
    }


# =============================================================================
# Brownian trait simulation (test and power-analysis utility)
# =============================================================================


def simulate_brownian_traits(
    tree: TreeLike,
    sigma: float = 1.0,
    *,
    n_traits: int = 1,
    trait_names: Sequence[str] | None = None,
    seed: int | None = None,
    rooted: bool | None = None,
    require_bifurcating_root: bool = True,
) -> pd.DataFrame:
    """Simulate traits evolving by Brownian motion on a validated tree.

    Draws from ``N(0, sigma^2 V)`` where ``V`` is the tree's Brownian tip
    covariance; species appear in the same depth-first order as
    ``brownian_vcv``. Deterministic for a fixed ``seed``.

    Raises:
        TreeInvariantError / ProvenanceError: On invalid tree input.
        ValueError: On invalid ``sigma``, ``n_traits``, or names.
    """
    if sigma is None or not np.isfinite(sigma) or sigma < 0.0:
        raise ValueError(f"sigma must be a finite non-negative number, got {sigma!r}")
    if not isinstance(n_traits, int) or isinstance(n_traits, bool) or n_traits < 1:
        raise ValueError(f"n_traits must be a positive integer, got {n_traits!r}")
    vcv = brownian_vcv(tree, rooted=rooted, require_bifurcating_root=require_bifurcating_root)
    names: list[str] | None = None
    if trait_names is not None:
        names = list(trait_names)
        if len(names) != n_traits:
            raise ValueError(
                f"trait_names has {len(names)} entries but n_traits={n_traits}"
            )
        if len(set(names)) != len(names):
            raise ValueError("trait_names must be unique")
        for name in names:
            if not isinstance(name, str) or not name.strip():
                raise ValueError(f"trait names must be non-empty strings, got {name!r}")
    else:
        names = [f"trait_{i}" for i in range(n_traits)]
    n = len(vcv.index)
    if sigma == 0.0:
        draws = np.zeros((n, n_traits), dtype=float)
    else:
        rng = np.random.default_rng(seed)
        draws = rng.multivariate_normal(
            np.zeros(n), (sigma**2) * vcv.to_numpy(dtype=float), size=n_traits
        ).T
    return pd.DataFrame(draws, index=vcv.index, columns=names)
