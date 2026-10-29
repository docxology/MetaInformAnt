"""Phenotype statistical analysis functions.

This module provides robust statistical testing and summary generation
for phenotypic traits, supporting both parametric and non-parametric
methods. Included are tools for ANOVA, Kruskal-Wallis, independent T-tests,
and comprehensive summary statistics.
"""

from __future__ import annotations

import math
from typing import Any, Dict

import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.formula.api import ols


def perform_linear_regression(df: pd.DataFrame, x_col: str, y_col: str) -> Dict[str, Any]:
    """Perform simple linear regression between two continuous features.

    Args:
        df: Input DataFrame
        x_col: Predictor variable column (continuous)
        y_col: Response variable column (continuous)

    Returns:
        Dictionary with slope, intercept, r_value, p_value, std_err. An ``{"error": ...}`` dictionary instead
        when a column is missing or not numeric, fewer than 3 rows are complete, or the predictor is constant
        (its slope is then not defined; the rank-deficient fit that was reported attributed the intercept's
        effect to the predictor and could call it significant).
    """
    if df.empty or x_col not in df.columns or y_col not in df.columns:
        return {"error": "Invalid input data or missing columns."}
    if not (pd.api.types.is_numeric_dtype(df[x_col]) and pd.api.types.is_numeric_dtype(df[y_col])):
        return {"error": "Both columns must be numeric for linear regression."}

    valid_data = df[[x_col, y_col]].dropna()
    if len(valid_data) < 3:
        return {"error": "Requires at least 3 valid observations for regression."}
    if valid_data[x_col].nunique() < 2:
        return {"error": "The predictor is constant, so its slope is not defined."}

    # Use statsmodels for more comprehensive regression results. Q() quotes the column names, which may hold
    # spaces or other characters a formula cannot parse bare.
    x_term = f"Q({x_col!r})"
    formula = f"Q({y_col!r}) ~ {x_term}"
    model = ols(formula, data=valid_data).fit()

    return {
        "test": "Linear Regression",
        "predictor": x_col,
        "response": y_col,
        "slope": float(model.params[x_term]),
        "intercept": float(model.params["Intercept"]),
        "r_squared": float(model.rsquared),
        "p_value": float(model.pvalues[x_term]),
        "std_err": float(model.bse[x_term]),
        "n_obs": len(valid_data),
        "significant": bool(model.pvalues[x_term] < 0.05),
    }


def perform_multifactor_anova(df: pd.DataFrame, formula: str) -> Dict[str, Any]:
    """Perform a multi-factor ANOVA using Ordinary Least Squares regression.

    Formula usually takes the form: 'dependent_var ~ C(factor1) * C(factor2)'
    Returns the ANOVA table containing F-statistics and P-values for
    main effects and interaction terms.

    Type II sums of squares (appropriate for unbalanced designs) are used. If they cannot be computed for the
    model, sequential Type I sums of squares are used instead, whose F and p depend on the order of the terms in
    the formula and can differ from Type II for unbalanced data; every term reports which one it got in
    ``"anova_type"`` (2 or 1) so the fallback is not silent.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", message="covariance of constraints does not have full rank")

        model = None
        anova_type = 2
        try:
            model = ols(formula, data=df).fit()
            # Type 2 ANOVA is standard for unbalanced designs
            anova_table = sm.stats.anova_lm(model, typ=2)
        except Exception:
            if model is None:
                return {"error": "Model could not be fit."}
            try:
                anova_type = 1
                anova_table = sm.stats.anova_lm(model, typ=1)
            except Exception:
                return {"error": "Model could not be fit."}

    results = {}
    for index, row in anova_table.iterrows():
        if index != "Residual":
            results[str(index)] = {"F": float(row["F"]), "P_value": float(row["PR(>F)"]), "anova_type": anova_type}
    return results


def get_comprehensive_pairwise_ttests(df: pd.DataFrame, value_col: str, group_col: str) -> list[Dict[str, Any]]:
    """Compute pairwise Welch's T-tests across all pairs of groups for annotation.

    Every pair with at least two observations in each group is tested; nothing is filtered by significance
    and no multiple-comparison correction is applied (compare each ``p_value`` with ``alpha / n_pairs``, or
    adjust them, before calling a pair different).

    Args:
        df: Input DataFrame
        value_col: Trait
        group_col: Categorical groups (labels of any type, e.g. integers, are matched as they are)

    Returns:
        List of result dictionaries from :func:`perform_ttest`, one per testable pair, sorted by p-value
        (pairs whose p-value is undefined, e.g. two constant groups, last).
    """
    groups = df[group_col].dropna().unique()
    results = []

    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            # The labels are passed on as they are: str() of an integer label never equals the column's values,
            # so every pair of a numerically labelled factor used to be dropped without a word
            test_res = perform_ttest(df, value_col, group_col, groups[i], groups[j])
            if "error" not in test_res:
                results.append(test_res)

    # Sort by p-value, undefined ones last
    results.sort(key=lambda x: (math.isnan(x["p_value"]), x["p_value"]))
    return results


def calculate_summary_stats(df: pd.DataFrame, value_col: str, group_col: str) -> pd.DataFrame:
    """Calculate descriptive statistics for a quantitative trait broken down by group.

    Args:
        df: Input DataFrame
        value_col: Name of the column containing the quantitative trait
        group_col: Name of the column defining the groups

    Returns:
        DataFrame containing mean, median, std, min, max, and count for each group
    """
    if df.empty or value_col not in df.columns or group_col not in df.columns:
        return pd.DataFrame()

    summary = (
        df.groupby(group_col)[value_col]
        .agg(mean="mean", median="median", std="std", min="min", max="max", count="count")
        .reset_index()
    )
    return summary


def perform_anova(df: pd.DataFrame, value_col: str, group_col: str) -> Dict[str, Any]:
    """Perform one-way ANOVA across groups for a given trait.

    Args:
        df: Input DataFrame
        value_col: Name of the column containing the quantitative trait
        group_col: Name of the column defining the groups

    Returns:
        Dictionary with f_statistic, p_value, and group specific data
    """
    if df.empty or value_col not in df.columns or group_col not in df.columns:
        return {"error": "Invalid input data or missing columns."}

    groups = df[group_col].dropna().unique()
    if len(groups) < 2:
        return {"error": "Requires at least 2 groups for ANOVA."}

    data = [df[df[group_col] == g][value_col].dropna().values for g in groups]
    # Filter out empty arrays
    data = [d for d in data if len(d) > 0]

    if len(data) < 2:
        return {"error": "Requires at least 2 groups with non-null values for ANOVA."}

    f_stat, p_val = stats.f_oneway(*data)

    return {
        "test": "One-way ANOVA",
        "f_statistic": float(f_stat),
        "p_value": float(p_val),
        "groups": len(data),
        "significant": bool(p_val < 0.05),
    }


def perform_kruskal(df: pd.DataFrame, value_col: str, group_col: str) -> Dict[str, Any]:
    """Perform non-parametric Kruskal-Wallis H-test across groups.

    Args:
        df: Input DataFrame
        value_col: Name of the column containing the quantitative trait
        group_col: Name of the column defining the groups

    Returns:
        Dictionary with h_statistic, p_value, and group specific data
    """
    if df.empty or value_col not in df.columns or group_col not in df.columns:
        return {"error": "Invalid input data or missing columns."}

    groups = df[group_col].dropna().unique()
    if len(groups) < 2:
        return {"error": "Requires at least 2 groups for Kruskal-Wallis."}

    data = [df[df[group_col] == g][value_col].dropna().values for g in groups]
    data = [d for d in data if len(d) > 0]

    if len(data) < 2:
        return {"error": "Requires at least 2 groups with non-null values."}

    h_stat, p_val = stats.kruskal(*data)

    return {
        "test": "Kruskal-Wallis H-test",
        "h_statistic": float(h_stat),
        "p_value": float(p_val),
        "groups": len(data),
        "significant": bool(p_val < 0.05),
    }


def perform_ttest(df: pd.DataFrame, value_col: str, group_col: str, group1: Any, group2: Any) -> Dict[str, Any]:
    """Perform independent T-test between two specific groups.

    Args:
        df: Input DataFrame
        value_col: Name of the column containing the quantitative trait
        group_col: Name of the column defining the groups
        group1: Label of the first group to compare (any type the group column holds)
        group2: Label of the second group to compare

    Returns:
        Dictionary with t_statistic, p_value, and relevance; an ``{"error": ...}`` dictionary when a column is
        missing or either group has fewer than two non-missing values (Welch's test needs a variance for each
        group; one value gave a NaN statistic and p-value reported as a result).
    """
    if df.empty or value_col not in df.columns or group_col not in df.columns:
        return {"error": "Invalid input data or missing columns."}

    g1_data = df[df[group_col] == group1][value_col].dropna()
    g2_data = df[df[group_col] == group2][value_col].dropna()

    if len(g1_data) < 2 or len(g2_data) < 2:
        return {
            "error": (
                f"Insufficient data for groups {group1} or {group2}: Welch's t-test needs at least 2 values in each."
            )
        }

    t_stat, p_val = stats.ttest_ind(g1_data, g2_data, equal_var=False)  # Welch's t-test by default

    return {
        "test": "Welch's t-test",
        "t_statistic": float(t_stat),
        "p_value": float(p_val),
        "group1": group1,
        "group2": group2,
        "group1_n": len(g1_data),
        "group2_n": len(g2_data),
        "significant": bool(p_val < 0.05),
    }


def correlate_phenotypes(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Calculate Pearson correlation matrix for numerical phenotype columns.

    Args:
        df: Input DataFrame
        cols: List of column names to correlate

    Returns:
        DataFrame containing the correlation matrix
    """
    valid_cols = [c for c in cols if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
    if len(valid_cols) < 2:
        return pd.DataFrame()

    return df[valid_cols].corr(method="pearson")
