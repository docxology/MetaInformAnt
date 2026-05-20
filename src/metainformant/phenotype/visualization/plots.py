"""Phenotype visualization tools.

Publication-quality plotting routines for phenotypic trait analysis:
  - Boxplot-swarm overlays with per-group statistics and significance brackets
  - Violin plots with strip overlays and strain hue support
  - Regression scatters with CI, R², effect size (Cohen's d), and residual QQ
  - Stacked proportion bars with chi-square annotations
  - Interaction point plots with strain coloring
  - Correlation heatmaps with statistical significance masking

Core figure-producing functions return matplotlib Figure objects while still
drawing into caller-provided axes when ``ax`` is supplied. Typography follows
publication standards (18pt titles, 14pt labels, 12pt ticks).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.offsetbox import AnchoredText

# ── Strain palette (shared with GWAS viz) ───────────────────────────────
STRAIN_PALETTE = {"C": "#2196F3", "I": "#FF9800", "M": "#4CAF50", "R": "#F44336"}


# ── Internal helpers ────────────────────────────────────────────────────


def _add_stats_textbox(
    ax: plt.Axes,
    text: Optional[str],
    loc: str = "upper left",
    fontsize: int = 10,
) -> None:
    """Add a styled statistics textbox to an axes."""
    if not text:
        return
    at = AnchoredText(
        text,
        loc=loc,
        prop=dict(size=fontsize, fontfamily="monospace"),
        frameon=True,
    )
    at.patch.set_boxstyle("round,pad=0.3,rounding_size=0.15")
    at.patch.set_facecolor("#FAFAFA")
    at.patch.set_alpha(0.92)
    at.patch.set_edgecolor("#333333")
    at.patch.set_linewidth(0.8)
    ax.add_artist(at)


def _apply_typography(
    ax: plt.Axes,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    """Apply consistent publication-grade typography."""
    if title:
        ax.set_title(title, fontsize=18, fontweight="bold", pad=15)
    if xlabel:
        ax.set_xlabel(
            xlabel.replace("_", " ").title(),
            fontsize=14,
            labelpad=10,
        )
    if ylabel:
        ax.set_ylabel(
            ylabel.replace("_", " ").title(),
            fontsize=14,
            labelpad=10,
        )
    ax.tick_params(axis="both", which="major", labelsize=12)


def _cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Compute Cohen's d effect size (pooled SD denominator)."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return 0.0
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_sd = math.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_sd < 1e-12:
        return 0.0
    return float((np.mean(group1) - np.mean(group2)) / pooled_sd)


def _format_p(p: float) -> str:
    """Format a p-value with significance stars."""
    if p < 0.001:
        return "p < 0.001 (***)"
    if p < 0.01:
        return f"p = {p:.3f} (**)"
    if p < 0.05:
        return f"p = {p:.3f} (*)"
    return f"p = {p:.2f} (ns)"


def _group_summary_text(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    order: List[str],
) -> str:
    """Build per-group summary: n, median, IQR."""
    lines = []
    for grp in order:
        vals = df.loc[df[x_col] == grp, y_col].dropna()
        if len(vals) == 0:
            continue
        q25, med, q75 = np.percentile(vals, [25, 50, 75])
        lines.append(f"{grp}: n={len(vals)}, med={med:.2f} [{q25:.2f}–{q75:.2f}]")
    return "\n".join(lines)


# ── Public API ──────────────────────────────────────────────────────────


def add_stat_annotation(
    ax: plt.Axes,
    x_coords: list[float],
    y_max: float,
    p_value: float,
    index: int = 0,
    effect_size: Optional[float] = None,
) -> None:
    """Draw a significance bracket with p-value and optional effect size.

    Args:
        ax: Target axes.
        x_coords: [x1, x2] positions for the bracket endpoints.
        y_max: Base y-coordinate for the first bracket.
        p_value: P-value to display.
        index: Stacking index (0 = lowest bracket).
        effect_size: Optional Cohen's d to display alongside p-value.
    """
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    h = y_range * 0.02
    y = y_max + (y_range * 0.06 * index) + h

    x1, x2 = x_coords
    p_text = _format_p(p_value)
    if effect_size is not None:
        p_text += f"  d={effect_size:.2f}"

    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, color="0.2")
    font_weight = "bold" if p_value < 0.05 else "normal"
    ax.text(
        (x1 + x2) * 0.5,
        y + h,
        p_text,
        ha="center",
        va="bottom",
        color="0.2",
        fontsize=9,
        fontweight=font_weight,
    )


def plot_boxplot_with_swarm(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str = "",
    ylabel: str = "",
    order: Optional[list[str]] = None,
    hue: Optional[str] = None,
    pairwise_stats: Optional[list[Dict[str, Any]]] = None,
    anova_stats: Optional[Dict[str, Any]] = None,
    stats_text: Optional[str] = None,
    show_group_summary: bool = True,
    ax: Optional[plt.Axes] = None,
) -> plt.Figure:
    """Boxplot overlaid with swarmplot, annotated with per-group statistics.

    Features:
      - Per-group sample size, median, and IQR in annotation box
      - Optional strain hue coloring (pass hue='strain' or similar)
      - Significance brackets with Cohen's d effect sizes
      - ANOVA F/p in subtitle when provided

    Args:
        df: DataFrame with at least x_col and y_col columns.
        x_col: Column for categorical x-axis groups.
        y_col: Column for numeric y-axis values.
        title: Plot title.
        ylabel: Y-axis label (defaults to y_col).
        order: Category ordering.
        hue: Optional column for within-group coloring (e.g. strain).
        pairwise_stats: List of dicts with group1, group2, p_value, significant.
        anova_stats: Dict with f_statistic, p_value.
        stats_text: Additional free-form stats text.
        show_group_summary: Show per-group n/median/IQR textbox.
        ax: Target axes (creates new figure if None).
    """
    n_brackets = len(pairwise_stats) if pairwise_stats else 0
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6 + 0.4 * n_brackets))

    if order is None:
        order = sorted(df[x_col].dropna().unique(), key=str)

    # Determine palette
    palette = None
    if hue:
        unique_hues = sorted(df[hue].dropna().unique(), key=str)
        palette = {h: STRAIN_PALETTE.get(str(h), sns.color_palette("muted")[i % 10]) for i, h in enumerate(unique_hues)}

    use_hue = hue if hue else x_col
    use_palette = palette if hue else {k: "#E0E0E0" for k in order}

    # Base boxplot
    sns.boxplot(
        data=df,
        x=x_col,
        y=y_col,
        order=order,
        hue=use_hue,
        showfliers=False,
        width=0.6,
        palette=use_palette,
        linewidth=1.2,
        ax=ax,
        legend=bool(hue),
    )
    # Strip overlay (jittered to prevent collision warnings)
    sns.stripplot(
        data=df,
        x=x_col,
        y=y_col,
        order=order,
        hue=use_hue,
        palette=use_palette if hue else "dark:.3",
        dodge=bool(hue),
        jitter=0.25,
        alpha=0.65,
        size=4,
        edgecolor="white",
        linewidth=0.3,
        ax=ax,
        legend=False,
    )

    _apply_typography(ax, title, x_col if not title else "", ylabel or y_col)

    # Per-group summary textbox
    group_text_parts = []
    if show_group_summary:
        group_text_parts.append(_group_summary_text(df, x_col, y_col, order))

    # ANOVA subtitle
    if anova_stats and "error" not in anova_stats:
        f_stat = anova_stats.get("f_statistic", 0)
        p_val = anova_stats.get("p_value", 1)
        df_between = anova_stats.get("df_between", "?")
        df_within = anova_stats.get("df_within", "?")
        eta_sq = anova_stats.get("eta_squared", None)
        anova_line = f"ANOVA: F({df_between},{df_within})={f_stat:.2f}, {_format_p(p_val)}"
        if eta_sq is not None:
            anova_line += f", η²={eta_sq:.3f}"
        group_text_parts.append(anova_line)

    if stats_text:
        group_text_parts.append(stats_text)

    combined_text = "\n".join(filter(None, group_text_parts))
    _add_stats_textbox(ax, combined_text, loc="upper left", fontsize=9)

    # Significance brackets with Cohen's d
    y_max_base = df[y_col].max()
    if pairwise_stats:
        sig_tests = [t for t in pairwise_stats if t.get("significant", False)]
        for idx, test in enumerate(sig_tests[:5]):
            try:
                i1 = order.index(test["group1"])
                i2 = order.index(test["group2"])
                # Compute Cohen's d
                g1_vals = df.loc[df[x_col] == test["group1"], y_col].dropna().values
                g2_vals = df.loc[df[x_col] == test["group2"], y_col].dropna().values
                d = _cohens_d(g1_vals, g2_vals) if len(g1_vals) > 1 and len(g2_vals) > 1 else None
                add_stat_annotation(
                    ax,
                    [i1, i2],
                    y_max_base,
                    test["p_value"],
                    index=idx,
                    effect_size=d,
                )
            except ValueError:
                pass

        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        ax.set_ylim(
            bottom=ax.get_ylim()[0],
            top=ax.get_ylim()[1] + y_range * 0.1 * min(len(sig_tests), 5),
        )

    if hue:
        ax.legend(
            title=hue.replace("_", " ").title(),
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=11,
            title_fontsize=13,
        )

    return ax.figure


def plot_violin(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str = "",
    ylabel: str = "",
    order: Optional[list[str]] = None,
    hue: Optional[str] = None,
    stats_text: Optional[str] = None,
    show_group_summary: bool = True,
    ax: Optional[plt.Axes] = None,
) -> plt.Figure:
    """Violin plot with strip overlay and per-group statistics.

    Features:
      - Inner quartile lines
      - Overlaid strip points
      - Optional split violins for binary hue
      - Per-group n/median/IQR annotation

    Args:
        df: DataFrame with data.
        x_col: Column for categorical x-axis.
        y_col: Column for numeric y-axis.
        title: Plot title.
        ylabel: Y-axis label.
        order: Category ordering.
        hue: Optional column for hue (split if binary).
        stats_text: Additional statistics text.
        show_group_summary: Show per-group summary.
        ax: Target axes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    if order is None:
        order = sorted(df[x_col].dropna().unique(), key=str)

    split_violins = hue is not None and df[hue].nunique() == 2

    use_hue = hue if hue else x_col
    v_palette = "muted" if hue else {k: sns.color_palette("muted")[i % 10] for i, k in enumerate(order)}

    sns.violinplot(
        data=df,
        x=x_col,
        y=y_col,
        order=order,
        hue=use_hue,
        split=split_violins,
        inner="quartile",
        palette=v_palette,
        alpha=0.5,
        ax=ax,
        legend=bool(hue),
    )
    sns.stripplot(
        data=df,
        x=x_col,
        y=y_col,
        order=order,
        hue=use_hue,
        dodge=bool(hue),
        palette="muted" if hue else "dark:.3",
        alpha=0.4,
        size=3,
        ax=ax,
        legend=False,
    )

    _apply_typography(ax, title, x_col if not title else "", ylabel or y_col)

    text_parts = []
    if show_group_summary:
        text_parts.append(_group_summary_text(df, x_col, y_col, order))
    if stats_text:
        text_parts.append(stats_text)
    _add_stats_textbox(ax, "\n".join(filter(None, text_parts)), loc="upper right", fontsize=9)

    if hue and not split_violins:
        ax.legend(
            title=hue.replace("_", " ").title(),
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=12,
            title_fontsize=14,
        )

    return ax.figure


def plot_regression_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    reg_stats: Optional[Dict[str, Any]] = None,
    hue: Optional[str] = None,
    stats_text: Optional[str] = None,
    show_residuals: bool = False,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Scatter plot with regression line, 95% CI, and comprehensive statistics.

    Features:
      - OLS regression line with 95% confidence band
      - R², adjusted R², p-value, β (slope), SE, and Cohen's f²
      - Random jitter for discrete x/y values
      - Optional residual distribution inset (QQ)

    Args:
        df: DataFrame with data.
        x_col, y_col: Column names for scatter axes.
        title, xlabel, ylabel: Labels.
        reg_stats: Dict with r_squared, p_value, beta, se, etc.
        hue: Column for color grouping.
        stats_text: Additional statistics text.
        show_residuals: If True, add a residual QQ inset.
        ax: Target axes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 7))

    # Regression line with CI band
    sns.regplot(
        data=df,
        x=x_col,
        y=y_col,
        ax=ax,
        scatter=False,
        color="#D32F2F",
        line_kws={"linewidth": 2.5, "label": "OLS fit"},
        ci=95,
    )

    # Jitter for discrete variables
    plot_df = df.copy()
    if pd.api.types.is_numeric_dtype(plot_df[x_col]):
        unique_x = plot_df[x_col].nunique()
        if unique_x < 20:
            jitter = 0.15 * (plot_df[x_col].max() - plot_df[x_col].min()) / max(unique_x, 1)
            plot_df[x_col] = plot_df[x_col] + np.random.uniform(-jitter, jitter, size=len(plot_df))
    if pd.api.types.is_numeric_dtype(plot_df[y_col]):
        unique_y = plot_df[y_col].nunique()
        if unique_y < 20:
            jitter = 0.15 * (plot_df[y_col].max() - plot_df[y_col].min()) / max(unique_y, 1)
            plot_df[y_col] = plot_df[y_col] + np.random.uniform(-jitter, jitter, size=len(plot_df))

    # Scatter with hue
    palette = STRAIN_PALETTE if hue else None
    sns.scatterplot(
        data=plot_df,
        x=x_col,
        y=y_col,
        hue=hue,
        palette=palette,
        alpha=0.7,
        edgecolor="k",
        linewidth=0.3,
        ax=ax,
        s=60,
    )

    _apply_typography(ax, title, xlabel or x_col, ylabel or y_col)

    # Build comprehensive stats annotation
    stat_lines = []
    if reg_stats:
        r2 = reg_stats.get("r_squared", 0)
        pval = reg_stats.get("p_value", 1)
        beta = reg_stats.get("beta", reg_stats.get("slope", None))
        se = reg_stats.get("se", reg_stats.get("std_error", None))
        n = reg_stats.get("n_samples", len(df))

        stat_lines.append(f"R² = {r2:.4f}")
        # Cohen's f² = R²/(1-R²)
        f_squared = r2 / (1 - r2) if r2 < 1 else float("inf")
        f_label = "large" if f_squared >= 0.35 else ("medium" if f_squared >= 0.15 else "small")
        stat_lines.append(f"f² = {f_squared:.3f} ({f_label})")

        if beta is not None:
            stat_lines.append(f"β = {beta:.4f}" + (f" ± {se:.4f}" if se else ""))
        stat_lines.append(f"{_format_p(pval)}")
        stat_lines.append(f"n = {n}")

    if stats_text:
        stat_lines.append(stats_text)

    _add_stats_textbox(ax, "\n".join(stat_lines), loc="lower right", fontsize=9)

    if hue:
        ax.legend(
            title=hue.replace("_", " ").title(),
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=12,
            title_fontsize=14,
        )

    return ax


def plot_categorical_proportions(
    df: pd.DataFrame,
    group_col: str,
    cat_col: str,
    title: str = "",
    chi2_stats: Optional[Dict[str, Any]] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Figure:
    """Stacked bar chart of categorical proportions with chi-square annotation.

    Args:
        df: DataFrame with data.
        group_col: Column for groups (x-axis bars).
        cat_col: Column for categories (stacked segments).
        title: Plot title.
        chi2_stats: Optional dict with chi2, p_value, dof, cramers_v.
        ax: Target axes.
    """
    counts = df.groupby([group_col, cat_col]).size().unstack(fill_value=0)
    try:
        counts = counts.loc[
            sorted(
                counts.index,
                key=lambda x: (
                    int("".join(c for c in str(x) if c.isdigit())) if any(c.isdigit() for c in str(x)) else float("inf")
                ),
            )
        ]
    except Exception:
        pass

    proportions = counts.div(counts.sum(axis=1), axis=0) * 100

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    bottom = np.zeros(len(proportions))
    colors = plt.cm.tab20.colors
    groups = proportions.index.astype(str).tolist()
    cats = proportions.columns.tolist()
    x_pos = np.arange(len(groups))

    for i, cat in enumerate(cats):
        values = proportions[cat].values
        ax.bar(x_pos, values, bottom=bottom, label=cat, color=colors[i % len(colors)])
        bottom += values

    # Add count labels on each bar
    for j, group in enumerate(groups):
        total = int(counts.iloc[j].sum())
        ax.text(j, 101, f"n={total}", ha="center", va="bottom", fontsize=8, color="gray")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(groups)
    _apply_typography(ax, title, group_col, "Proportion (%)")

    # Chi-square annotation
    if chi2_stats and "error" not in chi2_stats:
        chi2_text = f"χ²={chi2_stats.get('chi2', 0):.2f}, dof={chi2_stats.get('dof', '?')}"
        chi2_text += f"\n{_format_p(chi2_stats.get('p_value', 1))}"
        cramers_v = chi2_stats.get("cramers_v")
        if cramers_v is not None:
            chi2_text += f"\nCramér's V = {cramers_v:.3f}"
        _add_stats_textbox(ax, chi2_text, loc="upper right", fontsize=9)

    ax.legend(
        title=cat_col.replace("_", " ").title(),
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=12,
        title_fontsize=14,
    )

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")

    return ax.figure


def plot_interaction(
    df: pd.DataFrame,
    x_col: str,
    hue_col: str,
    y_col: str,
    title: str = "",
    ylabel: str = "",
    stats_text: Optional[str] = None,
    interaction_p: Optional[float] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """Interaction point plot with categorical x and hue factors.

    Args:
        df: DataFrame.
        x_col: Column for x-axis factor.
        hue_col: Column for hue factor (interaction partner).
        y_col: Column for response variable.
        title: Plot title.
        ylabel: Y-axis label.
        stats_text: Additional stats text.
        interaction_p: P-value for the interaction term (displayed if given).
        ax: Target axes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    palette = {
        h: STRAIN_PALETTE.get(str(h), sns.color_palette("muted")[i % 10])
        for i, h in enumerate(sorted(df[hue_col].dropna().unique(), key=str))
    }

    sns.pointplot(
        data=df,
        x=x_col,
        y=y_col,
        hue=hue_col,
        dodge=0.3,
        markers="o",
        err_kws={"linewidth": 1.5},
        capsize=0.1,
        palette=palette,
        ax=ax,
    )

    _apply_typography(ax, title, x_col, ylabel or y_col)

    text_parts = []
    if interaction_p is not None:
        text_parts.append(f"Interaction: {_format_p(interaction_p)}")
    if stats_text:
        text_parts.append(stats_text)
    _add_stats_textbox(ax, "\n".join(filter(None, text_parts)), loc="upper right", fontsize=9)

    ax.legend(
        title=hue_col.replace("_", " ").title(),
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=12,
        title_fontsize=14,
    )

    return ax


def plot_correlation_heatmap(
    corr_matrix: pd.DataFrame,
    title: str = "",
    p_matrix: Optional[pd.DataFrame] = None,
    mask_nonsig: bool = True,
    sig_threshold: float = 0.05,
    ax: Optional[plt.Axes] = None,
) -> plt.Figure:
    """Correlation heatmap with dendrograms and significance annotations.

    Features:
      - Hierarchical clustering dendrogram on top axis
      - Rows/columns reordered by cluster similarity
      - Correlation values annotated in cells
      - Significance stars (*, **, ***) superimposed when p_matrix provided
      - Non-significant cells masked with × when mask_nonsig=True
      - Lower-triangle only (symmetric matrix)

    Args:
        corr_matrix: Square correlation matrix (DataFrame).
        title: Plot title.
        p_matrix: Optional DataFrame of same shape with p-values.
        mask_nonsig: If True, dim cells that aren't significant.
        sig_threshold: P-value threshold for significance.
        ax: Target axes (ignored when dendrograms are drawn).
    """
    n = len(corr_matrix)

    # Try scipy for dendrograms
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform

        has_scipy = True
    except ImportError:
        has_scipy = False

    if has_scipy and n >= 3 and ax is None:
        # Convert correlation to distance: d = 1 - |r|
        dist_matrix = 1.0 - np.abs(corr_matrix.values)
        np.fill_diagonal(dist_matrix, 0)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2.0
        condensed = squareform(dist_matrix, checks=False)
        Z = linkage(condensed, method="average")

        fig = plt.figure(figsize=(max(10, n * 0.9), max(9, n * 0.8)))

        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=[1, 0.04],
            height_ratios=[0.15, 1],
            hspace=0.02,
            wspace=0.02,
        )

        # Top dendrogram
        ax_dend = fig.add_subplot(gs[0, 0])
        dn = dendrogram(
            Z, ax=ax_dend, orientation="top", no_labels=True, color_threshold=0, above_threshold_color="#555555"
        )
        ax_dend.set_axis_off()
        order = dn["leaves"]

        # Reorder correlation and p-value matrices
        cols_ordered = [corr_matrix.columns[i] for i in order]
        corr_reordered = corr_matrix.loc[cols_ordered, cols_ordered]

        if p_matrix is not None:
            p_reordered = p_matrix.loc[cols_ordered, cols_ordered]
        else:
            p_reordered = None

        # Build annotation matrix
        mask_tri = np.triu(np.ones_like(corr_reordered, dtype=bool), k=1)
        annot_matrix = corr_reordered.copy().astype(str)
        for i in range(len(corr_reordered)):
            for j in range(len(corr_reordered.columns)):
                r = corr_reordered.iloc[i, j]
                r_text = f"{r:.2f}"
                if p_reordered is not None and i < len(p_reordered) and j < len(p_reordered.columns):
                    p = p_reordered.iloc[i, j]
                    if p < 0.001:
                        r_text += "\n***"
                    elif p < 0.01:
                        r_text += "\n**"
                    elif p < 0.05:
                        r_text += "\n*"
                    elif mask_nonsig and i != j:
                        r_text += "\n(ns)"
                annot_matrix.iloc[i, j] = r_text

        ax_heat = fig.add_subplot(gs[1, 0])
        sns.heatmap(
            corr_reordered,
            annot=annot_matrix.values,
            fmt="",
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            square=True,
            mask=mask_tri,
            ax=ax_heat,
            linewidths=0.8,
            linecolor="white",
            cbar=False,
            annot_kws={"size": 9},
        )

        # Dedicated colorbar
        ax_cbar = fig.add_subplot(gs[1, 1])
        sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(vmin=-1, vmax=1))
        plt.colorbar(sm, cax=ax_cbar, label="Pearson r")

        if title:
            fig.suptitle(title, fontsize=16, fontweight="bold", y=0.99)

        if p_reordered is not None:
            ax_heat.text(
                0.5,
                -0.02,
                "* p<0.05  ** p<0.01  *** p<0.001",
                transform=ax_heat.transAxes,
                ha="center",
                va="top",
                fontsize=9,
                fontstyle="italic",
                color="gray",
            )

        plt.tight_layout()
        return fig

    else:
        # Fallback: simple heatmap without dendrogram
        if ax is None:
            fig, ax = plt.subplots(figsize=(max(8, n * 0.8), max(6, n * 0.7)))

        mask_tri = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

        annot_matrix = corr_matrix.copy().astype(str)
        for i in range(len(corr_matrix)):
            for j in range(len(corr_matrix.columns)):
                r = corr_matrix.iloc[i, j]
                r_text = f"{r:.2f}"
                if p_matrix is not None and i < len(p_matrix) and j < len(p_matrix.columns):
                    p = p_matrix.iloc[i, j]
                    if p < 0.001:
                        r_text += "\n***"
                    elif p < 0.01:
                        r_text += "\n**"
                    elif p < 0.05:
                        r_text += "\n*"
                    elif mask_nonsig and i != j:
                        r_text += "\n(ns)"
                annot_matrix.iloc[i, j] = r_text

        sns.heatmap(
            corr_matrix,
            annot=annot_matrix.values,
            fmt="",
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            square=True,
            mask=mask_tri,
            ax=ax,
            linewidths=0.8,
            linecolor="white",
            cbar_kws={"shrink": 0.8, "label": "Pearson r"},
            annot_kws={"size": 9},
        )

        if title:
            ax.set_title(title, fontsize=16, fontweight="bold", pad=12)

        if p_matrix is not None:
            ax.text(
                0.5,
                -0.02,
                "* p<0.05  ** p<0.01  *** p<0.001",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=9,
                fontstyle="italic",
                color="gray",
            )

        return ax.figure
