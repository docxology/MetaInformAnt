"""Atlas-style visualization module.

Implements the atlas figure grammar from the BeetleAtlas ontogenetic x tissue
transcriptomic atlas (Leader+ 2024, JMB 10.1016/j.jmb.2024.168520) and the
hoverfly stage/tissue transcriptome resource (Yuan+ 2026, Sci Data
10.1038/s41597-026-07148-9), adapted to cross-species tissue-specificity
displays:

- (a) species x tissue tau heatmap with deterministic (dendrogram-free) fixed
  ordering,
- (b) per-orthogroup cross-species expression-profile small multiples,
- (c) tau-vs-orthology-class summary strip plots with descriptive quantiles
  (no confidence intervals).

All statistics rendered here are DESCRIPTIVE ONLY. No inferential tests
(p-values, significance claims) are computed or displayed by this module;
those remain gated behind the campaign evidence-manifest freeze.

Conventions: matplotlib with the Agg backend, Okabe-Ito color-blind-safe
palette, hatching + text redundancy where fills encode categories, explicit
zero marks on logarithmic axes, and annotations placed in log space where a
log axis is in effect.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Okabe-Ito color-blind-safe palette (fixed order; deterministic).
OKABE_ITO: tuple[str, ...] = (
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#000000",  # black
)

# Hatching patterns paired with Okabe-Ito fills so category identity survives
# grayscale printing and color-vision deficiency (hatch + text redundancy).
ORTHOLOGY_CLASS_STYLE: dict[str, tuple[str, str]] = {
    # class label -> (color, hatch)
    "one_to_one": (OKABE_ITO[4], ""),
    "one_to_many": (OKABE_ITO[0], "//"),
    "many_to_many": (OKABE_ITO[5], "xx"),
    "many_to_one": (OKABE_ITO[2], "\\\\"),
}

_DEFAULT_FALLBACK_STYLE: tuple[str, str] = (OKABE_ITO[6], "..")

_TAU_VMIN = 0.0
_TAU_VMAX = 1.0


def compute_tau(expression: pd.DataFrame) -> pd.Series:
    """Compute tissue specificity index tau (Yanai 2005) per gene.

    Parameters
    ----------
    expression:
        Genes x tissues matrix of non-negative expression values
        (e.g. TPM). Rows (genes) with a non-positive mean are dropped,
        matching the lowest-mean-expression filtering convention used
        before tau computation in the tissue-specificity literature.

    Returns
    -------
    pd.Series
        tau in [0, 1] per retained gene, indexed by gene.
    """
    if expression.isna().any().any():
        raise ValueError("expression matrix contains NaN values")
    if (expression < 0).any().any():
        raise ValueError("expression values must be non-negative")
    means = expression.mean(axis=1)
    retained = means > 0
    if not retained.any():
        raise ValueError("all rows have non-positive mean expression")
    x = expression.loc[retained]
    x_hat = x.div(x.max(axis=1), axis=0)
    n = x.shape[1]
    tau = (
        ((1.0 - x_hat).sum(axis=1)) / (n - 1)
        if n > 1
        else pd.Series(0.0, index=x.index)
    )
    tau.name = "tau"
    return tau.astype(float)


def _validate_tau_frame(tau: pd.DataFrame) -> None:
    if tau.empty:
        raise ValueError("tau frame is empty")
    if not np.isfinite(tau.to_numpy(dtype=float)).all():
        raise ValueError("tau frame contains non-finite values")
    if (tau < 0).any().any() or (tau > 1).any().any():
        raise ValueError("tau values must lie in [0, 1]")


def _resolve_output_path(output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _sorted_labels(labels: Iterable) -> list[str]:
    """Deterministic, dendrogram-free ordering: plain lexicographic sort."""
    return sorted(labels, key=str)


def plot_tau_heatmap(
    tau: pd.DataFrame,
    output_path: str | Path,
    *,
    title: str = "Tissue specificity (tau) by species and tissue",
    annot: bool = True,
    dpi: int = 200,
) -> Path:
    """Render a species x tissue tau heatmap with deterministic fixed ordering.

    Parameters
    ----------
    tau:
        Species (index) x tissues (columns) matrix of mean tau values in
        [0, 1]. Ordering is a plain lexicographic sort of both axes — no
        clustering, so the same input frame always yields the same layout.
    """
    _validate_tau_frame(tau)
    ordered = tau.reindex(
        index=_sorted_labels(tau.index), columns=_sorted_labels(tau.columns)
    )
    path = _resolve_output_path(output_path)

    fig, ax = plt.subplots(
        figsize=(
            1.1 * max(len(ordered.columns), 4) + 2,
            0.6 * max(len(ordered.index), 4) + 2,
        )
    )
    im = ax.imshow(
        ordered.to_numpy(dtype=float),
        cmap="cividis",
        vmin=_TAU_VMIN,
        vmax=_TAU_VMAX,
        aspect="auto",
    )
    ax.set_xticks(range(len(ordered.columns)))
    ax.set_xticklabels(ordered.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(ordered.index)))
    ax.set_yticklabels(ordered.index)
    ax.set_xlabel("Tissue")
    ax.set_ylabel("Species")
    ax.set_title(title)

    if annot:
        data = ordered.to_numpy(dtype=float)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                # Text redundancy: every cell is annotated regardless of fill,
                # with a contrast-switching color across the cividis midpoint.
                color = "#FFFFFF" if data[i, j] < 0.5 else "#111111"
                ax.text(
                    j,
                    i,
                    f"{data[i, j]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=color,
                )

    cbar = fig.colorbar(im, ax=ax, ticks=np.linspace(_TAU_VMIN, _TAU_VMAX, 6))
    cbar.set_label("tau (0 = ubiquitous, 1 = tissue-specific)")

    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    logger.info("wrote tau heatmap to %s", path)
    return path


def plot_orthogroup_small_multiples(
    expression: pd.DataFrame,
    output_path: str | Path,
    *,
    tissues: Sequence[str] | None = None,
    max_panels: int | None = None,
    log_y: bool = True,
    dpi: int = 200,
) -> Path:
    """Render per-orthogroup cross-species expression-profile small multiples.

    Parameters
    ----------
    expression:
        Long-format frame with columns ``orthogroup``, ``species``,
        ``tissue``, ``expression`` (non-negative; log axis recommended).
    max_panels:
        If given, render only the first ``max_panels`` orthogroups in
        lexicographic order (deterministic selection).
    """
    required = {"orthogroup", "species", "tissue", "expression"}
    missing = required - set(expression.columns)
    if missing:
        raise ValueError(f"expression frame missing columns: {sorted(missing)}")
    if expression.empty:
        raise ValueError("expression frame is empty")
    if (expression["expression"] < 0).any() or not np.isfinite(
        expression["expression"].to_numpy(dtype=float)
    ).all():
        raise ValueError("expression values must be finite and non-negative")

    path = _resolve_output_path(output_path)
    orthogroups = _sorted_labels(expression["orthogroup"].unique())
    if max_panels is not None:
        orthogroups = orthogroups[:max_panels]
    if not orthogroups:
        raise ValueError("no orthogroups to render")

    species_order = _sorted_labels(expression["species"].unique())
    color_map = {
        sp: OKABE_ITO[i % len(OKABE_ITO)] for i, sp in enumerate(species_order)
    }
    tissue_order = (
        list(tissues)
        if tissues is not None
        else _sorted_labels(expression["tissue"].unique())
    )

    n = len(orthogroups)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False, sharex=True
    )
    x_pos = np.arange(len(tissue_order))

    for idx, og in enumerate(orthogroups):
        ax = axes[idx // ncols][idx % ncols]
        sub = expression[expression["orthogroup"] == og]
        for sp in species_order:
            sp_sub = sub[sub["species"] == sp].set_index("tissue")["expression"]
            values = [float(sp_sub.get(t, np.nan)) for t in tissue_order]
            # Redundant encoding: distinct color per species plus a distinct
            # marker, so profiles remain distinguishable in grayscale.
            marker = "o+sD^v<>"[species_order.index(sp) % 8]
            ax.plot(
                x_pos,
                values,
                color=color_map[sp],
                marker=marker,
                linestyle="-",
                linewidth=1.2,
                label=sp,
            )
        ax.set_title(str(og), fontsize=9)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(tissue_order, rotation=45, ha="right", fontsize=7)
        if log_y:
            ax.set_yscale("log")
            # Explicit zero mark: log axes silently drop true zeros, so draw a
            # reference line at the panel's minimum value as the floor.
            floor = float(np.nanmin(sub["expression"].to_numpy(dtype=float)))
            ax.axhline(max(floor, 1e-6), color="#888888", linewidth=0.6)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    handles = [
        plt.Line2D([], [], color=color_map[sp], marker="o+sD^v<>"[i % 8], linestyle="-")
        for i, sp in enumerate(species_order)
    ]
    fig.legend(
        handles,
        species_order,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.97),
        ncol=min(len(species_order), 4),
        fontsize=8,
    )
    fig.suptitle(
        "Cross-species expression profiles per orthogroup",
        y=1.0 - 0.02 * (1 / (nrows + 1)),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote orthogroup small multiples to %s", path)
    return path


def plot_tau_orthology_strips(
    tau_classes: pd.DataFrame,
    output_path: str | Path,
    *,
    tau_column: str = "tau",
    class_column: str = "orthology_class",
    dpi: int = 200,
) -> Path:
    """Render tau strip plots grouped by orthology class with descriptive quantiles.

    Displays per-gene tau as jittered points plus median / Q1 / Q3 quantile
    bars per class. Descriptive statistics only — no confidence intervals,
    no hypothesis tests, no significance annotation.

    Parameters
    ----------
    tau_classes:
        Frame with a tau column and an orthology-class column. Classes are
        rendered in lexicographic order with Okabe-Ito fill + hatch styling
        (hatch + text redundancy).
    """
    for column in (tau_column, class_column):
        if column not in tau_classes.columns:
            raise ValueError(f"frame missing required column: {column}")
    if tau_classes.empty:
        raise ValueError("frame is empty")
    values = tau_classes[tau_column].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("tau column contains non-finite values")
    if (values < 0).any() or (values > 1).any():
        raise ValueError("tau values must lie in [0, 1]")

    path = _resolve_output_path(output_path)
    classes = _sorted_labels(tau_classes[class_column].unique())
    fig, ax = plt.subplots(figsize=(1.8 * len(classes) + 2, 4.5))

    rng_seed = 0  # fixed seed: jitter is deterministic for a given input frame
    rng = np.random.default_rng(rng_seed)

    for i, cls in enumerate(classes):
        sub = tau_classes[tau_classes[class_column] == cls]
        y = sub[tau_column].to_numpy(dtype=float)
        color, hatch = ORTHOLOGY_CLASS_STYLE.get(str(cls), _DEFAULT_FALLBACK_STYLE)
        jitter = rng.uniform(-0.12, 0.12, size=len(y))
        ax.scatter(
            i + jitter, y, s=10, color=color, alpha=0.5, edgecolors="none", zorder=1
        )
        q1, median, q3 = (float(np.quantile(y, q)) for q in (0.25, 0.5, 0.75))
        ax.hlines([q1, q3], i - 0.22, i + 0.22, color=color, linewidth=2.0, zorder=3)
        ax.hlines([median], i - 0.3, i + 0.3, color="#111111", linewidth=2.0, zorder=4)
        # Text redundancy: print the median value next to each class strip.
        ax.text(
            i + 0.34, median, f"{median:.2f}", va="center", fontsize=7, color="#111111"
        )

    # Legend swatches carry both fill and hatch so the class key survives
    # grayscale reproduction.
    for cls in classes:
        color, hatch = ORTHOLOGY_CLASS_STYLE.get(str(cls), _DEFAULT_FALLBACK_STYLE)
        ax.scatter(
            [],
            [],
            s=40,
            color=color,
            hatch=hatch,
            label=str(cls),
            edgecolors="#333333",
            linewidths=0.4,
        )
    ax.legend(loc="lower right", fontsize=8, title="Orthology class", title_fontsize=8)

    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=20, ha="right")
    ax.set_ylabel("tau")
    ax.set_ylim(-0.02, 1.02)
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.set_title(
        "Tissue specificity by orthology class\n(descriptive quantiles; no inferential statistics)"
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    logger.info("wrote tau/orthology strip plot to %s", path)
    return path
