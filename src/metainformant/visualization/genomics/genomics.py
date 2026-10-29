"""Genomic data visualization functions for GWAS and genomic analysis.

This module provides specialized plotting functions for genomic data including
Manhattan plots, volcano plots, regional association plots, chromosome ideograms,
coverage plots, and variant visualizations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Hashable, Mapping, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

#: GRCh38 primary-assembly chromosome lengths in base pairs, used by :func:`chromosome_ideogram` by default.
GRCH38_CHROMOSOME_LENGTHS: dict[str, int] = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
}


def _chromosome_sort_key(label: Hashable) -> tuple[int, int, str]:
    """Natural genome order for a chromosome label: 1, 2, ..., 10, ..., X, Y, M/MT, then anything else by name.

    A ``chr`` prefix is ignored, so ``"chr10"`` follows ``"chr2"``. Plain string sorting put "10" and "11" before
    "2" and raised ``TypeError`` for a mix of integer and string labels.
    """
    name = str(label).strip()
    bare = name[3:] if name.lower().startswith("chr") else name
    if bare.isdigit():
        return (0, int(bare), "")
    special = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if bare.upper() in special:
        return (1, special[bare.upper()], "")
    return (2, 0, name)


def _check_pvalues(values: pd.Series, column: str) -> None:
    """Raise ``ValueError`` for a p-value outside [0, 1]; ``-log10`` of one is a fabricated peak or a negative height.

    A p-value below 0 used to be clipped to 1e-300 and drawn at a height of 300, and one above 1 drawn below zero.
    Missing values (NaN) are not drawn and are not an error.
    """
    invalid = values[(values < 0) | (values > 1)]
    if len(invalid):
        raise ValueError(
            f"p-values in column {column!r} must lie between 0 and 1; {len(invalid)} value(s) do not "
            f"(for example {invalid.iloc[0]!r})"
        )


def manhattan_plot(
    data: pd.DataFrame,
    *,
    chr_col: str = "CHR",
    pos_col: str = "BP",
    pval_col: str = "P",
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a Manhattan plot for GWAS results.

    Args:
        data: DataFrame with chromosome, position, and p-value columns
        chr_col: Column name for chromosome
        pos_col: Column name for position
        pval_col: Column name for p-value
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing or data is invalid
    """
    validation.validate_type(data, pd.DataFrame, "data")

    required_cols = [chr_col, pos_col, pval_col]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Prepare data
    if data.empty:
        raise ValueError("Manhattan plot requires at least one data row")
    _check_pvalues(data[pval_col], pval_col)
    plot_data = data.copy()
    plot_data["logP"] = -np.log10(plot_data[pval_col].clip(lower=1e-300))  # Avoid log(0)

    # Calculate cumulative positions for chromosomes (natural genome order: 1, 2, ..., 10, ..., X, Y)
    chromosomes = sorted(plot_data[chr_col].unique(), key=_chromosome_sort_key)
    cumulative_pos = 0

    for chr_num in chromosomes:
        chr_mask = plot_data[chr_col] == chr_num
        plot_data.loc[chr_mask, "cumulative_pos"] = plot_data.loc[chr_mask, pos_col] + cumulative_pos
        cumulative_pos += plot_data.loc[chr_mask, pos_col].max() + 1e7  # Add gap between chromosomes

    # Plot points colored by chromosome
    colors = plt.cm.tab10(np.linspace(0, 1, len(chromosomes)))
    point_kwargs = {"s": kwargs.pop("s", 1), "alpha": kwargs.pop("alpha", 0.8)}
    for i, chr_num in enumerate(chromosomes):
        chr_mask = plot_data[chr_col] == chr_num
        ax.scatter(
            plot_data.loc[chr_mask, "cumulative_pos"],
            plot_data.loc[chr_mask, "logP"],
            c=[colors[i]],
            **point_kwargs,
            **kwargs,
        )

    # Add chromosome labels at centers
    chr_centers = []
    chr_labels = []
    for chr_num in chromosomes:
        chr_mask = plot_data[chr_col] == chr_num
        if chr_mask.any():
            center_pos = (
                plot_data.loc[chr_mask, "cumulative_pos"].min() + plot_data.loc[chr_mask, "cumulative_pos"].max()
            ) / 2
            chr_centers.append(center_pos)
            chr_labels.append(str(chr_num))

    ax.set_xticks(chr_centers)
    ax.set_xticklabels(chr_labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title("Manhattan Plot")
    ax.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.7, label="Genome-wide significance")

    # Add significance threshold legend if line is visible
    if plot_data["logP"].max() > -np.log10(5e-8):
        ax.legend()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")

    return ax


def volcano_plot(
    data: pd.DataFrame,
    *,
    log2fc_col: str = "log2FoldChange",
    pval_col: str = "padj",
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a volcano plot for differential expression analysis.

    Args:
        data: DataFrame with log fold change and p-value columns
        log2fc_col: Column name for log2 fold change
        pval_col: Column name for adjusted p-value
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(data, pd.DataFrame, "data")

    required_cols = [log2fc_col, pval_col]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

    # Prepare data
    plot_data = data.copy()
    plot_data["logP"] = -np.log10(plot_data[pval_col].clip(lower=1e-300))

    # Color points based on significance and fold change
    colors = []
    for _, row in plot_data.iterrows():
        if row["logP"] >= -np.log10(0.05) and abs(row[log2fc_col]) >= 1:
            colors.append("red")  # Significant and large fold change
        elif row["logP"] >= -np.log10(0.05):
            colors.append("orange")  # Significant
        elif abs(row[log2fc_col]) >= 1:
            colors.append("blue")  # Large fold change
        else:
            colors.append("gray")  # Not significant

    ax.scatter(
        plot_data[log2fc_col],
        plot_data["logP"],
        c=colors,
        s=kwargs.pop("s", 20),
        alpha=kwargs.pop("alpha", 0.6),
        **kwargs,
    )

    # Add threshold lines
    ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", alpha=0.7)
    ax.axvline(x=1, color="black", linestyle="--", alpha=0.7)
    ax.axvline(x=-1, color="black", linestyle="--", alpha=0.7)

    ax.set_xlabel("log₂(Fold Change)")
    ax.set_ylabel("-log₁₀(adjusted p-value)")
    ax.set_title("Volcano Plot")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Volcano plot saved to {output_path}")

    return ax


def regional_plot(
    data: pd.DataFrame,
    *,
    chr: str,
    start: int,
    end: int,
    recombination_rate: pd.DataFrame | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a regional association plot for a specific genomic region.

    Args:
        data: DataFrame with genomic data
        chr: Chromosome name
        start: Start position
        end: End position
        recombination_rate: Optional recombination map with ``BP`` (position) and ``RATE`` (cM/Mb) columns
            and, optionally, ``CHR``. The rows on ``chr`` inside the region are drawn as a line on a
            secondary axis, in the same positions-from-``start`` coordinates as the p-values. Without it no
            recombination track is drawn: no map is bundled and none is simulated.
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If region parameters are invalid, or ``recombination_rate`` lacks ``BP`` or ``RATE``
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if start >= end:
        raise ValueError("Start position must be less than end position")

    if recombination_rate is not None:
        validation.validate_type(recombination_rate, pd.DataFrame, "recombination_rate")
        missing_columns = [col for col in ("BP", "RATE") if col not in recombination_rate.columns]
        if missing_columns:
            raise ValueError(f"recombination_rate is missing required columns: {missing_columns}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Filter data for the specified region
    region_mask = (data["CHR"].astype(str) == str(chr)) & (data["BP"] >= start) & (data["BP"] <= end)
    region_data = data[region_mask].copy()

    if region_data.empty:
        logger.warning(f"No data found in region {chr}:{start}-{end}")
        ax.text(0.5, 0.5, f"No data in region {chr}:{start}-{end}", transform=ax.transAxes, ha="center", va="center")
        ax.set_xlim(start, end)
        ax.set_xlabel(f"Position on chromosome {chr}")
        ax.set_ylabel("-log₁₀(p-value)")
        ax.set_title(f"Regional Plot - {chr}:{start}-{end}")
        # an empty region is still a figure the caller asked to have written
        if output_path:
            paths.ensure_directory(Path(output_path).parent)
            save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
            logger.info(f"Regional plot saved to {output_path}")
        return ax

    _check_pvalues(region_data["P"], "P")

    # Calculate relative positions
    region_data["relative_pos"] = region_data["BP"] - start
    region_data["logP"] = -np.log10(region_data["P"].clip(lower=1e-300))

    # Plot association results
    ax.scatter(
        region_data["relative_pos"],
        region_data["logP"],
        s=kwargs.pop("s", 30),
        alpha=kwargs.pop("alpha", 0.7),
        **kwargs,
    )

    # Recombination rate: drawn from the supplied map only (nothing is simulated), on a secondary axis
    if recombination_rate is not None:
        rate_rows = recombination_rate
        if "CHR" in rate_rows.columns:
            rate_rows = rate_rows[rate_rows["CHR"].astype(str) == str(chr)]
        rate_rows = rate_rows[(rate_rows["BP"] >= start) & (rate_rows["BP"] <= end)].sort_values("BP")
        if rate_rows.empty:
            logger.warning(f"No recombination-rate rows in region {chr}:{start}-{end}; no recombination track drawn")
        else:
            rate_ax = ax.twinx()
            rate_ax.plot(rate_rows["BP"] - start, rate_rows["RATE"], color="tab:blue", alpha=0.7, linewidth=1.5)
            rate_ax.set_ylabel("Recombination rate (cM/Mb)", color="tab:blue")
            rate_ax.tick_params(axis="y", labelcolor="tab:blue")

    ax.set_xlabel(f"Position on chromosome {chr} (bp from {start:,})")
    ax.set_ylabel("-log₁₀(p-value)")
    ax.set_title(f"Regional Association Plot - {chr}:{start:,}-{end:,}")
    ax.axhline(y=-np.log10(5e-8), color="red", linestyle="--", alpha=0.7)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Regional plot saved to {output_path}")

    return ax


def circular_manhattan_plot(
    data: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a circular Manhattan plot.

    Args:
        data: DataFrame with chromosome, position, and p-value columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 8)), subplot_kw={"projection": "polar"})

    # Prepare data similar to linear Manhattan plot
    plot_data = data.copy()
    plot_data["logP"] = -np.log10(plot_data["P"].clip(lower=1e-300))

    # Calculate angles for chromosomes (divide circle evenly), in natural genome order
    chromosomes = sorted(plot_data["CHR"].unique(), key=_chromosome_sort_key)
    chr_angles = np.linspace(0, 2 * np.pi, len(chromosomes) + 1)[:-1]

    # Marker size and opacity are read once: popping them inside the loop gave the caller's values to the
    # first chromosome only and the defaults to every later one
    point_size = kwargs.pop("s", 2)
    point_alpha = kwargs.pop("alpha", 0.7)

    # Plot each chromosome
    for i, (chr_num, angle) in enumerate(zip(chromosomes, chr_angles)):
        chr_mask = plot_data["CHR"] == chr_num
        chr_data = plot_data[chr_mask]

        # Normalize positions within chromosome to angle range
        chr_start = chr_data["BP"].min()
        chr_end = chr_data["BP"].max()
        chr_range = chr_end - chr_start

        sector_width = 2 * np.pi / len(chromosomes)
        if chr_range > 0:
            relative_angles = angle + (chr_data["BP"] - chr_start) / chr_range * sector_width
        else:
            # Every variant of this chromosome is at one position (a lone variant, say): put them in the
            # middle of the chromosome's sector rather than dropping the chromosome from the plot
            relative_angles = np.full(len(chr_data), angle + sector_width / 2)

        ax.scatter(
            relative_angles,
            chr_data["logP"],
            s=point_size,
            alpha=point_alpha,
            label=f"Chr {chr_num}" if i < 5 else "",  # Label first few chromosomes
            **kwargs,
        )

    # Add chromosome labels
    for angle, chr_num in zip(chr_angles, chromosomes):
        ax.text(angle, ax.get_ylim()[1] * 1.1, f"{chr_num}", ha="center", va="bottom", fontsize=8)

    ax.set_title("Circular Manhattan Plot")
    cast("Any", ax).set_rlabel_position(90)
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular Manhattan plot saved to {output_path}")

    return ax


def chromosome_ideogram(
    *,
    chromosome_lengths: Mapping[Any, int] | None = None,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Draw one horizontal bar per chromosome, with its length relative to the longest chromosome.

    Only the lengths are drawn (no cytogenetic bands, which would need real banding data). By default the bars
    are the 24 human chromosomes of GRCh38 (:data:`GRCH38_CHROMOSOME_LENGTHS`); pass ``chromosome_lengths`` for
    any other genome, such as the linkage groups of an ant assembly.

    Args:
        chromosome_lengths: Mapping from chromosome name to length in base pairs, drawn top to bottom in mapping
            order. Defaults to GRCh38.
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib barh() (``alpha``, ``color``, ``edgecolor``, ...);
            ``figsize`` sizes a new figure.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If ``chromosome_lengths`` is empty or has a length that is not positive
    """
    lengths = dict(GRCH38_CHROMOSOME_LENGTHS if chromosome_lengths is None else chromosome_lengths)
    if not lengths:
        raise ValueError("chromosome_lengths cannot be empty")
    if any(length <= 0 for length in lengths.values()):
        raise ValueError("chromosome lengths must be positive")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 4)))

    # Plot chromosomes as horizontal bars (caller keyword arguments override the defaults)
    max_length = max(lengths.values())
    bar_kwargs: dict[str, Any] = {"height": 0.8, "alpha": 0.7, "color": "lightblue", "edgecolor": "black", **kwargs}
    ax.barh(range(len(lengths)), [length / max_length for length in lengths.values()], **bar_kwargs)

    ax.set_yticks(range(len(lengths)))
    ax.set_yticklabels([str(name) for name in lengths])
    ax.set_xlabel("Relative Chromosome Length")
    ax.set_title("Human Chromosome Lengths (GRCh38)" if chromosome_lengths is None else "Chromosome Lengths")
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Chromosome ideogram saved to {output_path}")

    return ax


def coverage_plot(
    coverage: np.ndarray,
    positions: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a sequencing coverage plot.

    Args:
        coverage: Array of coverage depths
        positions: Array of genomic positions
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If arrays have mismatched lengths
    """
    validation.validate_type(coverage, np.ndarray, "coverage")
    validation.validate_type(positions, np.ndarray, "positions")

    if len(coverage) != len(positions):
        raise ValueError("Coverage and positions arrays must have same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 4)))

    line_kwargs = dict(kwargs)
    fill_kwargs = {"alpha": kwargs.pop("alpha", 0.3), **kwargs}
    ax.plot(positions, coverage, **line_kwargs)
    ax.fill_between(positions, coverage, **fill_kwargs)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Coverage Depth")
    ax.set_title("Sequencing Coverage")
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Coverage plot saved to {output_path}")

    return ax


def variant_plot(
    variants: pd.DataFrame, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a variant visualization plot.

    Args:
        variants: DataFrame with variant information
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(variants, pd.DataFrame, "variants")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 6)))

    # Assume variants has POS and some value column
    if "POS" not in variants.columns:
        raise ValueError("Variants DataFrame must contain 'POS' column")

    # Simple variant visualization - plot positions
    positions = variants["POS"].values
    y_values = np.ones(len(positions))  # Default y-position

    # Add deterministic jitter for visibility (fixed seed: reproducible output)
    y_values += np.random.default_rng(0).uniform(-0.1, 0.1, len(positions))

    ax.scatter(positions, y_values, s=kwargs.pop("s", 10), alpha=kwargs.pop("alpha", 0.7), **kwargs)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Variants")
    ax.set_title("Variant Positions")
    ax.set_yticks([])  # Hide y-axis ticks

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Variant plot saved to {output_path}")

    return ax
