"""Regional visualization functions for GWAS.

This module provides plots for regional association analysis and gene annotation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core.utils import logging
from metainformant.gwas.visualization.utils import chromosome_label, get_chromosome, standardize_result_columns

logger = logging.get_logger(__name__)


def _on_chromosome(chrom_column: Any, chrom: Any) -> Any:
    """Boolean mask of the rows of ``chrom_column`` that lie on ``chrom``.

    Chromosomes are compared as text (:func:`chromosome_label`), so an integer ``1``, a float ``1.0`` and
    the string ``"1"`` all select the same rows instead of silently selecting none.
    """
    labels = {value: chromosome_label(value) for value in chrom_column.unique()}
    return chrom_column.map(labels) == chromosome_label(chrom)


def _region_label(chrom: Any, start: Optional[int], end: Optional[int]) -> str:
    """Human-readable region such as ``1:500-2,500``; an open end is left blank instead of formatted."""
    label = str(chrom) if chrom is not None else "unknown"
    if start is None and end is None:
        return label
    low = f"{start:,}" if start is not None else ""
    high = f"{end:,}" if end is not None else ""
    return f"{label}:{low}-{high}"


def gene_annotation_plot(
    results: Any,
    output_file: Optional[str | Path] = None,
    *,
    chrom: Optional[str | int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    gene_data: Optional[Dict[str, Any]] = None,
    title: str = "Regional Association Plot",
) -> Dict[str, Any]:
    """Create a regional association plot with gene annotations.

    Args:
        results: GWAS results as list of dicts or DataFrame (columns: CHROM/CHR, POS/BP, p_value/P)
        output_file: Optional output file path
        chrom: Chromosome to plot
        start: Start position (bp)
        end: End position (bp)
        gene_data: Optional dictionary with gene annotations
        title: Plot title

    Returns:
        Dictionary with status and metadata

    Example:
        >>> result = gene_annotation_plot(results, "out.png", chrom="chr1", start=100000, end=200000)
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle  # noqa: F401
    except ImportError:
        logger.warning("matplotlib not available for regional plot")
        return {"status": "skipped", "reason": "matplotlib not available"}

    # Convert list of dicts to internal format
    if isinstance(results, list):
        import pandas as pd

        # Normalize column names: chrom / CHROM / ... become CHR, pos / POS / ... become BP, p_value / ... become P
        results_df = standardize_result_columns(pd.DataFrame(results))
    elif hasattr(results, "columns"):
        results_df = standardize_result_columns(results)
    else:
        logger.error("Input data must be a list of dicts or DataFrame")
        return {"status": "failed", "reason": "Invalid input data format"}

    # Filter data for region if specified
    if chrom is not None:
        mask = _on_chromosome(results_df["CHR"], chrom)
        if start is not None:
            mask = mask & (results_df["BP"] >= start)
        if end is not None:
            mask = mask & (results_df["BP"] <= end)
        region_data = results_df[mask].copy()
    else:
        region_data = results_df.copy()

    if region_data.empty:
        logger.warning(f"No data found for region {chrom}:{start}-{end}")
        return {"status": "success", "reason": "No data in region", "n_variants": 0}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={"height_ratios": [3, 1]})

    # Plot 1: Manhattan plot for region
    region_data["neg_log_p"] = -np.log10(region_data["P"].clip(lower=1e-50))

    ax1.scatter(region_data["BP"], region_data["neg_log_p"], c="blue", s=20, alpha=0.7)
    ax1.axhline(
        y=-np.log10(5e-8),
        color="red",
        linestyle="--",
        alpha=0.8,
        label="Genome-wide significance",
    )
    ax1.set_ylabel("-log10(P-value)", fontsize=12)
    ax1.set_title(title, fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Gene annotation track
    if gene_data and chrom is not None and start is not None and end is not None:
        plot_gene_track(ax2, chrom, start, end, gene_data)
    else:
        # No track can be drawn; say why instead of leaving the panel blank
        track_note = (
            "Gene annotations not provided"
            if not gene_data
            else "Gene annotations cannot be placed without chrom, start and end"
        )
        ax2.text(
            0.5,
            0.5,
            track_note,
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )
        if start is not None and end is not None:
            ax2.set_xlim(start, end)

    chrom_label = chrom or "unknown"
    ax2.set_xlabel(f"Position on {chrom_label} (bp)", fontsize=12)
    ax2.set_yticks([])

    # Format x-axis labels
    def format_bp(x: float, pos: int | None) -> str:
        if x >= 1e6:
            return f"{x / 1e6:.1f}M"
        elif x >= 1e3:
            return f"{x / 1e3:.0f}K"
        else:
            return str(int(x))

    ax1.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    plt.tight_layout()

    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "n_variants": len(region_data),
        "output_path": str(output_file) if output_file else None,
    }


def plot_gene_track(ax: Any, chrom: str | int, start: int, end: int, gene_data: Dict[str, Any]) -> None:
    """Draw the genes of ``gene_data["genes"]`` that overlap ``[start, end]`` on ``ax``.

    Each gene is a dict with integer ``start`` and ``end``, an optional ``name`` and an optional chromosome
    (``chrom``, ``chromosome``, ``chr``, ``CHROM`` or ``CHR``); genes are clipped to the region. A gene that
    names a chromosome is drawn only on that chromosome (compared as text, so ``1`` and ``"1"`` agree); one
    that names none is assumed to lie on ``chrom``. ``gene_data`` is the only annotation source (none is
    bundled), so the track never invents genes: a gene without both coordinates is not drawn (it used to be
    stretched across the whole region), and when no gene can be drawn the track says why instead of staying
    blank.
    """
    from matplotlib.patches import Rectangle

    genes = gene_data.get("genes") or []
    track_chromosome = chromosome_label(chrom)
    y_pos = 0.5
    n_drawn = 0
    n_without_coordinates = 0

    for gene in genes:
        gene_chromosome = get_chromosome(gene)
        if gene_chromosome is not None and chromosome_label(gene_chromosome) != track_chromosome:
            continue  # a gene on another chromosome is not on this track
        gene_start = gene.get("start")
        gene_end = gene.get("end")
        if gene_start is None or gene_end is None:
            n_without_coordinates += 1
            continue
        gene_start = max(start, gene_start)
        gene_end = min(end, gene_end)

        if gene_start < gene_end:
            n_drawn += 1
            # Draw gene rectangle
            width = gene_end - gene_start
            rect = Rectangle(
                (gene_start, y_pos - 0.1),
                width,
                0.2,
                facecolor="lightblue",
                alpha=0.7,
            )
            ax.add_patch(rect)

            # Add gene name if space allows
            if width > (end - start) * 0.05:  # If gene covers >5% of region
                ax.text(
                    (gene_start + gene_end) / 2,
                    y_pos,
                    gene.get("name", ""),
                    ha="center",
                    va="center",
                    fontsize=8,
                    rotation=45,
                )

    if n_without_coordinates:
        logger.warning(f"{n_without_coordinates} gene annotation(s) without start/end coordinates were not drawn")

    if n_drawn == 0:
        if not genes:
            message = "Gene annotations not provided"
        else:
            message = "No annotated genes overlap this region"
            if n_without_coordinates:
                noun = "gene" if n_without_coordinates == 1 else "genes"
                message += f" ({n_without_coordinates} {noun} without coordinates not drawn)"
        ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes)

    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Genes", fontsize=10)


def recombination_rate_plot(
    results: Any,
    output_file: Optional[str | Path] = None,
    *,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    recombination_rates: Optional[Dict[str, List[float]]] = None,
    positions: Optional[Dict[str, List[int]]] = None,
    title: str = "Recombination Rate Profile",
) -> Dict[str, Any]:
    """Create a plot showing recombination rates across genomic regions.

    GWAS results carry no recombination rates and no recombination map is bundled or simulated, so
    the rates must be supplied: pass **both** ``recombination_rates`` and ``positions``. Without them
    nothing is plotted and no file is written; the result has status ``"skipped"`` and
    ``"output_path": None`` rather than a placeholder figure.

    Args:
        results: GWAS results as list of dicts or DataFrame (not used to derive rates; kept for
            call-compatibility)
        output_file: Optional output file path
        chrom: Chromosome to plot (unused without a recombination map)
        start: Start position (bp) (unused without a recombination map)
        end: End position (bp) (unused without a recombination map)
        recombination_rates: Dictionary mapping chromosomes to recombination rate lists (cM/Mb)
        positions: Dictionary mapping chromosomes to position lists (same lengths as the rates)
        title: Plot title

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError:
        logger.warning("matplotlib not available for recombination rate plot")
        return {"status": "skipped", "reason": "matplotlib not available"}

    # If called with explicit recombination data, use the legacy path
    if recombination_rates is not None and positions is not None:
        return _recombination_rate_plot_legacy(recombination_rates, positions, output_file, title)

    # Without an explicit recombination map there is nothing true to draw: report that, write no file
    logger.warning(
        "Recombination rate plot skipped: pass both recombination_rates and positions (a recombination map); "
        "no file was written"
    )
    return {
        "status": "skipped",
        "reason": "Recombination map data not available: both recombination_rates and positions are required",
        "output_path": None,
    }


#: Polynomial order of the Savitzky-Golay line drawn over each recombination-rate series.
_SMOOTHING_POLYORDER = 3
#: Smallest odd window longer than that order (``savgol_filter`` needs ``window_length > polyorder``).
_SMOOTHING_MIN_WINDOW = 5


def _smoothed_rates(rates: List[float]) -> Optional[np.ndarray]:
    """Savitzky-Golay smoothing of one recombination-rate series, or None when no line should be drawn.

    A series of 10 or fewer points is not smoothed, and neither is any series when SciPy is not installed. The
    (odd) window grows with the series, ``len(rates) // 10 * 2 + 1`` capped at 21, but never drops below
    ``_SMOOTHING_MIN_WINDOW``: the plain formula gives 3 for 11 to 19 points, which is not longer than the cubic
    polynomial order, and ``savgol_filter`` raised ``ValueError`` for those chromosomes.
    """
    if len(rates) <= 10:
        return None
    try:
        from scipy import signal
    except ImportError:
        return None
    window_size = max(_SMOOTHING_MIN_WINDOW, min(21, len(rates) // 10 * 2 + 1))
    return np.asarray(signal.savgol_filter(rates, window_size, _SMOOTHING_POLYORDER))


def _recombination_rate_plot_legacy(
    recombination_rates: Dict[str, List[float]],
    positions: Dict[str, List[int]],
    output_file: Optional[str | Path] = None,
    title: str = "Recombination Rate Profile",
) -> Dict[str, Any]:
    """Legacy implementation for recombination rate plotting with explicit data."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {"status": "skipped", "reason": "matplotlib not available"}

    if not recombination_rates or not positions:
        logger.error("No recombination rate or position data provided")
        return {"status": "failed", "reason": "No data provided"}

    # Create figure with subplots for each chromosome
    chromosomes = list(recombination_rates.keys())
    n_chroms = len(chromosomes)

    if n_chroms <= 3:
        n_cols = n_chroms
        n_rows = 1
    else:
        n_cols = 3
        n_rows = (n_chroms + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False, sharey=True)
    fig.suptitle(title, fontsize=16, y=0.95)

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    for i, chrom in enumerate(chromosomes):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]

        rates = recombination_rates[chrom]
        pos = positions[chrom]

        if len(rates) != len(pos):
            logger.warning(f"Mismatched data lengths for {chrom}, skipping")
            continue

        # Plot recombination rate
        ax.plot(pos, rates, color=colors[i % len(colors)], linewidth=1.5, alpha=0.8)

        # Add smoothing line
        smoothed = _smoothed_rates(rates)
        if smoothed is not None:
            ax.plot(pos, smoothed, color="red", linewidth=2, alpha=0.6, label="Smoothed")

        ax.set_title(f"Chromosome {chrom}", fontsize=12)
        ax.set_xlabel("Position (bp)", fontsize=10)
        if col == 0:
            ax.set_ylabel("Recombination Rate (cM/Mb)", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend()

    # Hide unused subplots
    for i in range(n_chroms, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)

    plt.tight_layout()

    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved recombination rate plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "n_chromosomes": n_chroms,
        "output_path": str(output_file) if output_file else None,
    }


def effect_direction_plot(
    effect_sizes: List[float],
    positions: List[int],
    chrom: str,
    output_file: Optional[str | Path] = None,
    title: str = "Effect Direction Plot",
) -> Optional[Any]:
    """Create a plot showing effect directions (positive/negative) across a region.

    Args:
        effect_sizes: List of effect sizes (can be positive or negative)
        positions: List of genomic positions
        chrom: Chromosome name
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for effect direction plot")
        return None

    if not effect_sizes or not positions:
        logger.error("No effect size or position data provided")
        return None

    if len(effect_sizes) != len(positions):
        logger.error("Effect sizes and positions must have same length")
        return None

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f"{title} - {chrom}", fontsize=14)

    # Plot 1: Effect sizes with color coding for direction
    colors = ["red" if x < 0 else "blue" for x in effect_sizes]
    ax1.scatter(positions, effect_sizes, c=colors, alpha=0.7, s=20)
    ax1.axhline(y=0, color="black", linestyle="--", alpha=0.5)
    ax1.set_ylabel("Effect Size", fontsize=12)
    ax1.set_title("Effect Sizes by Direction", fontsize=12)
    ax1.grid(True, alpha=0.3)

    # Add legend
    ax1.scatter([], [], c="red", alpha=0.7, s=20, label="Negative effect")
    ax1.scatter([], [], c="blue", alpha=0.7, s=20, label="Positive effect")
    ax1.legend()

    # Plot 2: Effect size distribution
    ax2.hist(effect_sizes, bins=30, alpha=0.7, color="skyblue", edgecolor="navy")
    ax2.axvline(x=0, color="red", linestyle="--", linewidth=2, label="No effect")
    ax2.axvline(
        x=np.mean(effect_sizes),
        color="green",
        linestyle="-",
        linewidth=2,
        label=f"Mean: {np.mean(effect_sizes):.3f}",
    )
    ax2.set_xlabel("Effect Size", fontsize=12)
    ax2.set_ylabel("Frequency", fontsize=12)
    ax2.set_title("Effect Size Distribution", fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Add statistics text
    positive_count = sum(1 for x in effect_sizes if x > 0)
    negative_count = sum(1 for x in effect_sizes if x < 0)
    positive_percent = positive_count / len(effect_sizes) * 100
    negative_percent = negative_count / len(effect_sizes) * 100
    stats_text = f"""Statistics:
Mean: {np.mean(effect_sizes):.3f}
Median: {np.median(effect_sizes):.3f}
SD: {np.std(effect_sizes):.3f}
Positive effects: {positive_count} ({positive_percent:.1f}%)
Negative effects: {negative_count} ({negative_percent:.1f}%)"""

    fig.text(
        0.02,
        0.02,
        stats_text,
        fontsize=10,
        verticalalignment="bottom",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved effect direction plot to {output_file}")

    return fig


def regional_plot(
    results: Any,
    output_file: Optional[str | Path] = None,
    *,
    chrom: Optional[str | int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    figsize: tuple[int, int] = (12, 6),
    title: Optional[str] = None,
    significance_threshold: float = 5e-8,
    recombination_rate: Any = None,
) -> Dict[str, Any]:
    """Create a regional association plot.

    Args:
        results: GWAS results as list of dicts or DataFrame (columns: CHROM/CHR, POS/BP, p_value/P)
        output_file: Optional output file path
        chrom: Chromosome to plot; compared as text, so ``1``, ``1.0`` and ``"1"`` all select the same rows
        start: Start position (bp); omit for an open lower end
        end: End position (bp); omit for an open upper end
        figsize: Figure size (width, height)
        title: Plot title (auto-generated if None)
        significance_threshold: P-value threshold for significance line
        recombination_rate: Optional recombination rate DataFrame with ``BP`` and ``RATE`` (cM/Mb) columns
            and, optionally, ``CHR``; only rows on ``chrom`` and inside the region are drawn. No track is
            drawn (and none is invented) when it is not supplied or has no rows in the region.

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for regional plot")
        return {"status": "skipped", "reason": "matplotlib not available"}

    # Convert list of dicts to DataFrame
    if isinstance(results, list):
        import pandas as pd

        # Normalize column names: chrom / CHROM / ... become CHR, pos / POS / ... become BP, p_value / ... become P
        results_df = standardize_result_columns(pd.DataFrame(results))
    elif hasattr(results, "columns"):
        results_df = standardize_result_columns(results)
    else:
        return {"status": "failed", "reason": "Invalid input data format"}

    # Filter data for the specified region
    if chrom is not None:
        mask = _on_chromosome(results_df["CHR"], chrom)
        if start is not None:
            mask = mask & (results_df["BP"] >= start)
        if end is not None:
            mask = mask & (results_df["BP"] <= end)
        region_data = results_df[mask].copy()
    else:
        region_data = results_df.copy()

    if region_data.empty:
        chrom_label = chrom or "unknown"
        logger.warning(f"No data found for chr {chrom_label}:{start}-{end}")
        return {"status": "success", "reason": "No data in region", "n_variants": 0}

    # Convert p-values to -log10 scale
    region_data["NEG_LOG_P"] = -np.log10(region_data["P"].clip(lower=1e-50))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={"height_ratios": [3, 1]})

    # Main plot: -log10(p) vs position
    ax1.scatter(
        region_data["BP"],
        region_data["NEG_LOG_P"],
        c=region_data["NEG_LOG_P"],
        cmap="Reds",
        alpha=0.7,
        s=20,
    )

    # Add significance threshold line
    ax1.axhline(
        y=-np.log10(significance_threshold),
        color="red",
        linestyle="--",
        alpha=0.7,
        label=f"p = {significance_threshold}",
    )

    # Add recombination rate if provided (an open end of the region leaves that side unbounded)
    if recombination_rate is not None and hasattr(recombination_rate, "columns"):
        recomb_data = recombination_rate
        if chrom is not None and "CHR" in recomb_data.columns:
            recomb_data = recomb_data[_on_chromosome(recomb_data["CHR"], chrom)]
        if start is not None:
            recomb_data = recomb_data[recomb_data["BP"] >= start]
        if end is not None:
            recomb_data = recomb_data[recomb_data["BP"] <= end]
        if recomb_data.empty:
            logger.warning("No recombination-rate rows fall inside the plotted region; no recombination track drawn")
        else:
            ax2_twin = ax2.twinx()
            ax2_twin.plot(
                recomb_data["BP"],
                recomb_data["RATE"],
                color="blue",
                alpha=0.7,
                linewidth=2,
            )
            ax2_twin.set_ylabel("Recombination Rate (cM/Mb)", color="blue")
            ax2_twin.tick_params(axis="y", labelcolor="blue")

    # Format axes
    if start is not None and end is not None:
        ax1.set_xlim(start, end)
    ax1.set_ylabel("-log10(p)")
    ax1.set_title(title or f"GWAS Regional Plot - Chr {_region_label(chrom, start, end)}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Bottom plot: SNP density
    ax2.hist(region_data["BP"], bins=50, alpha=0.7, color="lightblue", edgecolor="navy")
    ax2.set_xlabel("Position (bp)")
    ax2.set_ylabel("SNP Density")
    if start is not None and end is not None:
        ax2.set_xlim(start, end)
    ax2.grid(True, alpha=0.3)

    # Format x-axis labels
    def format_bp(x: float, pos: int | None) -> str:
        if x >= 1e6:
            return f"{x / 1e6:.1f}M"
        elif x >= 1e3:
            return f"{x / 1e3:.0f}K"
        else:
            return f"{x:.0f}"

    ax2.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    plt.tight_layout()

    if output_file:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional plot to {output_file}")

    plt.close(fig)

    return {
        "status": "success",
        "n_variants": len(region_data),
        "output_path": str(output_file) if output_file else None,
    }


def regional_ld_plot(
    vcf_path: str | Path,
    output_file: Optional[str | Path] = None,
    chrom: str = "chr1",
    start: int = 1000,
    end: int = 2000,
    lead_snp_pos: Optional[int] = None,
    ld_threshold: float = 0.8,
    figsize: tuple[int, int] = (12, 8),
) -> Dict[str, Any]:
    """Create a regional LD plot showing linkage disequilibrium around a locus.

    Args:
        vcf_path: Path to VCF file
        output_file: Optional output file path
        chrom: Chromosome to plot
        start: Start position (bp)
        end: End position (bp)
        lead_snp_pos: Position of lead SNP to highlight
        ld_threshold: LD threshold for visualization
        figsize: Figure size (width, height)

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for regional LD plot")
        return {"status": "skipped", "reason": "matplotlib not available"}

    # Check if VCF file exists
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        logger.warning(f"VCF file not found: {vcf_path}")
        return {"status": "skipped", "reason": "VCF file not found"}

    # Parse VCF file to extract SNP positions and genotypes for LD calculation
    logger.info(f"Parsing VCF file for LD calculation: {vcf_path}")

    try:
        positions = []
        genotypes = []

        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue

                vcf_chrom = fields[0]
                vcf_pos = int(fields[1])

                # Filter to region
                if str(vcf_chrom) != str(chrom):
                    continue
                if not (start <= vcf_pos <= end):
                    continue

                positions.append(vcf_pos)

                # Parse genotypes (simplified: count alt alleles)
                sample_genotypes = []
                for gt_field in fields[9:]:
                    gt = gt_field.split(":")[0]
                    if "/" in gt:
                        alleles = gt.split("/")
                    elif "|" in gt:
                        alleles = gt.split("|")
                    else:
                        alleles = [".", "."]

                    # Count alt alleles (0, 1, or 2)
                    alt_count = sum(1 for a in alleles if a not in [".", "0"])
                    sample_genotypes.append(alt_count)

                genotypes.append(sample_genotypes)

        if len(positions) < 2:
            logger.warning(f"Insufficient SNPs in region {chrom}:{start}-{end} for LD calculation")
            return {
                "status": "skipped",
                "reason": f"Only {len(positions)} SNPs found in region (need at least 2)",
            }

        n_snps = len(positions)
        positions_arr = np.array(positions)
        genotype_arr = np.array(genotypes, dtype=float)

        # Calculate LD matrix (r² values) using correlation
        # Handle missing data by replacing with mean
        for i in range(genotype_arr.shape[0]):
            col = genotype_arr[i]
            valid_mask = ~np.isnan(col)
            if valid_mask.sum() > 0:
                col[~valid_mask] = np.mean(col[valid_mask])

        # Calculate pairwise correlations (r²)
        ld_matrix = np.corrcoef(genotype_arr)
        ld_matrix = np.square(ld_matrix)  # Convert r to r²

        # Handle NaN values (e.g., constant genotypes)
        ld_matrix = np.nan_to_num(ld_matrix, nan=0.0)

        # Ensure diagonal is 1.0
        np.fill_diagonal(ld_matrix, 1.0)

    except Exception as e:
        logger.error(f"Failed to parse VCF for LD calculation: {e}")
        return {"status": "error", "reason": f"VCF parsing error: {e}"}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={"height_ratios": [1, 3]})

    # Top plot: SNP positions
    ax1.scatter(positions_arr, [1] * n_snps, c="blue", s=50, alpha=0.7)
    if lead_snp_pos:
        # Highlight lead SNP
        lead_idx = np.argmin(np.abs(positions_arr - lead_snp_pos))
        ax1.scatter(
            positions_arr[lead_idx],
            1,
            c="red",
            s=100,
            marker="*",
            label=f"Lead SNP (pos {lead_snp_pos})",
        )

    ax1.set_xlim(start, end)
    ax1.set_ylim(0.5, 1.5)
    ax1.set_xlabel("Position (bp)")
    ax1.set_title(f"Regional LD Plot - {chrom}:{start:,}-{end:,}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Format x-axis
    def format_bp(x: float, pos: int | None) -> str:
        if x >= 1e6:
            return f"{x / 1e6:.1f}M"
        elif x >= 1e3:
            return f"{x / 1e3:.0f}K"
        else:
            return f"{x:.0f}"

    ax1.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    # Bottom plot: LD heatmap
    im = ax2.imshow(
        ld_matrix,
        cmap="RdYlBu_r",
        aspect="auto",
        extent=[start, end, end, start],
        vmin=0,
        vmax=1,
    )

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
    cbar.set_label("LD (r²)")

    # Highlight regions above LD threshold
    if ld_threshold > 0:
        high_ld_mask = ld_matrix >= ld_threshold
        if np.any(high_ld_mask):
            # Add contour lines for high LD regions
            ax2.contour(
                ld_matrix,
                levels=[ld_threshold],
                colors="red",
                linewidths=2,
                extent=[start, end, end, start],
            )

    ax2.set_xlabel("Position 1 (bp)")
    ax2.set_ylabel("Position 2 (bp)")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved regional LD plot to {output_file}")

    return {
        "status": "success",
        "n_snps": int(n_snps),
        "chromosome": chrom,
        "region_start": start,
        "region_end": end,
        "lead_snp_pos": lead_snp_pos,
        "ld_threshold": ld_threshold,
        "output_path": str(output_file) if output_file else None,
    }
