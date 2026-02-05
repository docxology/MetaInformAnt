"""Structural variant visualization.

Provides publication-quality plots for SV analysis including Circos-style
genome-wide views, coverage tracks with SV overlays, size distribution
histograms, SV type summaries, breakpoint detail views, and genome-wide
CNV profiles.

All plot functions save to files and return the matplotlib Figure object
for further customization.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.figure import Figure
    from matplotlib.colors import Normalize

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    Figure = Any  # type: ignore[assignment, misc]

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

logger = get_logger(__name__)

# Color scheme for SV types
SV_COLORS: dict[str, str] = {
    "DEL": "#E74C3C",  # Red
    "DUP": "#3498DB",  # Blue
    "INV": "#2ECC71",  # Green
    "TRA": "#9B59B6",  # Purple
    "BND": "#9B59B6",  # Purple
    "INS": "#F39C12",  # Orange
    "UNKNOWN": "#95A5A6",  # Gray
    "NEUTRAL": "#BDC3C7",  # Light gray
    "AMP": "#1A5276",  # Dark blue
    "HOMODEL": "#C0392B",  # Dark red
}

# CNV state colors
CNV_COLORS: dict[str, str] = {
    "HOMODEL": "#C0392B",
    "DEL": "#E74C3C",
    "NEUTRAL": "#BDC3C7",
    "DUP": "#3498DB",
    "AMP": "#1A5276",
}


def _check_matplotlib() -> None:
    """Raise an import error if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for visualization. Install with: uv pip install matplotlib")


def plot_circos(
    variants: list[dict[str, Any]],
    chromosomes: dict[str, int],
    output_path: str | Path,
    title: str = "Structural Variants - Circos View",
    figsize: tuple[float, float] = (12, 12),
    show_labels: bool = True,
) -> Any:
    """Generate a Circos-style genome-wide structural variant plot.

    Draws chromosomes as arcs around a circle, with SVs represented as
    colored chords connecting breakpoint pairs. Intra-chromosomal SVs
    are drawn as arcs within the chromosome. Colors indicate SV type.

    Args:
        variants: List of variant dictionaries with:
            - 'chrom': Chromosome
            - 'start': Start position
            - 'end': End position
            - 'sv_type': SV type string
            - 'chrom2': Second chromosome (for translocations, optional)
        chromosomes: Dictionary mapping chromosome names to sizes (in bp).
        output_path: Path to save the plot image.
        title: Plot title.
        figsize: Figure size in inches.
        show_labels: Whether to show chromosome labels.

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={"projection": "polar"})

    # Calculate total genome size and chromosome angular positions
    total_size = sum(chromosomes.values())
    gap_fraction = 0.02  # 2% gap between chromosomes
    n_chroms = len(chromosomes)
    total_gap = gap_fraction * n_chroms
    available = 2 * math.pi * (1 - total_gap)
    gap_angle = 2 * math.pi * gap_fraction

    chrom_angles: dict[str, tuple[float, float]] = {}
    current_angle = 0.0
    chrom_order = sorted(chromosomes.keys(), key=_chrom_sort_key)

    for chrom in chrom_order:
        size = chromosomes[chrom]
        arc_length = available * size / total_size
        chrom_angles[chrom] = (current_angle, current_angle + arc_length)
        current_angle += arc_length + gap_angle

    # Draw chromosome arcs
    outer_radius = 0.9
    inner_radius = 0.85
    _draw_chromosome_arcs(ax, chrom_angles, outer_radius, inner_radius, chrom_order, show_labels)

    # Draw SV chords
    chord_radius = inner_radius - 0.02

    for v in variants:
        chrom1 = v.get("chrom", "")
        start1 = v.get("start", 0)
        chrom2 = v.get("chrom2", chrom1)
        if not chrom2:
            chrom2 = chrom1
        end2 = v.get("end", 0)
        sv_type = v.get("sv_type", "UNKNOWN")
        if hasattr(sv_type, "value"):
            sv_type = sv_type.value

        if chrom1 not in chrom_angles or chrom2 not in chrom_angles:
            continue

        color = SV_COLORS.get(sv_type, SV_COLORS["UNKNOWN"])

        # Convert genomic positions to angles
        angle1 = _pos_to_angle(chrom1, start1, chrom_angles, chromosomes)
        angle2 = _pos_to_angle(chrom2, end2, chrom_angles, chromosomes)

        if angle1 is None or angle2 is None:
            continue

        # Draw chord (Bezier curve through center region)
        _draw_chord(ax, angle1, angle2, chord_radius, color, alpha=0.5)

    ax.set_title(title, pad=20, fontsize=14, fontweight="bold")
    ax.set_ylim(0, 1.1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)

    # Legend
    legend_patches = [mpatches.Patch(color=SV_COLORS[t], label=t) for t in ["DEL", "DUP", "INV", "TRA", "INS"]]
    ax.legend(
        handles=legend_patches,
        loc="lower right",
        bbox_to_anchor=(1.15, 0),
        fontsize=10,
    )

    plt.tight_layout()
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output), dpi=150, bbox_inches="tight")
    logger.info(f"Circos plot saved to {output}")
    plt.close(fig)
    return fig


def plot_coverage_track(
    coverage: list[float] | Any,
    variants: list[dict[str, Any]],
    region: dict[str, Any] | tuple[str, int, int],
    output_path: str | Path = "",
    title: str = "",
    figsize: tuple[float, float] = (14, 5),
    bin_size: int = 1000,
) -> Any:
    """Plot read depth coverage with structural variant overlay.

    Shows a coverage track (read depth) for a genomic region with
    colored rectangles indicating structural variant positions and types.

    Args:
        coverage: Array of coverage values (one per bin).
        variants: List of variant dictionaries (same format as other functions).
        region: Genomic region as dict with 'chrom', 'start', 'end' or
            tuple (chrom, start, end).
        output_path: Path to save the plot. If empty, plot is not saved.
        title: Plot title. Auto-generated from region if empty.
        figsize: Figure size in inches.
        bin_size: Size of each coverage bin in base pairs.

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    if isinstance(region, tuple):
        chrom, reg_start, reg_end = region
    else:
        chrom = region.get("chrom", "")
        reg_start = region.get("start", 0)
        reg_end = region.get("end", 0)

    if not title:
        title = f"Coverage - {chrom}:{reg_start:,}-{reg_end:,}"

    if np is not None:
        cov = np.asarray(coverage, dtype=np.float64)
    else:
        cov = list(coverage)

    n_bins = len(cov)
    positions = [reg_start + i * bin_size for i in range(n_bins)]

    fig, (ax_cov, ax_sv) = plt.subplots(2, 1, figsize=figsize, height_ratios=[4, 1], sharex=True)

    # Coverage track
    if np is not None:
        median_cov = float(np.median(cov[cov > 0])) if np.any(cov > 0) else 1.0
    else:
        positive = [c for c in cov if c > 0]
        median_cov = sorted(positive)[len(positive) // 2] if positive else 1.0

    ax_cov.fill_between(positions, cov, alpha=0.6, color="#3498DB")
    ax_cov.plot(positions, cov, linewidth=0.5, color="#2C3E50")
    ax_cov.axhline(y=median_cov, color="#E74C3C", linestyle="--", linewidth=1, label=f"Median: {median_cov:.0f}x")
    ax_cov.set_ylabel("Read Depth", fontsize=11)
    ax_cov.set_title(title, fontsize=13, fontweight="bold")
    ax_cov.legend(fontsize=9)

    # SV overlay track
    for v in variants:
        v_chrom = v.get("chrom", "")
        if v_chrom != chrom:
            continue

        v_start = max(v.get("start", 0), reg_start)
        v_end = min(v.get("end", 0), reg_end)
        sv_type = v.get("sv_type", "UNKNOWN")
        if hasattr(sv_type, "value"):
            sv_type = sv_type.value

        if v_start >= reg_end or v_end <= reg_start:
            continue

        color = SV_COLORS.get(sv_type, SV_COLORS["UNKNOWN"])
        width = v_end - v_start

        rect = mpatches.Rectangle(
            (v_start, 0),
            width,
            1,
            linewidth=1,
            edgecolor=color,
            facecolor=color,
            alpha=0.6,
        )
        ax_sv.add_patch(rect)
        ax_sv.text(
            v_start + width / 2,
            0.5,
            sv_type,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )

    ax_sv.set_xlim(reg_start, reg_end)
    ax_sv.set_ylim(0, 1)
    ax_sv.set_ylabel("SVs", fontsize=11)
    ax_sv.set_xlabel(f"Position on {chrom} (bp)", fontsize=11)
    ax_sv.set_yticks([])

    plt.tight_layout()

    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output), dpi=150, bbox_inches="tight")
        logger.info(f"Coverage track saved to {output}")

    plt.close(fig)
    return fig


def plot_sv_size_distribution(
    variants: list[dict[str, Any]],
    sv_type: str | None = None,
    output_path: str | Path = "",
    title: str = "SV Size Distribution",
    figsize: tuple[float, float] = (10, 6),
    log_scale: bool = True,
) -> Any:
    """Plot size distribution of structural variants.

    Creates a histogram showing the distribution of SV sizes, optionally
    filtered by SV type. Uses log-scale x-axis by default since SV sizes
    span several orders of magnitude.

    Args:
        variants: List of variant dictionaries with 'size' or 'start'/'end',
            and 'sv_type' keys.
        sv_type: If specified, only plot variants of this type.
        output_path: Path to save the plot. If empty, plot is not saved.
        title: Plot title.
        figsize: Figure size in inches.
        log_scale: Whether to use log scale for x-axis (default True).

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    # Collect sizes by type
    sizes_by_type: dict[str, list[float]] = {}

    for v in variants:
        vtype = v.get("sv_type", "UNKNOWN")
        if hasattr(vtype, "value"):
            vtype = vtype.value

        if sv_type is not None and vtype != sv_type:
            continue

        size = v.get("size", 0)
        if size == 0:
            size = abs(v.get("end", 0) - v.get("start", 0))
        if size <= 0:
            continue

        if vtype not in sizes_by_type:
            sizes_by_type[vtype] = []
        sizes_by_type[vtype].append(float(size))

    if not sizes_by_type:
        logger.warning("No variants to plot in size distribution")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No variants to display", ha="center", va="center", transform=ax.transAxes)
        plt.close(fig)
        return fig

    fig, ax = plt.subplots(figsize=figsize)

    if log_scale and np is not None:
        # Use log-spaced bins
        all_sizes = [s for sizes in sizes_by_type.values() for s in sizes]
        min_size = max(1, min(all_sizes))
        max_size = max(all_sizes)
        bins = np.logspace(np.log10(min_size), np.log10(max_size), 50)
    else:
        bins = 50  # type: ignore[assignment]

    for vtype in sorted(sizes_by_type.keys()):
        color = SV_COLORS.get(vtype, SV_COLORS["UNKNOWN"])
        ax.hist(
            sizes_by_type[vtype],
            bins=bins,
            alpha=0.6,
            color=color,
            label=f"{vtype} (n={len(sizes_by_type[vtype])})",
            edgecolor="white",
            linewidth=0.5,
        )

    if log_scale:
        ax.set_xscale("log")

    ax.set_xlabel("SV Size (bp)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)

    # Add size category annotations
    for size, label in [(50, "50bp"), (1000, "1kb"), (100_000, "100kb"), (1_000_000, "1Mb")]:
        ax.axvline(x=size, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()

    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output), dpi=150, bbox_inches="tight")
        logger.info(f"Size distribution plot saved to {output}")

    plt.close(fig)
    return fig


def plot_sv_type_summary(
    variants: list[dict[str, Any]],
    output_path: str | Path = "",
    title: str = "Structural Variant Type Summary",
    figsize: tuple[float, float] = (12, 5),
) -> Any:
    """Plot summary of SV types as bar chart and pie chart.

    Creates a two-panel figure showing (left) a bar chart of SV counts
    by type and (right) a pie chart of the same data.

    Args:
        variants: List of variant dictionaries with 'sv_type' keys.
        output_path: Path to save the plot.
        title: Plot title.
        figsize: Figure size.

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    # Count by type
    type_counts: dict[str, int] = {}
    for v in variants:
        vtype = v.get("sv_type", "UNKNOWN")
        if hasattr(vtype, "value"):
            vtype = vtype.value
        type_counts[vtype] = type_counts.get(vtype, 0) + 1

    if not type_counts:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No variants to display", ha="center", va="center", transform=ax.transAxes)
        plt.close(fig)
        return fig

    types = sorted(type_counts.keys())
    counts = [type_counts[t] for t in types]
    colors = [SV_COLORS.get(t, SV_COLORS["UNKNOWN"]) for t in types]

    fig, (ax_bar, ax_pie) = plt.subplots(1, 2, figsize=figsize)

    # Bar chart
    bars = ax_bar.bar(types, counts, color=colors, edgecolor="white", linewidth=1)
    ax_bar.set_xlabel("SV Type", fontsize=12)
    ax_bar.set_ylabel("Count", fontsize=12)
    ax_bar.set_title("Count by Type", fontsize=12)

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.02,
            str(count),
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    # Pie chart
    ax_pie.pie(
        counts,
        labels=[f"{t}\n({c})" for t, c in zip(types, counts)],
        colors=colors,
        autopct="%1.1f%%",
        startangle=90,
        textprops={"fontsize": 10},
    )
    ax_pie.set_title("Proportion by Type", fontsize=12)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output), dpi=150, bbox_inches="tight")
        logger.info(f"SV type summary saved to {output}")

    plt.close(fig)
    return fig


def plot_breakpoint_detail(
    variant: dict[str, Any],
    reads: list[dict[str, Any]],
    flanking: int = 500,
    output_path: str | Path = "",
    figsize: tuple[float, float] = (14, 8),
) -> Any:
    """Plot detailed view of a structural variant breakpoint.

    Creates a multi-panel figure showing:
    - Top: Read alignment pileup around the breakpoint
    - Middle: Coverage profile
    - Bottom: Split/discordant read evidence

    Args:
        variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
        reads: List of read alignment dictionaries.
        flanking: Number of base pairs to show flanking the breakpoint.
        output_path: Path to save the plot.
        figsize: Figure size.

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    chrom = variant.get("chrom", "")
    bp1 = variant.get("start", 0)
    bp2 = variant.get("end", 0)
    sv_type = variant.get("sv_type", "UNKNOWN")
    if hasattr(sv_type, "value"):
        sv_type = sv_type.value

    region_start = max(0, bp1 - flanking)
    region_end = bp2 + flanking

    fig, (ax_align, ax_cov, ax_evidence) = plt.subplots(3, 1, figsize=figsize, height_ratios=[3, 2, 2], sharex=True)

    # Filter reads in region
    region_reads: list[dict[str, Any]] = []
    for read in reads:
        r_chrom = read.get("chrom", "")
        r_pos = read.get("pos", 0)
        r_end = r_pos + read.get("read_length", 150)
        if r_chrom == chrom and r_pos < region_end and r_end > region_start:
            region_reads.append(read)

    # Panel 1: Read alignment pileup
    y_offset = 0
    for read in region_reads[:100]:  # Limit for readability
        r_pos = read.get("pos", 0)
        r_len = read.get("read_length", 150)
        r_end = r_pos + r_len
        name = read.get("name", "")

        # Color by evidence type
        evidence = variant.get("evidence")
        evidence_reads: list[str] = []
        if evidence is not None:
            if hasattr(evidence, "evidence_reads"):
                evidence_reads = evidence.evidence_reads
            elif isinstance(evidence, dict):
                evidence_reads = evidence.get("evidence_reads", [])

        if name in evidence_reads:
            color = SV_COLORS.get(sv_type, "#E74C3C")
        else:
            color = "#3498DB"

        ax_align.barh(y_offset, r_len, left=r_pos, height=0.8, color=color, alpha=0.6, edgecolor="none")
        y_offset += 1

    ax_align.set_ylabel("Reads", fontsize=11)
    ax_align.set_title(
        f"Breakpoint Detail - {chrom}:{bp1:,}-{bp2:,} ({sv_type})",
        fontsize=13,
        fontweight="bold",
    )

    # Draw breakpoint lines
    for ax in (ax_align, ax_cov, ax_evidence):
        ax.axvline(x=bp1, color="#E74C3C", linestyle="--", linewidth=1.5, alpha=0.8)
        ax.axvline(x=bp2, color="#E74C3C", linestyle="--", linewidth=1.5, alpha=0.8)

    # Panel 2: Coverage profile
    bin_size = max(1, (region_end - region_start) // 200)
    n_bins = (region_end - region_start) // bin_size + 1
    coverage_bins = [0.0] * n_bins

    for read in region_reads:
        r_pos = read.get("pos", 0)
        r_len = read.get("read_length", 150)
        for pos in range(max(r_pos, region_start), min(r_pos + r_len, region_end)):
            bin_idx = (pos - region_start) // bin_size
            if 0 <= bin_idx < n_bins:
                coverage_bins[bin_idx] += 1

    positions = [region_start + i * bin_size for i in range(n_bins)]
    ax_cov.fill_between(positions, coverage_bins, alpha=0.6, color="#2ECC71")
    ax_cov.plot(positions, coverage_bins, linewidth=0.5, color="#27AE60")
    ax_cov.set_ylabel("Depth", fontsize=11)

    # Panel 3: Evidence type track
    split_count = 0
    discordant_count = 0

    for read in region_reads:
        r_pos = read.get("pos", 0)
        name = read.get("name", "")
        cigar = read.get("cigar", "")

        if name in evidence_reads:
            if "S" in str(cigar):
                ax_evidence.barh(
                    split_count,
                    read.get("read_length", 150),
                    left=r_pos,
                    height=0.8,
                    color="#E74C3C",
                    alpha=0.7,
                )
                split_count += 1
            else:
                ax_evidence.barh(
                    split_count + discordant_count + 1,
                    read.get("read_length", 150),
                    left=r_pos,
                    height=0.8,
                    color="#9B59B6",
                    alpha=0.7,
                )
                discordant_count += 1

    ax_evidence.set_ylabel("Evidence", fontsize=11)
    ax_evidence.set_xlabel(f"Position on {chrom} (bp)", fontsize=11)

    # Legend
    evidence_patches = [
        mpatches.Patch(color="#E74C3C", label="Split reads"),
        mpatches.Patch(color="#9B59B6", label="Discordant pairs"),
        mpatches.Patch(color="#3498DB", label="Normal reads"),
    ]
    ax_align.legend(handles=evidence_patches, loc="upper right", fontsize=9)

    plt.tight_layout()

    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output), dpi=150, bbox_inches="tight")
        logger.info(f"Breakpoint detail plot saved to {output}")

    plt.close(fig)
    return fig


def plot_cnv_profile(
    segments: list[dict[str, Any]],
    chromosomes: dict[str, int],
    output_path: str | Path = "",
    title: str = "Genome-wide CNV Profile",
    figsize: tuple[float, float] = (16, 5),
    ploidy: int = 2,
) -> Any:
    """Plot genome-wide CNV profile.

    Displays log2 ratio values across the genome with colored segments
    indicating copy number states. Chromosomes are laid out linearly
    with alternating background shading.

    Args:
        segments: List of segment dictionaries with:
            - 'chrom': Chromosome
            - 'start': Start position
            - 'end': End position
            - 'mean_log2ratio': Mean log2 ratio
            - 'state': CNV state (DEL, DUP, NEUTRAL, etc.)
        chromosomes: Dictionary mapping chromosome names to sizes.
        output_path: Path to save the plot.
        title: Plot title.
        figsize: Figure size.
        ploidy: Expected ploidy for reference line.

    Returns:
        matplotlib Figure object.
    """
    _check_matplotlib()

    fig, ax = plt.subplots(figsize=figsize)

    # Calculate cumulative chromosome positions
    chrom_order = sorted(chromosomes.keys(), key=_chrom_sort_key)
    chrom_offsets: dict[str, int] = {}
    cumulative = 0

    for chrom in chrom_order:
        chrom_offsets[chrom] = cumulative
        cumulative += chromosomes[chrom]

    total_genome = cumulative

    # Draw chromosome backgrounds
    for i, chrom in enumerate(chrom_order):
        offset = chrom_offsets[chrom]
        size = chromosomes[chrom]
        color = "#F5F5F5" if i % 2 == 0 else "#FFFFFF"
        ax.axvspan(offset, offset + size, alpha=0.3, color=color)

        # Chromosome label
        ax.text(
            offset + size / 2,
            -2.8,
            chrom.replace("chr", ""),
            ha="center",
            va="top",
            fontsize=8,
            rotation=45,
        )

    # Plot segments
    for seg in segments:
        chrom = seg.get("chrom", "")
        if chrom not in chrom_offsets:
            continue

        offset = chrom_offsets[chrom]
        start = seg.get("start", 0) + offset
        end = seg.get("end", 0) + offset
        lr = seg.get("mean_log2ratio", 0.0)
        state = seg.get("state", "NEUTRAL")
        if hasattr(state, "value"):
            state = state.value

        color = CNV_COLORS.get(state, CNV_COLORS["NEUTRAL"])
        width = end - start

        # Draw segment as a horizontal bar
        ax.barh(
            lr,
            width,
            left=start,
            height=0.05,
            color=color,
            alpha=0.8,
            edgecolor="none",
        )

        # Also plot as a scatter point at the segment midpoint
        mid = (start + end) / 2
        ax.scatter(mid, lr, s=max(3, width / total_genome * 3000), color=color, alpha=0.6, edgecolors="none")

    # Reference lines
    ax.axhline(y=0, color="black", linewidth=1, linestyle="-", alpha=0.5)
    ax.axhline(
        y=math.log2(1 / ploidy + 1e-10) if ploidy > 0 else -1,
        color="#E74C3C",
        linewidth=0.5,
        linestyle="--",
        alpha=0.5,
        label="DEL",
    )
    ax.axhline(
        y=math.log2((ploidy + 1) / ploidy) if ploidy > 0 else 0.58,
        color="#3498DB",
        linewidth=0.5,
        linestyle="--",
        alpha=0.5,
        label="DUP",
    )

    ax.set_xlim(0, total_genome)
    ax.set_ylim(-3, 3)
    ax.set_xlabel("Genomic Position", fontsize=12)
    ax.set_ylabel("Log2 Ratio", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Legend
    legend_patches = [
        mpatches.Patch(color=CNV_COLORS["HOMODEL"], label="Homozygous DEL"),
        mpatches.Patch(color=CNV_COLORS["DEL"], label="Deletion"),
        mpatches.Patch(color=CNV_COLORS["NEUTRAL"], label="Neutral"),
        mpatches.Patch(color=CNV_COLORS["DUP"], label="Duplication"),
        mpatches.Patch(color=CNV_COLORS["AMP"], label="Amplification"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    plt.tight_layout()

    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output), dpi=150, bbox_inches="tight")
        logger.info(f"CNV profile saved to {output}")

    plt.close(fig)
    return fig


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Sort key for chromosome names (numeric first, then alphabetic).

    Args:
        chrom: Chromosome name.

    Returns:
        Sort key tuple.
    """
    name = chrom.replace("chr", "").replace("Chr", "")
    try:
        return (0, f"{int(name):05d}")
    except ValueError:
        return (1, name)


def _pos_to_angle(
    chrom: str,
    pos: int,
    chrom_angles: dict[str, tuple[float, float]],
    chromosomes: dict[str, int],
) -> float | None:
    """Convert a genomic position to a polar angle.

    Args:
        chrom: Chromosome.
        pos: Genomic position.
        chrom_angles: Chromosome angular ranges.
        chromosomes: Chromosome sizes.

    Returns:
        Polar angle in radians, or None if chromosome not found.
    """
    if chrom not in chrom_angles:
        return None

    start_angle, end_angle = chrom_angles[chrom]
    chrom_size = chromosomes.get(chrom, 1)

    fraction = min(1.0, max(0.0, pos / chrom_size))
    return start_angle + fraction * (end_angle - start_angle)


def _draw_chromosome_arcs(
    ax: Any,
    chrom_angles: dict[str, tuple[float, float]],
    outer_r: float,
    inner_r: float,
    chrom_order: list[str],
    show_labels: bool,
) -> None:
    """Draw chromosome arcs around the Circos plot.

    Args:
        ax: Polar axes.
        chrom_angles: Angular ranges per chromosome.
        outer_r: Outer radius.
        inner_r: Inner radius.
        chrom_order: Ordered chromosome names.
        show_labels: Whether to show labels.
    """
    # Alternating colors for chromosomes
    alt_colors = ["#3498DB", "#2ECC71"]

    for i, chrom in enumerate(chrom_order):
        if chrom not in chrom_angles:
            continue

        start_angle, end_angle = chrom_angles[chrom]
        color = alt_colors[i % 2]

        # Draw arc
        n_points = max(10, int((end_angle - start_angle) * 50))
        angles = [start_angle + j * (end_angle - start_angle) / n_points for j in range(n_points + 1)]

        for j in range(len(angles) - 1):
            theta = [angles[j], angles[j + 1], angles[j + 1], angles[j]]
            r = [inner_r, inner_r, outer_r, outer_r]
            ax.fill(theta, r, color=color, alpha=0.7)

        # Label
        if show_labels:
            mid_angle = (start_angle + end_angle) / 2
            label = chrom.replace("chr", "")
            ax.text(
                mid_angle,
                outer_r + 0.05,
                label,
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
            )


def _draw_chord(
    ax: Any,
    angle1: float,
    angle2: float,
    radius: float,
    color: str,
    alpha: float = 0.5,
) -> None:
    """Draw a chord between two angles in polar coordinates.

    Uses a quadratic Bezier curve through the center for visual clarity.

    Args:
        ax: Polar axes.
        angle1: First angle in radians.
        angle2: Second angle in radians.
        radius: Radius at which the chord endpoints are placed.
        color: Color of the chord.
        alpha: Transparency.
    """
    # Number of interpolation points
    n_points = 50

    # Create a curved path through the center
    t_values = [i / n_points for i in range(n_points + 1)]

    angles_path = []
    radii_path = []

    for t in t_values:
        # Interpolate angle
        angle = angle1 + t * (angle2 - angle1)
        # Radius follows a parabolic curve (dips toward center)
        r = radius * (1 - 4 * t * (1 - t) * 0.7)  # Dips to 30% of radius at midpoint
        angles_path.append(angle)
        radii_path.append(r)

    ax.plot(angles_path, radii_path, color=color, alpha=alpha, linewidth=1)
