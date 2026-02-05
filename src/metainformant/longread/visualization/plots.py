"""Visualization tools for long-read sequencing data.

Generates publication-quality plots for read length distributions,
quality vs length scatter plots, sequence dotplots, alignment views,
methylation tracks, and phasing block diagrams.

All functions save plots to disk and return the output path. Uses
matplotlib and seaborn for rendering.

Optional dependencies:
    - matplotlib: For plot generation
    - seaborn: For enhanced styling
    - numpy: For data processing
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import matplotlib  # type: ignore[import-untyped]

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt  # type: ignore[import-untyped]
    import matplotlib.patches as mpatches  # type: ignore[import-untyped]
    from matplotlib.collections import LineCollection  # type: ignore[import-untyped]

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import seaborn as sns  # type: ignore[import-untyped]

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    import numpy as np  # type: ignore[import-untyped]

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def _ensure_matplotlib() -> None:
    """Ensure matplotlib is available."""
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib is required for visualization. " "Install it with: uv pip install matplotlib")


def _apply_style() -> None:
    """Apply consistent plot styling."""
    if HAS_SEABORN:
        sns.set_theme(style="whitegrid", palette="deep")
    else:
        plt.style.use("seaborn-v0_8-whitegrid") if "seaborn-v0_8-whitegrid" in plt.style.available else None


def plot_read_length_histogram(
    reads: Sequence[int | dict[str, Any]],
    output_path: str | Path,
    log_scale: bool = True,
    bins: int = 100,
    title: str = "Read Length Distribution",
    color: str = "#2196F3",
    figsize: tuple[float, float] = (10, 6),
) -> Path:
    """Plot a histogram of read length distribution.

    Generates a publication-quality histogram showing the distribution of
    read lengths. Optionally uses log scale for the x-axis to better
    visualize the typical long-tail distribution of long reads.

    Adds N50 and median lines as vertical annotations.

    Args:
        reads: Sequence of read lengths (int) or read dicts with 'length' or 'sequence'.
        output_path: Path to save the plot.
        log_scale: Whether to use log scale for the x-axis.
        bins: Number of histogram bins.
        title: Plot title.
        color: Bar color.
        figsize: Figure size as (width, height).

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    # Extract lengths
    lengths = _extract_lengths(reads)
    if not lengths:
        logger.warning("No read lengths to plot")
        return Path(output_path)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    if log_scale and all(l > 0 for l in lengths):
        # Log-spaced bins
        log_min = math.log10(min(lengths))
        log_max = math.log10(max(lengths))
        if HAS_NUMPY:
            bin_edges = np.logspace(log_min, log_max, bins + 1)
        else:
            step = (log_max - log_min) / bins
            bin_edges = [10 ** (log_min + i * step) for i in range(bins + 1)]
        ax.hist(lengths, bins=bin_edges, color=color, alpha=0.8, edgecolor="white", linewidth=0.5)
        ax.set_xscale("log")
    else:
        ax.hist(lengths, bins=bins, color=color, alpha=0.8, edgecolor="white", linewidth=0.5)

    # Add N50 and median lines
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumulative = 0
    n50 = sorted_lengths[-1]
    for l in sorted_lengths:
        cumulative += l
        if cumulative >= total * 0.5:
            n50 = l
            break

    median = sorted(lengths)[len(lengths) // 2]

    ax.axvline(n50, color="#FF5722", linestyle="--", linewidth=2, label=f"N50 = {n50:,} bp")
    ax.axvline(median, color="#4CAF50", linestyle="--", linewidth=2, label=f"Median = {median:,} bp")

    ax.set_xlabel("Read Length (bp)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=10)

    # Add summary text
    summary_text = (
        f"Reads: {len(lengths):,}\n"
        f"Total: {total / 1e6:.1f} Mb\n"
        f"Mean: {total / len(lengths):,.0f} bp\n"
        f"Max: {max(lengths):,} bp"
    )
    ax.text(
        0.98,
        0.95,
        summary_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved read length histogram to %s", output_path)
    return output_path


def plot_quality_vs_length(
    reads: Sequence[dict[str, Any]],
    output_path: str | Path,
    title: str = "Quality vs Read Length",
    figsize: tuple[float, float] = (10, 6),
) -> Path:
    """Plot a scatter plot of read quality vs read length.

    Creates a 2D scatter plot (or hexbin for large datasets) showing the
    relationship between read length and mean quality score.

    Args:
        reads: Sequence of read dicts with 'sequence' (or 'length') and
            'quality_string' (or 'quality') keys.
        output_path: Path to save the plot.
        title: Plot title.
        figsize: Figure size.

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    lengths: list[int] = []
    qualities: list[float] = []

    for read in reads:
        length = 0
        qual = 0.0

        if isinstance(read, dict):
            if "length" in read:
                length = int(read["length"])
            elif "sequence" in read and read["sequence"]:
                length = len(read["sequence"])

            qual_str = read.get("quality_string") or read.get("quality")
            if isinstance(qual_str, str):
                scores = [ord(c) - 33 for c in qual_str]
                qual = sum(scores) / len(scores) if scores else 0
            elif isinstance(qual_str, (int, float)):
                qual = float(qual_str)

        if length > 0 and qual > 0:
            lengths.append(length)
            qualities.append(qual)

    if not lengths:
        logger.warning("No reads with both length and quality data")
        return Path(output_path)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    if len(lengths) > 10000 and HAS_NUMPY:
        # Use hexbin for large datasets
        hb = ax.hexbin(lengths, qualities, gridsize=50, cmap="YlOrRd", mincnt=1)
        fig.colorbar(hb, ax=ax, label="Count")
    else:
        ax.scatter(lengths, qualities, alpha=0.3, s=8, c="#2196F3", edgecolors="none")

    ax.set_xlabel("Read Length (bp)", fontsize=12)
    ax.set_ylabel("Mean Quality Score (Phred)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Add quality threshold lines
    for q_threshold in [7, 10, 20]:
        ax.axhline(q_threshold, color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
        ax.text(max(lengths) * 0.98, q_threshold + 0.3, f"Q{q_threshold}", fontsize=8, color="gray", ha="right")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved quality vs length plot to %s", output_path)
    return output_path


def plot_dotplot(
    seq1: str,
    seq2: str,
    output_path: str | Path,
    word_size: int = 11,
    title: str = "Sequence Dotplot",
    figsize: tuple[float, float] = (10, 10),
    seq1_name: str = "Sequence 1",
    seq2_name: str = "Sequence 2",
) -> Path:
    """Plot a sequence dotplot comparing two sequences.

    Generates a dotplot by finding shared k-mers (words) between two sequences
    and plotting their positions. Forward matches appear on the main diagonal
    and reverse complement matches appear on the anti-diagonal.

    Args:
        seq1: First sequence.
        seq2: Second sequence.
        output_path: Path to save the plot.
        word_size: K-mer size for matches (default 11).
        title: Plot title.
        figsize: Figure size.
        seq1_name: Label for sequence 1 (x-axis).
        seq2_name: Label for sequence 2 (y-axis).

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seq1_upper = seq1.upper()
    seq2_upper = seq2.upper()

    # Build k-mer index for seq1
    seq1_kmers: dict[str, list[int]] = {}
    for i in range(len(seq1_upper) - word_size + 1):
        kmer = seq1_upper[i : i + word_size]
        if "N" not in kmer:
            if kmer not in seq1_kmers:
                seq1_kmers[kmer] = []
            seq1_kmers[kmer].append(i)

    # Find matches in seq2
    fwd_x: list[int] = []
    fwd_y: list[int] = []
    rev_x: list[int] = []
    rev_y: list[int] = []

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    for j in range(len(seq2_upper) - word_size + 1):
        kmer = seq2_upper[j : j + word_size]
        if "N" in kmer:
            continue

        # Forward matches
        if kmer in seq1_kmers:
            for i in seq1_kmers[kmer]:
                fwd_x.append(i)
                fwd_y.append(j)

        # Reverse complement matches
        rc_kmer = "".join(complement.get(c, "N") for c in reversed(kmer))
        if rc_kmer in seq1_kmers:
            for i in seq1_kmers[rc_kmer]:
                rev_x.append(i)
                rev_y.append(j)

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    if fwd_x:
        ax.scatter(fwd_x, fwd_y, s=0.5, c="#2196F3", alpha=0.6, label="Forward")
    if rev_x:
        ax.scatter(rev_x, rev_y, s=0.5, c="#FF5722", alpha=0.6, label="Reverse complement")

    ax.set_xlabel(f"{seq1_name} (bp)", fontsize=12)
    ax.set_ylabel(f"{seq2_name} (bp)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlim(0, len(seq1))
    ax.set_ylim(0, len(seq2))

    if fwd_x or rev_x:
        ax.legend(fontsize=10, markerscale=10)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved dotplot to %s (%d fwd, %d rev matches)", output_path, len(fwd_x), len(rev_x))
    return output_path


def plot_alignment_view(
    alignments: Sequence[dict[str, Any]],
    region: str | dict[str, Any],
    output_path: str | Path,
    title: str | None = None,
    figsize: tuple[float, float] = (14, 8),
    max_reads: int = 100,
) -> Path:
    """Plot an IGV-style alignment view for a genomic region.

    Displays reads as horizontal bars stacked in rows, with colors
    indicating alignment quality and arrows showing read direction.

    Args:
        alignments: Sequence of alignment dicts with keys:
            reference_start, reference_end, is_reverse, mapping_quality,
            read_name, is_supplementary (optional).
        region: Region specification as "chr:start-end" string or dict with
            chromosome, start, end keys.
        output_path: Path to save the plot.
        title: Optional plot title. Defaults to the region string.
        figsize: Figure size.
        max_reads: Maximum number of reads to display.

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Parse region
    if isinstance(region, str):
        parts = region.replace("-", ":").split(":")
        chrom = parts[0]
        reg_start = int(parts[1]) if len(parts) > 1 else 0
        reg_end = int(parts[2]) if len(parts) > 2 else reg_start + 10000
    else:
        chrom = region.get("chromosome", "")
        reg_start = int(region.get("start", 0))
        reg_end = int(region.get("end", reg_start + 10000))

    if title is None:
        title = f"Alignment View: {chrom}:{reg_start:,}-{reg_end:,}"

    # Filter alignments to region
    region_alns: list[dict[str, Any]] = []
    for aln in alignments:
        aln_start = aln.get("reference_start", 0)
        aln_end = aln.get("reference_end", aln_start)
        aln_chrom = aln.get("reference_name", "")

        if aln_chrom == chrom and aln_start < reg_end and aln_end > reg_start:
            region_alns.append(aln)

    # Sort and limit
    region_alns.sort(key=lambda a: a.get("reference_start", 0))
    region_alns = region_alns[:max_reads]

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    # Pack reads into rows (greedy interval scheduling)
    rows: list[int] = []  # End position of the last read in each row
    row_assignments: list[int] = []

    for aln in region_alns:
        aln_start = max(aln.get("reference_start", 0), reg_start)
        aln_end = min(aln.get("reference_end", aln_start), reg_end)

        # Find the first row where this read fits
        placed = False
        for row_idx, row_end in enumerate(rows):
            if aln_start > row_end + 50:  # Small gap between reads
                rows[row_idx] = aln_end
                row_assignments.append(row_idx)
                placed = True
                break

        if not placed:
            rows.append(aln_end)
            row_assignments.append(len(rows) - 1)

    # Draw reads
    for i, aln in enumerate(region_alns):
        aln_start = max(aln.get("reference_start", 0), reg_start)
        aln_end = min(aln.get("reference_end", aln_start), reg_end)
        row = row_assignments[i]
        is_reverse = aln.get("is_reverse", False)
        mapq = aln.get("mapping_quality", 60)
        is_supp = aln.get("is_supplementary", False)

        # Color based on mapping quality
        if is_supp:
            color = "#FF9800"  # Orange for supplementary
        elif mapq >= 60:
            color = "#BBDEFB" if not is_reverse else "#FFCDD2"
        elif mapq >= 20:
            color = "#90CAF9" if not is_reverse else "#EF9A9A"
        else:
            color = "#64B5F6" if not is_reverse else "#E57373"

        # Draw read as rectangle
        height = 0.8
        y = row * 1.0
        width = aln_end - aln_start

        rect = mpatches.FancyBboxPatch(
            (aln_start, y),
            width,
            height,
            boxstyle="round,pad=0",
            facecolor=color,
            edgecolor="gray",
            linewidth=0.3,
        )
        ax.add_patch(rect)

    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(-0.5, len(rows) + 0.5)
    ax.set_xlabel(f"Position on {chrom} (bp)", fontsize=12)
    ax.set_ylabel("Reads", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Add legend
    legend_patches = [
        mpatches.Patch(color="#BBDEFB", label="Forward"),
        mpatches.Patch(color="#FFCDD2", label="Reverse"),
        mpatches.Patch(color="#FF9800", label="Supplementary"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    # Remove y-axis ticks (row numbers not meaningful)
    ax.set_yticks([])

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved alignment view to %s (%d reads)", output_path, len(region_alns))
    return output_path


def plot_methylation_track(
    methylation_data: Sequence[dict[str, Any]],
    region: str | dict[str, Any],
    output_path: str | Path,
    title: str | None = None,
    figsize: tuple[float, float] = (14, 4),
    threshold: float = 0.5,
) -> Path:
    """Plot a methylation track showing methylation levels along a genomic region.

    Displays per-CpG methylation levels as colored bars, with unmethylated
    sites in blue and methylated sites in red.

    Args:
        methylation_data: Sequence of methylation call dicts with keys:
            position (int), probability (float), chromosome (str, optional).
        region: Region specification as "chr:start-end" or dict.
        output_path: Path to save the plot.
        title: Optional plot title.
        figsize: Figure size.
        threshold: Methylation probability threshold for coloring.

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Parse region
    if isinstance(region, str):
        parts = region.replace("-", ":").split(":")
        chrom = parts[0]
        reg_start = int(parts[1]) if len(parts) > 1 else 0
        reg_end = int(parts[2]) if len(parts) > 2 else reg_start + 10000
    else:
        chrom = region.get("chromosome", "")
        reg_start = int(region.get("start", 0))
        reg_end = int(region.get("end", reg_start + 10000))

    if title is None:
        title = f"Methylation: {chrom}:{reg_start:,}-{reg_end:,}"

    # Filter data to region
    positions: list[int] = []
    levels: list[float] = []

    for call in methylation_data:
        pos = call.get("position", 0)
        prob = float(call.get("probability", call.get("methylation_level", 0.0)))
        call_chrom = call.get("chromosome", chrom)

        if call_chrom == chrom and reg_start <= pos < reg_end:
            positions.append(pos)
            levels.append(prob)

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    if positions:
        # Color by methylation level
        colors = ["#2196F3" if l < threshold else "#F44336" for l in levels]
        ax.bar(
            positions,
            levels,
            width=max(1, (reg_end - reg_start) // len(positions) // 2),
            color=colors,
            alpha=0.8,
            edgecolor="none",
        )

        ax.axhline(threshold, color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
    else:
        ax.text(
            0.5,
            0.5,
            "No methylation data in region",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="gray",
        )

    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(f"Position on {chrom} (bp)", fontsize=12)
    ax.set_ylabel("Methylation Level", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Legend
    legend_patches = [
        mpatches.Patch(color="#F44336", label="Methylated"),
        mpatches.Patch(color="#2196F3", label="Unmethylated"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved methylation track to %s (%d CpGs)", output_path, len(positions))
    return output_path


def plot_phasing_blocks(
    phase_blocks: Sequence[dict[str, Any]],
    output_path: str | Path,
    title: str = "Haplotype Phase Blocks",
    figsize: tuple[float, float] = (14, 6),
) -> Path:
    """Plot haplotype phase blocks as colored horizontal bars.

    Displays phase blocks along the genome, with alternating colors for
    haplotype 1 and haplotype 2, and block boundaries clearly marked.

    Args:
        phase_blocks: Sequence of phase block dicts with keys:
            chromosome (str), start (int), end (int),
            num_variants (int, optional), quality (float, optional).
        output_path: Path to save the plot.
        title: Plot title.
        figsize: Figure size.

    Returns:
        Path to the saved plot file.
    """
    _ensure_matplotlib()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not phase_blocks:
        logger.warning("No phase blocks to plot")
        # Create empty plot
        _apply_style()
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5, 0.5, "No phase blocks", transform=ax.transAxes, ha="center", va="center", fontsize=14, color="gray"
        )
        ax.set_title(title, fontsize=14, fontweight="bold")
        fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
        return output_path

    # Normalize block data
    blocks: list[dict[str, Any]] = []
    for b in phase_blocks:
        if isinstance(b, dict):
            blocks.append(b)
        elif hasattr(b, "chromosome"):
            blocks.append(
                {
                    "chromosome": b.chromosome,
                    "start": b.start,
                    "end": b.end,
                    "num_variants": getattr(b, "num_variants", 0),
                    "quality": getattr(b, "quality", 0.0),
                }
            )

    # Group by chromosome
    chroms: dict[str, list[dict[str, Any]]] = {}
    for b in blocks:
        chrom = b.get("chromosome", "unknown")
        if chrom not in chroms:
            chroms[chrom] = []
        chroms[chrom].append(b)

    # Sort blocks within each chromosome
    for chrom in chroms:
        chroms[chrom].sort(key=lambda b: b.get("start", 0))

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    colors = ["#2196F3", "#FF9800"]
    y_pos = 0
    y_labels: list[str] = []
    y_positions: list[float] = []

    for chrom, chrom_blocks in sorted(chroms.items()):
        y_labels.append(chrom)
        y_positions.append(y_pos)

        for i, block in enumerate(chrom_blocks):
            start = block.get("start", 0)
            end = block.get("end", start)
            width = end - start
            num_vars = block.get("num_variants", 0)
            quality = block.get("quality", 0.0)

            color = colors[i % 2]
            alpha = min(1.0, 0.5 + quality * 0.5) if quality > 0 else 0.8

            rect = mpatches.FancyBboxPatch(
                (start, y_pos - 0.35),
                width,
                0.7,
                boxstyle="round,pad=0",
                facecolor=color,
                alpha=alpha,
                edgecolor="black",
                linewidth=0.5,
            )
            ax.add_patch(rect)

            # Annotate with variant count
            if num_vars > 0 and width > 0:
                ax.text(
                    start + width / 2,
                    y_pos,
                    f"{num_vars}",
                    fontsize=7,
                    ha="center",
                    va="center",
                    color="white",
                    fontweight="bold",
                )

        y_pos += 1

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=10)
    ax.set_xlabel("Genomic Position (bp)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Auto-scale x-axis
    all_starts = [b.get("start", 0) for b in blocks]
    all_ends = [b.get("end", 0) for b in blocks]
    if all_starts and all_ends:
        ax.set_xlim(min(all_starts) - 1000, max(all_ends) + 1000)

    ax.set_ylim(-0.5, y_pos - 0.5)

    # Legend
    legend_patches = [
        mpatches.Patch(color=colors[0], label="Block (odd)"),
        mpatches.Patch(color=colors[1], label="Block (even)"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved phase block plot to %s (%d blocks)", output_path, len(blocks))
    return output_path


def _extract_lengths(reads: Sequence[int | dict[str, Any]]) -> list[int]:
    """Extract read lengths from various representations."""
    lengths: list[int] = []
    for r in reads:
        if isinstance(r, int):
            if r > 0:
                lengths.append(r)
        elif isinstance(r, dict):
            if "length" in r:
                lengths.append(int(r["length"]))
            elif "sequence" in r and r["sequence"]:
                lengths.append(len(r["sequence"]))
        elif hasattr(r, "sequence") and r.sequence:
            lengths.append(len(r.sequence))
        elif hasattr(r, "query_length"):
            lengths.append(int(r.query_length))
    return lengths
