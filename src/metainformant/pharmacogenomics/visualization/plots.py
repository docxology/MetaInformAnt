"""Pharmacogenomics visualization functions.

Provides publication-quality plots for metabolizer status distributions,
allele frequencies, activity score distributions, drug-gene response heatmaps,
cross-population allele frequency comparisons, and ACMG criteria summaries.

All plot functions accept an output_path parameter and save figures to disk.
They use matplotlib/seaborn with graceful degradation when these optional
dependencies are not available.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend for server/CLI use
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None  # type: ignore
    mpatches = None  # type: ignore

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None  # type: ignore


def _ensure_matplotlib() -> None:
    """Raise ImportError if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for pharmacogenomics visualization. " "Install with: uv pip install matplotlib"
        )


def _ensure_output_dir(output_path: str | Path) -> Path:
    """Ensure the output directory exists and return the Path."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


# ── Phenotype color palettes ──────────────────────────────────────────────────

_PHENOTYPE_COLORS = {
    "Poor Metabolizer": "#D32F2F",
    "PM": "#D32F2F",
    "Intermediate Metabolizer": "#F57C00",
    "IM": "#F57C00",
    "Normal Metabolizer": "#388E3C",
    "NM": "#388E3C",
    "Rapid Metabolizer": "#1976D2",
    "RM": "#1976D2",
    "Ultrarapid Metabolizer": "#7B1FA2",
    "UM": "#7B1FA2",
    "Indeterminate": "#9E9E9E",
    "IND": "#9E9E9E",
}

_SEVERITY_COLORS = {
    "Major": "#D32F2F",
    "Moderate": "#F57C00",
    "Minor": "#388E3C",
    "None": "#9E9E9E",
}


def plot_metabolizer_status(
    phenotypes: dict[str, str | int | float],
    output_path: str | Path,
    title: str = "Metabolizer Status Distribution",
    figsize: tuple[float, float] = (10, 6),
) -> Path:
    """Plot metabolizer status distribution as a bar chart.

    Creates a bar chart showing the distribution of metabolizer phenotypes,
    color-coded by category (PM=red, IM=orange, NM=green, RM=blue, UM=purple).

    Args:
        phenotypes: Dictionary mapping phenotype name/abbreviation to count or frequency.
            Example: {"PM": 5, "IM": 20, "NM": 65, "RM": 8, "UM": 2}
        output_path: Path to save the figure
        title: Plot title
        figsize: Figure size in inches (width, height)

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    categories = list(phenotypes.keys())
    values = [float(v) for v in phenotypes.values()]
    colors = [_PHENOTYPE_COLORS.get(c, "#607D8B") for c in categories]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(categories, values, color=colors, edgecolor="white", linewidth=0.5)

    # Add value labels on bars
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + max(values) * 0.01,
            f"{val:.1f}" if isinstance(val, float) and val != int(val) else str(int(val)),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_xlabel("Metabolizer Phenotype", fontsize=12)
    ax.set_ylabel("Count / Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved metabolizer status plot to %s", path)
    return path


def plot_allele_frequencies(
    allele_freqs: dict[str, float],
    gene: str,
    output_path: str | Path,
    plot_type: str = "bar",
    title: str | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> Path:
    """Plot allele frequency distribution for a pharmacogene.

    Creates either a bar chart or pie chart showing the frequency of each
    star allele in a population.

    Args:
        allele_freqs: Dictionary mapping allele name to frequency.
            Example: {"*1": 0.45, "*2": 0.25, "*4": 0.15, "*10": 0.10, "*17": 0.05}
        gene: Gene symbol for the title
        output_path: Path to save the figure
        plot_type: "bar" for bar chart, "pie" for pie chart
        title: Optional custom title
        figsize: Figure size in inches

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    alleles = list(allele_freqs.keys())
    freqs = list(allele_freqs.values())
    plot_title = title or f"{gene} Allele Frequencies"

    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == "pie":
        # Filter out very small slices for readability
        threshold = 0.02
        labels_display = []
        sizes_display = []
        other_sum = 0.0
        for a, f in zip(alleles, freqs):
            if f >= threshold:
                labels_display.append(f"{a}\n({f:.1%})")
                sizes_display.append(f)
            else:
                other_sum += f
        if other_sum > 0:
            labels_display.append(f"Other\n({other_sum:.1%})")
            sizes_display.append(other_sum)

        colors_pie = plt.cm.Set3(range(len(labels_display)))
        wedges, texts, autotexts = ax.pie(
            sizes_display,
            labels=labels_display,
            colors=colors_pie,
            autopct="%1.1f%%",
            startangle=90,
            pctdistance=0.75,
        )
        for text in autotexts:
            text.set_fontsize(8)
        ax.set_title(plot_title, fontsize=14, fontweight="bold")
    else:
        # Bar chart
        colors_bar = plt.cm.viridis([i / max(len(alleles) - 1, 1) for i in range(len(alleles))])
        bars = ax.bar(alleles, freqs, color=colors_bar, edgecolor="white", linewidth=0.5)

        for bar, freq in zip(bars, freqs):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + max(freqs) * 0.01,
                f"{freq:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        ax.set_xlabel("Star Allele", fontsize=12)
        ax.set_ylabel("Frequency", fontsize=12)
        ax.set_title(plot_title, fontsize=14, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if len(alleles) > 8:
            plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved allele frequency plot to %s", path)
    return path


def plot_activity_score_distribution(
    scores: list[float],
    gene: str,
    output_path: str | Path,
    title: str | None = None,
    figsize: tuple[float, float] = (10, 6),
) -> Path:
    """Plot activity score distribution as a histogram.

    Shows the distribution of activity scores across a cohort, with vertical
    lines indicating phenotype threshold boundaries.

    Args:
        scores: List of activity score values
        gene: Gene symbol
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size in inches

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    plot_title = title or f"{gene} Activity Score Distribution"

    fig, ax = plt.subplots(figsize=figsize)

    # Histogram
    if HAS_NUMPY:
        bins = np.arange(0, max(scores) + 0.5, 0.25) if scores else [0, 0.5, 1.0]
    else:
        max_score = max(scores) if scores else 1.0
        bins = [i * 0.25 for i in range(int(max_score / 0.25) + 3)]

    ax.hist(scores, bins=bins, color="#1976D2", edgecolor="white", alpha=0.8, linewidth=0.5)

    # Add phenotype threshold lines (gene-specific)
    from ..alleles.phenotype import _PHENOTYPE_THRESHOLDS

    gene_upper = gene.upper()
    thresholds = _PHENOTYPE_THRESHOLDS.get(gene_upper, [])

    threshold_values_added: set[float] = set()
    for low, high, phenotype in thresholds:
        if low > 0 and low not in threshold_values_added:
            ax.axvline(
                x=low, color=_PHENOTYPE_COLORS.get(phenotype.value, "#888"), linestyle="--", alpha=0.7, linewidth=1.5
            )
            ax.text(
                low,
                ax.get_ylim()[1] * 0.95,
                phenotype.abbreviation,
                ha="center",
                fontsize=9,
                color=_PHENOTYPE_COLORS.get(phenotype.value, "#888"),
                fontweight="bold",
            )
            threshold_values_added.add(low)

    ax.set_xlabel("Activity Score", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(plot_title, fontsize=14, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved activity score distribution plot to %s", path)
    return path


def plot_drug_response_heatmap(
    drugs: list[str],
    genes: list[str],
    phenotypes: dict[tuple[str, str], str],
    output_path: str | Path,
    title: str = "Drug-Gene Pharmacogenomic Response Matrix",
    figsize: tuple[float, float] | None = None,
) -> Path:
    """Plot drug-gene response heatmap.

    Creates a matrix visualization showing the interaction severity or
    metabolizer impact for each drug-gene combination.

    Args:
        drugs: List of drug names (rows)
        genes: List of gene symbols (columns)
        phenotypes: Dict mapping (drug, gene) tuple to phenotype/severity string.
            Example: {("codeine", "CYP2D6"): "Major", ("warfarin", "CYP2C9"): "Moderate"}
        output_path: Path to save the figure
        title: Plot title
        figsize: Optional figure size (auto-calculated if None)

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    if figsize is None:
        figsize = (max(8, len(genes) * 1.5), max(6, len(drugs) * 0.8))

    # Build numeric matrix for coloring
    severity_to_num = {"Major": 3, "Moderate": 2, "Minor": 1, "None": 0, "": 0}
    phenotype_to_num = {
        "Poor Metabolizer": 0,
        "PM": 0,
        "Intermediate Metabolizer": 1,
        "IM": 1,
        "Normal Metabolizer": 2,
        "NM": 2,
        "Rapid Metabolizer": 3,
        "RM": 3,
        "Ultrarapid Metabolizer": 4,
        "UM": 4,
    }

    matrix: list[list[float]] = []
    annotations: list[list[str]] = []

    for drug in drugs:
        row_vals: list[float] = []
        row_annots: list[str] = []
        for gene in genes:
            val = phenotypes.get((drug, gene), phenotypes.get((drug.lower(), gene.upper()), ""))
            # Try severity mapping first, then phenotype
            num = severity_to_num.get(val, phenotype_to_num.get(val, -1))
            row_vals.append(float(num) if num >= 0 else float("nan"))
            row_annots.append(val if val else "-")
        matrix.append(row_vals)
        annotations.append(row_annots)

    fig, ax = plt.subplots(figsize=figsize)

    if HAS_SEABORN and HAS_NUMPY:
        import numpy as np_local

        matrix_arr = np_local.array(matrix)
        sns.heatmap(
            matrix_arr,
            xticklabels=genes,
            yticklabels=drugs,
            annot=[[str(a) for a in row] for row in annotations],
            fmt="",
            cmap="RdYlGn_r",
            ax=ax,
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "Severity / Impact"},
        )
    else:
        # Fallback: simple colored grid
        for i, drug in enumerate(drugs):
            for j, gene in enumerate(genes):
                val = annotations[i][j]
                color = _SEVERITY_COLORS.get(val, _PHENOTYPE_COLORS.get(val, "#E0E0E0"))
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color=color, alpha=0.7))
                ax.text(j + 0.5, i + 0.5, val, ha="center", va="center", fontsize=8)

        ax.set_xlim(0, len(genes))
        ax.set_ylim(0, len(drugs))
        ax.set_xticks([i + 0.5 for i in range(len(genes))])
        ax.set_xticklabels(genes)
        ax.set_yticks([i + 0.5 for i in range(len(drugs))])
        ax.set_yticklabels(drugs)
        ax.invert_yaxis()

    ax.set_title(title, fontsize=14, fontweight="bold")
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved drug response heatmap to %s", path)
    return path


def plot_population_comparison(
    allele_freqs_by_pop: dict[str, dict[str, float]],
    gene: str,
    output_path: str | Path,
    title: str | None = None,
    figsize: tuple[float, float] = (12, 7),
) -> Path:
    """Plot cross-population allele frequency comparison.

    Creates a grouped bar chart comparing allele frequencies across
    different populations (e.g., African, European, East Asian).

    Args:
        allele_freqs_by_pop: Nested dict: population -> allele -> frequency.
            Example: {"European": {"*1": 0.6, "*2": 0.3}, "East Asian": {"*1": 0.4, "*10": 0.4}}
        gene: Gene symbol
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size in inches

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    plot_title = title or f"{gene} Allele Frequencies by Population"

    # Collect all unique alleles across populations
    all_alleles: set[str] = set()
    for pop_freqs in allele_freqs_by_pop.values():
        all_alleles.update(pop_freqs.keys())
    alleles_sorted = sorted(all_alleles)
    populations = list(allele_freqs_by_pop.keys())

    fig, ax = plt.subplots(figsize=figsize)

    n_alleles = len(alleles_sorted)
    n_pops = len(populations)
    bar_width = 0.8 / n_pops

    if HAS_NUMPY:
        x_positions = np.arange(n_alleles)
    else:
        x_positions = list(range(n_alleles))

    colors = plt.cm.Set2(range(n_pops))

    for pop_idx, population in enumerate(populations):
        pop_freqs = allele_freqs_by_pop[population]
        values = [pop_freqs.get(a, 0.0) for a in alleles_sorted]

        if HAS_NUMPY:
            positions = x_positions + pop_idx * bar_width
        else:
            positions = [x + pop_idx * bar_width for x in x_positions]

        ax.bar(
            positions,
            values,
            width=bar_width,
            label=population,
            color=colors[pop_idx],
            edgecolor="white",
            linewidth=0.5,
        )

    if HAS_NUMPY:
        ax.set_xticks(x_positions + bar_width * (n_pops - 1) / 2)
    else:
        ax.set_xticks([x + bar_width * (n_pops - 1) / 2 for x in x_positions])

    ax.set_xticklabels(alleles_sorted, rotation=45, ha="right")
    ax.set_xlabel("Star Allele", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(plot_title, fontsize=14, fontweight="bold")
    ax.legend(title="Population", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved population comparison plot to %s", path)
    return path


def plot_acmg_criteria(
    criteria_results: dict[str, bool],
    output_path: str | Path,
    title: str = "ACMG Variant Classification Criteria",
    figsize: tuple[float, float] = (14, 6),
) -> Path:
    """Plot ACMG criteria evaluation summary.

    Creates a horizontal bar chart showing which ACMG criteria are met (green)
    or not met (gray), organized by evidence strength category.

    Args:
        criteria_results: Dictionary mapping criterion name to met/not-met.
            Example: {"PVS1": True, "PS1": False, "PM2": True, "PP3": True, ...}
        output_path: Path to save the figure
        title: Plot title
        figsize: Figure size in inches

    Returns:
        Path to saved figure
    """
    _ensure_matplotlib()
    path = _ensure_output_dir(output_path)

    # Organize criteria by category
    categories = {
        "Very Strong Pathogenic": [c for c in criteria_results if c.startswith("PVS")],
        "Strong Pathogenic": [c for c in criteria_results if c.startswith("PS")],
        "Moderate Pathogenic": [c for c in criteria_results if c.startswith("PM")],
        "Supporting Pathogenic": [c for c in criteria_results if c.startswith("PP")],
        "Stand-alone Benign": [c for c in criteria_results if c.startswith("BA")],
        "Strong Benign": [c for c in criteria_results if c.startswith("BS")],
        "Supporting Benign": [c for c in criteria_results if c.startswith("BP")],
    }

    # Flatten into ordered list
    all_criteria: list[str] = []
    category_boundaries: list[tuple[int, str]] = []

    for cat_name, cat_criteria in categories.items():
        if cat_criteria:
            category_boundaries.append((len(all_criteria), cat_name))
            all_criteria.extend(sorted(cat_criteria))

    if not all_criteria:
        logger.warning("No criteria to plot")
        # Create empty plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No ACMG criteria evaluated", ha="center", va="center", fontsize=14)
        ax.set_title(title, fontsize=14, fontweight="bold")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return path

    n_criteria = len(all_criteria)
    met_values = [1 if criteria_results.get(c, False) else 0 for c in all_criteria]
    colors = ["#4CAF50" if v else "#E0E0E0" for v in met_values]

    fig, ax = plt.subplots(figsize=figsize)

    y_positions = list(range(n_criteria))
    bars = ax.barh(y_positions, [1] * n_criteria, color=colors, edgecolor="white", linewidth=0.5, height=0.7)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(all_criteria, fontsize=9)
    ax.set_xlim(0, 1.5)
    ax.set_xticks([])

    # Add met/not met labels
    for i, (criterion, is_met) in enumerate(zip(all_criteria, met_values)):
        label = "MET" if is_met else "NOT MET"
        color = "white" if is_met else "#888"
        ax.text(0.5, i, label, ha="center", va="center", fontsize=8, fontweight="bold", color=color)

    # Add category labels
    category_colors = {
        "Very Strong Pathogenic": "#B71C1C",
        "Strong Pathogenic": "#D32F2F",
        "Moderate Pathogenic": "#F57C00",
        "Supporting Pathogenic": "#FFA726",
        "Stand-alone Benign": "#1B5E20",
        "Strong Benign": "#388E3C",
        "Supporting Benign": "#66BB6A",
    }

    for idx, cat_name in category_boundaries:
        ax.text(
            1.1,
            idx,
            cat_name,
            ha="left",
            va="center",
            fontsize=8,
            color=category_colors.get(cat_name, "#333"),
            fontweight="bold",
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # Legend
    met_patch = mpatches.Patch(color="#4CAF50", label="Criterion Met")
    not_met_patch = mpatches.Patch(color="#E0E0E0", label="Criterion Not Met")
    ax.legend(handles=[met_patch, not_met_patch], loc="lower right")

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved ACMG criteria plot to %s", path)
    return path
