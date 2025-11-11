"""Amalgkit progress visualization utilities.

This module provides visualization functions for RNA-seq workflow progress tracking,
generating bar charts showing sample counts in each category across all species.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Use non-interactive backend by default
matplotlib.use("Agg", force=True)

try:
    from ..core.logging import get_logger
except ImportError:
    from metainformant.core.logging import get_logger

from .basic import bar_plot
from .export import save_figure
from .layout import adjust_spacing, create_multi_panel
from .style import apply_publication_style

logger = get_logger(__name__)

# Category color palette
CATEGORY_COLORS = {
    "need_download": "#CCCCCC",  # Light gray
    "ongoing_download": "#4A90E2",  # Blue
    "failed_download": "#E24A4A",  # Red
    "needs_quant": "#F5A623",  # Orange
    "needs_delete": "#F8E71C",  # Yellow
    "completed": "#7ED321",  # Green
}

# Category display names
CATEGORY_LABELS = {
    "need_download": "Need Download",
    "ongoing_download": "Downloading",
    "failed_download": "Failed",
    "needs_quant": "Needs Quant",
    "needs_delete": "Needs Delete",
    "completed": "Completed",
}

# Category order for consistent display
CATEGORY_ORDER = [
    "need_download",
    "ongoing_download",
    "failed_download",
    "needs_quant",
    "needs_delete",
    "completed",
]


def load_progress_state(state_file: Path | str) -> dict[str, Any]:
    """Load progress state from JSON file.
    
    Args:
        state_file: Path to progress_state.json file
        
    Returns:
        Dictionary mapping species to their state data
        
    Raises:
        FileNotFoundError: If state file doesn't exist
        json.JSONDecodeError: If file is invalid JSON
    """
    state_file = Path(state_file)
    if not state_file.exists():
        raise FileNotFoundError(f"Progress state file not found: {state_file}")
    
    with open(state_file) as f:
        data = json.load(f)
    
    return data


def plot_progress_dashboard(
    progress_state: dict[str, Any] | Path | str,
    *,
    output_path: Path | str | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 300,
) -> Path:
    """Generate multi-panel bar chart visualization of amalgkit progress.
    
    Creates a dashboard showing sample counts in each category (need_download,
    ongoing_download, failed_download, needs_quant, needs_delete, completed)
    for each species, plus an overall summary panel.
    
    Args:
        progress_state: Either a dictionary of progress state data, or path to
            progress_state.json file
        output_path: Path to save PNG file (if None, returns figure without saving)
        figsize: Figure size tuple (auto-calculated if None)
        dpi: Resolution for saved figure
        
    Returns:
        Path to saved file (if output_path provided), or None
        
    Example:
        >>> from metainformant.visualization.amalgkit_visualization import plot_progress_dashboard
        >>> path = plot_progress_dashboard("output/amalgkit/progress_state.json")
    """
    # Load state if path provided
    if isinstance(progress_state, (str, Path)):
        state_data = load_progress_state(progress_state)
    else:
        state_data = progress_state
    
    if not state_data:
        logger.warning("No species data in progress state")
        if output_path:
            # Create empty figure
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, "No species data available", 
                   ha='center', va='center', fontsize=14)
            ax.set_xticks([])
            ax.set_yticks([])
            save_figure(fig, output_path, dpi=dpi)
            plt.close(fig)
            return Path(output_path)
        return None
    
    # Apply publication style
    apply_publication_style()
    
    # Get species list and calculate layout
    species_list = sorted(state_data.keys())
    n_species = len(species_list)
    
    # Calculate number of panels: 1 summary + n_species species panels
    n_panels = n_species + 1
    
    # Calculate optimal grid layout
    # Summary panel spans full width, species panels in grid below
    # Use 3-4 columns for species panels
    ncols = min(4, max(3, int(np.ceil(np.sqrt(n_species)))))
    nrows_summary = 1
    nrows_species = int(np.ceil(n_species / ncols))
    total_rows = nrows_summary + nrows_species
    
    # Calculate figure size
    if figsize is None:
        width = max(12, 3 * ncols)
        height = 4 + 2.5 * nrows_species  # Summary panel + species panels
        figsize = (width, height)
    
    # Create figure with custom layout
    # We'll use subplot2grid for more control
    fig = plt.figure(figsize=figsize, dpi=dpi)
    
    # Create summary panel at top (spans full width)
    gs = fig.add_gridspec(total_rows, ncols, hspace=0.4, wspace=0.3,
                         height_ratios=[1.2] + [1.0] * nrows_species)
    
    ax_summary = fig.add_subplot(gs[0, :])
    
    # Calculate summary totals
    summary_counts = {cat: 0 for cat in CATEGORY_ORDER}
    for species_data in state_data.values():
        for cat in CATEGORY_ORDER:
            if cat in species_data:
                # Handle both list and set formats
                count = len(species_data[cat]) if isinstance(species_data[cat], (list, set)) else 0
                summary_counts[cat] += count
    
    # Plot summary bar chart
    categories = [CATEGORY_LABELS[cat] for cat in CATEGORY_ORDER]
    counts = [summary_counts[cat] for cat in CATEGORY_ORDER]
    colors = [CATEGORY_COLORS[cat] for cat in CATEGORY_ORDER]
    
    bars = ax_summary.bar(categories, counts, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax_summary.set_ylabel("Sample Count", fontsize=11, fontweight='bold')
    ax_summary.set_title("Overall Progress Summary", fontsize=14, fontweight='bold', pad=10)
    ax_summary.grid(axis='y', alpha=0.3, linestyle='--')
    ax_summary.set_axisbelow(True)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax_summary.text(bar.get_x() + bar.get_width()/2., height,
                          f'{int(height)}',
                          ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Rotate x-axis labels if needed
    ax_summary.tick_params(axis='x', rotation=45, labelsize=9)
    
    # Create species panels
    species_axes = []
    for i, species in enumerate(species_list):
        row = 1 + (i // ncols)
        col = i % ncols
        ax = fig.add_subplot(gs[row, col])
        species_axes.append((species, ax))
        
        # Get species data
        species_data = state_data[species]
        
        # Extract counts for each category
        species_counts = []
        for cat in CATEGORY_ORDER:
            if cat in species_data:
                count = len(species_data[cat]) if isinstance(species_data[cat], (list, set)) else 0
            else:
                count = 0
            species_counts.append(count)
        
        # Plot bar chart
        bars = ax.bar(categories, species_counts, color=colors, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
        ax.set_ylabel("Count", fontsize=9)
        
        # Format species name for display
        species_display = species.replace("_", " ").title()
        ax.set_title(species_display, fontsize=10, fontweight='bold', pad=5)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        # Add value labels on bars (only if count > 0 and space allows)
        max_count = max(species_counts) if species_counts else 0
        for bar, count in zip(bars, species_counts):
            height = bar.get_height()
            if height > 0 and max_count < 1000:  # Only label if not too crowded
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}',
                       ha='center', va='bottom', fontsize=7)
        
        # Rotate x-axis labels
        ax.tick_params(axis='x', rotation=45, labelsize=7)
    
    # Hide unused subplots
    total_species_panels = nrows_species * ncols
    for i in range(n_species, total_species_panels):
        row = 1 + (i // ncols)
        col = i % ncols
        ax_empty = fig.add_subplot(gs[row, col])
        ax_empty.set_visible(False)
    
    # Add overall title
    fig.suptitle("Amalgkit RNA-seq Workflow Progress Dashboard", 
                fontsize=16, fontweight='bold', y=0.995)
    
    # Add legend at bottom
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=CATEGORY_COLORS[cat], edgecolor='black', 
              label=CATEGORY_LABELS[cat], alpha=0.8)
        for cat in CATEGORY_ORDER
    ]
    fig.legend(handles=legend_elements, loc='lower center', 
              ncol=len(CATEGORY_ORDER), fontsize=9, frameon=True,
              bbox_to_anchor=(0.5, 0.01))
    
    # Adjust layout (use subplots_adjust instead of tight_layout for custom gridspec)
    plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.12, hspace=0.4, wspace=0.3)
    
    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        save_figure(fig, output_path, dpi=dpi)
        logger.info(f"Progress dashboard saved to {output_path}")
        plt.close(fig)
        return output_path
    
    return fig

