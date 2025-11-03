#!/usr/bin/env python3
"""Enhanced visualization functions for GWAS population structure analysis.

Creates comprehensive visualizations for:
- Kinship matrix heatmaps
- PCA scatter plots (multiple PC combinations)
- PCA scree plots (variance explained)
- Population structure summaries
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Any


def plot_kinship_heatmap(
    kinship_matrix: np.ndarray,
    output_path: Path,
    sample_ids: list[str] | None = None,
    title: str = "Kinship Matrix",
) -> dict[str, Any]:
    """Create a heatmap visualization of the kinship matrix.
    
    Args:
        kinship_matrix: NxN kinship/relatedness matrix
        output_path: Path to save the plot
        sample_ids: Optional sample IDs for axis labels
        title: Plot title
    
    Returns:
        Dictionary with status and statistics
    """
    try:
        n_samples = len(kinship_matrix)
        
        # Compute statistics
        # Exclude diagonal for relatedness statistics
        mask = ~np.eye(n_samples, dtype=bool)
        off_diagonal = kinship_matrix[mask]
        
        mean_relatedness = np.mean(off_diagonal)
        std_relatedness = np.std(off_diagonal)
        min_relatedness = np.min(off_diagonal)
        max_relatedness = np.max(off_diagonal)
        median_relatedness = np.median(off_diagonal)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Create heatmap
        im = ax.imshow(kinship_matrix, cmap='RdYlBu_r', aspect='auto',
                      vmin=min_relatedness, vmax=max_relatedness)
        
        # Add colorbar with label
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20, fontsize=11)
        
        # Set title with statistics
        stats_text = (f'Mean: {mean_relatedness:.4f}, SD: {std_relatedness:.4f}\n'
                     f'Range: [{min_relatedness:.4f}, {max_relatedness:.4f}], '
                     f'Median: {median_relatedness:.4f}')
        ax.set_title(f'{title}\n{stats_text}', fontsize=13, pad=15)
        
        # Configure axes
        ax.set_xlabel('Sample Index', fontsize=11)
        ax.set_ylabel('Sample Index', fontsize=11)
        
        # Add grid for better readability
        if n_samples <= 50:
            # Show tick labels for small matrices
            ax.set_xticks(np.arange(n_samples))
            ax.set_yticks(np.arange(n_samples))
            if sample_ids:
                ax.set_xticklabels(sample_ids, rotation=90, fontsize=6)
                ax.set_yticklabels(sample_ids, fontsize=6)
        else:
            # Show fewer ticks for large matrices
            step = max(1, n_samples // 20)
            ticks = np.arange(0, n_samples, step)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            if sample_ids:
                ax.set_xticklabels([sample_ids[i] for i in ticks], rotation=90, fontsize=8)
                ax.set_yticklabels([sample_ids[i] for i in ticks], fontsize=8)
        
        # Add text annotation for highly related pairs
        high_related = np.where(kinship_matrix > mean_relatedness + 2*std_relatedness)
        n_high_related = len(high_related[0])
        if n_high_related > 0 and n_high_related < 100:
            for i, j in zip(high_related[0], high_related[1]):
                if i != j:  # Skip diagonal
                    ax.plot(j, i, 'r*', markersize=3, alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return {
            'status': 'success',
            'output_path': str(output_path),
            'n_samples': n_samples,
            'mean_relatedness': float(mean_relatedness),
            'std_relatedness': float(std_relatedness),
            'min_relatedness': float(min_relatedness),
            'max_relatedness': float(max_relatedness),
            'median_relatedness': float(median_relatedness),
        }
        
    except Exception as e:
        return {
            'status': 'failed',
            'error': str(e),
        }


def plot_pca_scatter(
    pca_components: np.ndarray,
    variance_explained: np.ndarray,
    output_path: Path,
    pc_x: int = 0,
    pc_y: int = 1,
    phenotypes: dict[str, list] | None = None,
    color_by: str | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Create PCA scatter plot.
    
    Args:
        pca_components: PC matrix (samples x components)
        variance_explained: Variance explained by each PC
        output_path: Path to save the plot
        pc_x: PC index for X-axis (0-based)
        pc_y: PC index for Y-axis (0-based)
        phenotypes: Optional phenotype data for coloring
        color_by: Phenotype name to color points by
        title: Optional plot title
    
    Returns:
        Dictionary with status
    """
    try:
        n_samples, n_components = pca_components.shape
        
        if pc_x >= n_components or pc_y >= n_components:
            return {
                'status': 'failed',
                'error': f'PC indices out of range (max: {n_components-1})',
            }
        
        # Extract PCs
        pc_x_vals = pca_components[:, pc_x]
        pc_y_vals = pca_components[:, pc_y]
        
        var_x = variance_explained[pc_x] * 100
        var_y = variance_explained[pc_y] * 100
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Determine coloring
        if phenotypes and color_by and color_by in phenotypes:
            pheno_vals = phenotypes[color_by]
            # Check if categorical or continuous
            unique_vals = np.unique(pheno_vals)
            if len(unique_vals) <= 10:  # Categorical
                scatter = ax.scatter(pc_x_vals, pc_y_vals, c=pheno_vals,
                                   cmap='tab10', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(color_by, rotation=270, labelpad=20)
            else:  # Continuous
                scatter = ax.scatter(pc_x_vals, pc_y_vals, c=pheno_vals,
                                   cmap='viridis', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(color_by, rotation=270, labelpad=20)
        else:
            # Default: no coloring
            ax.scatter(pc_x_vals, pc_y_vals, c='steelblue', s=50, alpha=0.7,
                      edgecolors='black', linewidth=0.5)
        
        # Labels and title
        ax.set_xlabel(f'PC{pc_x+1} ({var_x:.2f}% variance)', fontsize=12)
        ax.set_ylabel(f'PC{pc_y+1} ({var_y:.2f}% variance)', fontsize=12)
        
        if title:
            plot_title = title
        else:
            plot_title = f'PCA: PC{pc_x+1} vs PC{pc_y+1}'
        
        # Add statistics to title
        stats_text = f'n={n_samples} samples'
        ax.set_title(f'{plot_title}\n{stats_text}', fontsize=13, pad=15)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Add reference lines at zero
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return {
            'status': 'success',
            'output_path': str(output_path),
            'pc_x': pc_x + 1,
            'pc_y': pc_y + 1,
            'var_x': float(var_x),
            'var_y': float(var_y),
        }
        
    except Exception as e:
        return {
            'status': 'failed',
            'error': str(e),
        }


def plot_pca_scree(
    variance_explained: np.ndarray,
    output_path: Path,
    title: str = "PCA Scree Plot",
) -> dict[str, Any]:
    """Create scree plot showing variance explained by each PC.
    
    Args:
        variance_explained: Variance explained by each PC (as proportions)
        output_path: Path to save the plot
        title: Plot title
    
    Returns:
        Dictionary with status
    """
    try:
        n_components = len(variance_explained)
        variance_pct = variance_explained * 100
        cumulative_var = np.cumsum(variance_pct)
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Plot 1: Individual variance explained
        pcs = np.arange(1, n_components + 1)
        ax1.bar(pcs, variance_pct, color='steelblue', alpha=0.7, edgecolor='black')
        ax1.plot(pcs, variance_pct, 'ro-', linewidth=2, markersize=6, label='Variance per PC')
        ax1.set_xlabel('Principal Component', fontsize=11)
        ax1.set_ylabel('Variance Explained (%)', fontsize=11)
        ax1.set_title('Variance Explained per PC', fontsize=12)
        ax1.legend()
        ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax1.set_xticks(pcs)
        
        # Plot 2: Cumulative variance
        ax2.plot(pcs, cumulative_var, 'go-', linewidth=2, markersize=6, label='Cumulative variance')
        ax2.axhline(y=80, color='red', linestyle='--', linewidth=1, alpha=0.7, label='80% threshold')
        ax2.axhline(y=95, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='95% threshold')
        ax2.set_xlabel('Principal Component', fontsize=11)
        ax2.set_ylabel('Cumulative Variance Explained (%)', fontsize=11)
        ax2.set_title('Cumulative Variance Explained', fontsize=12)
        ax2.legend()
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.set_xticks(pcs)
        ax2.set_ylim([0, 105])
        
        # Add overall title with statistics
        total_var = cumulative_var[-1]
        pcs_80 = np.searchsorted(cumulative_var, 80) + 1
        pcs_95 = np.searchsorted(cumulative_var, 95) + 1
        stats_text = f'Total: {total_var:.1f}% | 80% reached by PC{pcs_80} | 95% reached by PC{pcs_95}'
        fig.suptitle(f'{title}\n{stats_text}', fontsize=13, y=1.02)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return {
            'status': 'success',
            'output_path': str(output_path),
            'n_components': n_components,
            'total_variance': float(total_var),
            'pcs_for_80': int(pcs_80),
            'pcs_for_95': int(pcs_95),
        }
        
    except Exception as e:
        return {
            'status': 'failed',
            'error': str(e),
        }


def plot_pca_loadings(
    loadings: np.ndarray,
    output_path: Path,
    pc_index: int = 0,
    top_n: int = 20,
    title: str | None = None,
) -> dict[str, Any]:
    """Plot top variant loadings for a specific PC.
    
    Args:
        loadings: Variant loadings matrix (variants x PCs)
        output_path: Path to save the plot
        pc_index: Which PC to plot (0-based)
        top_n: Number of top loadings to show
        title: Optional plot title
    
    Returns:
        Dictionary with status
    """
    try:
        n_variants, n_components = loadings.shape
        
        if pc_index >= n_components:
            return {
                'status': 'failed',
                'error': f'PC index out of range (max: {n_components-1})',
            }
        
        # Get loadings for this PC
        pc_loadings = loadings[:, pc_index]
        
        # Get top N by absolute value
        abs_loadings = np.abs(pc_loadings)
        top_indices = np.argsort(abs_loadings)[-top_n:][::-1]
        top_values = pc_loadings[top_indices]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create horizontal bar plot
        y_pos = np.arange(len(top_indices))
        colors = ['red' if v < 0 else 'blue' for v in top_values]
        ax.barh(y_pos, top_values, color=colors, alpha=0.7, edgecolor='black')
        
        # Labels
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f'Variant {i+1}' for i in top_indices], fontsize=9)
        ax.set_xlabel(f'PC{pc_index+1} Loading', fontsize=11)
        ax.set_ylabel('Variant', fontsize=11)
        
        if title:
            plot_title = title
        else:
            plot_title = f'Top {top_n} Variant Loadings for PC{pc_index+1}'
        ax.set_title(plot_title, fontsize=12, pad=15)
        
        # Add vertical line at zero
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--', axis='x')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return {
            'status': 'success',
            'output_path': str(output_path),
            'pc_index': pc_index + 1,
            'top_n': top_n,
        }
        
    except Exception as e:
        return {
            'status': 'failed',
            'error': str(e),
        }


if __name__ == '__main__':
    print("Enhanced GWAS visualization functions")
    print("Functions: plot_kinship_heatmap, plot_pca_scatter, plot_pca_scree, plot_pca_loadings")

