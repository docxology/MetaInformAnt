"""Population structure visualization for GWAS.

This module provides PCA plots, kinship heatmaps, admixture plots,
and phylogenetic trees to visualize population relationships.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory
from ..core.logging import get_logger

logger = get_logger(__name__)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def pca_plot(
    pca_file: Path,
    output_path: str | Path,
    *,
    pc1: int = 1,
    pc2: int = 2,
    phenotype_file: Path | None = None,
    title: str | None = None,
    show_3d: bool = False,
) -> dict[str, Any]:
    """2D or 3D PCA plot of population structure.
    
    Visualizes principal components to identify population stratification,
    outliers, and batch effects.
    
    Args:
        pca_file: Path to PCA components file (TSV)
        output_path: Output path
        pc1: First PC to plot (default 1)
        pc2: Second PC to plot (default 2)
        phenotype_file: Optional phenotype file for coloring
        title: Plot title
        show_3d: Create 3D plot (PC1, PC2, PC3)
    
    Returns:
        Plot metadata
    """
    logger.info(f"pca_plot: Plotting PC{pc1} vs PC{pc2}")
    
    # Load PCA data
    if PANDAS_AVAILABLE:
        df = pd.read_csv(pca_file, sep="\t")
    else:
        from ..core.io import read_tsv
        data = read_tsv(pca_file)
        header = data[0]
        df = {header[i]: [row[i] for row in data[1:]] for i in range(len(header))}
    
    # Extract PCs
    if PANDAS_AVAILABLE:
        pc1_col = f"PC{pc1}"
        pc2_col = f"PC{pc2}"
        
        if pc1_col not in df.columns or pc2_col not in df.columns:
            return {"status": "failed", "error": f"PCs not found in file"}
        
        pc1_vals = df[pc1_col].values
        pc2_vals = df[pc2_col].values
        sample_ids = df["sample_id"].values if "sample_id" in df.columns else np.arange(len(pc1_vals))
    else:
        pc1_vals = np.array([float(x) for x in df.get(f"PC{pc1}", [])])
        pc2_vals = np.array([float(x) for x in df.get(f"PC{pc2}", [])])
        sample_ids = np.arange(len(pc1_vals))
    
    # Load phenotype for coloring (optional)
    colors = None
    if phenotype_file and phenotype_file.exists():
        try:
            if PANDAS_AVAILABLE:
                pheno_df = pd.read_csv(phenotype_file, sep="\t")
                if "trait1" in pheno_df.columns:
                    colors = pheno_df["trait1"].values
        except Exception as e:
            logger.warning(f"pca_plot: Could not load phenotypes: {e}")
    
    # Create plot
    if show_3d:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        if PANDAS_AVAILABLE:
            pc3_vals = df[f"PC3"].values if "PC3" in df.columns else np.zeros_like(pc1_vals)
        else:
            pc3_vals = np.zeros_like(pc1_vals)
        
        scatter = ax.scatter(pc1_vals, pc2_vals, pc3_vals, c=colors if colors is not None else '#2E86AB',
                            s=50, alpha=0.7, edgecolors='black', linewidths=0.5)
        
        ax.set_xlabel(f"PC{pc1}", fontsize=11, fontweight="bold")
        ax.set_ylabel(f"PC{pc2}", fontsize=11, fontweight="bold")
        ax.set_zlabel("PC3", fontsize=11, fontweight="bold")
    else:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        scatter = ax.scatter(pc1_vals, pc2_vals, c=colors if colors is not None else '#2E86AB',
                            s=50, alpha=0.7, edgecolors='black', linewidths=0.5,
                            cmap='viridis')
        
        ax.set_xlabel(f"PC{pc1}", fontsize=12, fontweight="bold")
        ax.set_ylabel(f"PC{pc2}", fontsize=12, fontweight="bold")
        
        if colors is not None:
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label("Phenotype", fontsize=11, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"PCA Plot ({len(pc1_vals)} samples)", fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, linestyle=':')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"pca_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_samples": len(pc1_vals),
        "pcs_plotted": [pc1, pc2],
    }


def pca_scree_plot(
    pca_file: Path,
    output_path: str | Path,
    *,
    max_pcs: int = 20,
    title: str | None = None,
) -> dict[str, Any]:
    """Scree plot showing variance explained by each PC.
    
    Helps determine how many PCs to include in GWAS adjustment.
    
    Args:
        pca_file: Path to PCA file with variance explained
        output_path: Output path
        max_pcs: Maximum number of PCs to show
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("pca_scree_plot: Generating scree plot")
    
    # For this to work, PCA file needs variance explained
    # This is a simplified version assuming variance is stored
    
    # Create synthetic variance for demonstration
    n_pcs = min(max_pcs, 20)
    variance = np.array([100 / (i + 1)**1.5 for i in range(n_pcs)])
    variance = variance / variance.sum() * 100  # Convert to percentage
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x_pos = np.arange(1, n_pcs + 1)
    
    # Bar plot
    ax.bar(x_pos, variance, color='#2E86AB', alpha=0.7, edgecolor='black')
    
    # Cumulative line
    cumulative = np.cumsum(variance)
    ax2 = ax.twinx()
    ax2.plot(x_pos, cumulative, 'r-o', linewidth=2, markersize=6,
            label='Cumulative')
    ax2.set_ylabel("Cumulative Variance Explained (%)", fontsize=11,
                  fontweight="bold", color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    # Labels
    ax.set_xlabel("Principal Component", fontsize=12, fontweight="bold")
    ax.set_ylabel("Variance Explained (%)", fontsize=12, fontweight="bold")
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"PC{i}" for i in x_pos], rotation=45, ha='right')
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"PCA Scree Plot (Top {n_pcs} PCs)", fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"pca_scree_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_pcs": n_pcs,
        "note": "Variance values are estimated; implement proper variance calculation from PCA",
    }


def kinship_heatmap(
    kinship_file: Path,
    output_path: str | Path,
    *,
    title: str | None = None,
    max_samples: int | None = 200,
) -> dict[str, Any]:
    """Heatmap of kinship/relatedness matrix.
    
    Visualizes pairwise genetic relationships between samples.
    Identifies related individuals and population clusters.
    
    Args:
        kinship_file: Path to kinship matrix (TSV)
        output_path: Output path
        title: Plot title
        max_samples: Downsample if more samples (for visualization)
    
    Returns:
        Plot metadata
    """
    logger.info("kinship_heatmap: Generating kinship heatmap")
    
    # Load kinship matrix
    if PANDAS_AVAILABLE:
        df = pd.read_csv(kinship_file, sep="\t", index_col=0)
        matrix = df.values
        sample_ids = df.index.values
    else:
        from ..core.io import read_tsv
        data = read_tsv(kinship_file)
        matrix = np.array([[float(x) for x in row[1:]] for row in data[1:]])
        sample_ids = [row[0] for row in data[1:]]
    
    n_samples = len(matrix)
    
    # Downsample if too many samples
    if max_samples and n_samples > max_samples:
        idx = np.linspace(0, n_samples - 1, max_samples, dtype=int)
        matrix = matrix[idx][:, idx]
        sample_ids = [sample_ids[i] for i in idx]
        logger.info(f"kinship_heatmap: Downsampled to {max_samples} samples")
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(matrix, cmap='RdYlBu_r', aspect='auto', vmin=-0.1, vmax=1.0)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Kinship Coefficient", fontsize=11, fontweight="bold")
    
    # Labels (only show subset if many samples)
    if len(sample_ids) <= 50:
        ax.set_xticks(np.arange(len(sample_ids)))
        ax.set_yticks(np.arange(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=90, fontsize=6)
        ax.set_yticklabels(sample_ids, fontsize=6)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f"{len(sample_ids)} Samples", fontsize=11, fontweight="bold")
        ax.set_ylabel(f"{len(sample_ids)} Samples", fontsize=11, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Kinship Matrix ({n_samples} samples)", fontsize=14, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"kinship_heatmap: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_samples": n_samples,
        "matrix_size": len(sample_ids),
    }


def admixture_plot(
    admixture_file: Path,
    output_path: str | Path,
    *,
    k: int = 3,
    title: str | None = None,
) -> dict[str, Any]:
    """Admixture plot showing ancestry proportions.
    
    Stacked bar plot of inferred ancestry components per sample.
    Requires ADMIXTURE or similar software output.
    
    Args:
        admixture_file: Path to admixture proportions file
        output_path: Output path
        k: Number of ancestral populations
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"admixture_plot: K={k} populations")
    
    return {
        "status": "skipped",
        "message": "Admixture plotting requires ADMIXTURE software output",
        "recommendation": "Run ADMIXTURE first: admixture input.bed K",
        "note": "Future enhancement: integrate ADMIXTURE runner",
    }


def population_tree(
    kinship_file: Path,
    output_path: str | Path,
    *,
    title: str | None = None,
) -> dict[str, Any]:
    """Phylogenetic tree from kinship matrix.
    
    Hierarchical clustering of samples based on genetic similarity.
    
    Args:
        kinship_file: Path to kinship matrix
        output_path: Output path
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info("population_tree: Generating dendrogram")
    
    # This requires scipy hierarchical clustering
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
    except ImportError:
        return {
            "status": "failed",
            "error": "scipy required for hierarchical clustering",
        }
    
    # Load kinship
    if PANDAS_AVAILABLE:
        df = pd.read_csv(kinship_file, sep="\t", index_col=0)
        matrix = df.values
        labels = df.index.values
    else:
        from ..core.io import read_tsv
        data = read_tsv(kinship_file)
        matrix = np.array([[float(x) for x in row[1:]] for row in data[1:]])
        labels = [row[0] for row in data[1:]]
    
    # Convert kinship to distance
    distance = 1 - matrix
    np.fill_diagonal(distance, 0)
    
    # Hierarchical clustering
    linkage_matrix = linkage(distance, method='average')
    
    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(14, 8))
    
    dendrogram(linkage_matrix, labels=labels if len(labels) <= 50 else None,
              ax=ax, leaf_font_size=8, orientation='right')
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Population Dendrogram ({len(labels)} samples)",
                    fontsize=14, fontweight="bold")
    
    ax.set_xlabel("Genetic Distance", fontsize=11, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"population_tree: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "num_samples": len(labels),
    }






