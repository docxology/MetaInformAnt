#!/usr/bin/env python3
"""Single-cell genomics workflow orchestrator.

This script provides comprehensive orchestration for single-cell analysis workflows,
including data loading, QC, normalization, dimensionality reduction, and clustering.

Usage:
    python3 scripts/singlecell/run_singlecell_analysis.py --input counts.csv --output output/singlecell/results
    python3 scripts/singlecell/run_singlecell_analysis.py --input counts.csv --qc --normalize --cluster
    python3 scripts/singlecell/run_singlecell_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json, read_csv
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Single-cell genomics workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline with QC and clustering
  %(prog)s --input counts.csv --qc --normalize --cluster --output output/singlecell/full

  # QC and dimensionality reduction only
  %(prog)s --input counts.csv --qc --normalize --pca --n-components 50

  # Clustering with marker genes
  %(prog)s --input counts.csv --cluster --find-markers
        """
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input count matrix (CSV/TSV format, cells x genes)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/singlecell"),
        help="Output directory (default: output/singlecell)",
    )
    parser.add_argument(
        "--qc",
        action="store_true",
        help="Perform quality control and filtering",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize counts",
    )
    parser.add_argument(
        "--pca",
        action="store_true",
        help="Perform PCA dimensionality reduction",
    )
    parser.add_argument(
        "--umap",
        action="store_true",
        help="Perform UMAP dimensionality reduction",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="Perform clustering (Leiden algorithm)",
    )
    parser.add_argument(
        "--find-markers",
        action="store_true",
        help="Find marker genes for clusters",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=50,
        help="Number of PCA components (default: 50)",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Clustering resolution (default: 0.5)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )
    return parser.parse_args()


def run_workflow(args):
    """Execute single-cell analysis workflow."""
    logger.info("Starting single-cell analysis workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load count matrix
    logger.info(f"Loading count matrix from {args.input}")
    try:
        df = read_csv(args.input)
        logger.info(f"Loaded count matrix: {df.shape[0]} cells, {df.shape[1]} genes")
        
        # Convert to SingleCellData format
        from metainformant.singlecell import SingleCellData
        import numpy as np
        
        # Assuming genes are columns, cells are rows
        counts = df.values
        gene_names = list(df.columns)
        cell_names = list(df.index) if hasattr(df.index, 'tolist') else [f"cell_{i}" for i in range(len(df))]
        
        data = SingleCellData(
            counts=counts,
            gene_names=gene_names,
            cell_names=cell_names
        )
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise
    
    workflow_results = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "n_cells": data.n_cells,
        "n_genes": data.n_genes,
        "analyses": {},
    }
    
    # Quality control
    if args.qc:
        try:
            logger.info("Performing quality control...")
            from metainformant.singlecell import calculate_qc_metrics, filter_cells
            
            data = calculate_qc_metrics(data)
            data = filter_cells(data, min_genes=200, max_pct_mt=20.0)
            
            workflow_results["analyses"]["qc"] = {
                "n_cells_after_filtering": data.n_cells,
                "n_genes_after_filtering": data.n_genes,
            }
            logger.info(f"After QC: {data.n_cells} cells, {data.n_genes} genes")
        except Exception as e:
            logger.error(f"QC failed: {e}", exc_info=True)
            workflow_results["analyses"]["qc"] = {"error": str(e)}
    
    # Normalization
    if args.normalize:
        try:
            logger.info("Normalizing counts...")
            from metainformant.singlecell import normalize_counts, log_transform
            
            data = normalize_counts(data)
            data = log_transform(data)
            
            workflow_results["analyses"]["normalization"] = {"status": "completed"}
        except Exception as e:
            logger.error(f"Normalization failed: {e}", exc_info=True)
            workflow_results["analyses"]["normalization"] = {"error": str(e)}
    
    # PCA
    if args.pca:
        try:
            logger.info(f"Computing PCA with {args.n_components} components...")
            from metainformant.singlecell import select_hvgs, compute_pca
            
            data = select_hvgs(data, n_top_genes=2000)
            data = compute_pca(data, n_components=args.n_components)
            
            workflow_results["analyses"]["pca"] = {
                "n_components": args.n_components,
                "status": "completed",
            }
        except Exception as e:
            logger.error(f"PCA failed: {e}", exc_info=True)
            workflow_results["analyses"]["pca"] = {"error": str(e)}
    
    # UMAP
    if args.umap:
        try:
            logger.info("Computing UMAP...")
            from metainformant.singlecell import compute_neighbors, compute_umap
            
            if not hasattr(data, 'pca') or data.pca is None:
                # Need PCA first
                data = select_hvgs(data, n_top_genes=2000)
                data = compute_pca(data, n_components=50)
            
            data = compute_neighbors(data, n_neighbors=15)
            data = compute_umap(data, min_dist=0.1)
            
            workflow_results["analyses"]["umap"] = {"status": "completed"}
        except Exception as e:
            logger.error(f"UMAP failed: {e}", exc_info=True)
            workflow_results["analyses"]["umap"] = {"error": str(e)}
    
    # Clustering
    if args.cluster:
        try:
            logger.info(f"Performing clustering (resolution={args.resolution})...")
            from metainformant.singlecell import leiden_clustering
            
            data = leiden_clustering(data, resolution=args.resolution, random_state=42)
            
            # Count clusters
            if hasattr(data, 'obs') and 'cluster' in data.obs:
                n_clusters = len(data.obs['cluster'].unique())
            else:
                n_clusters = "unknown"
            
            workflow_results["analyses"]["clustering"] = {
                "resolution": args.resolution,
                "n_clusters": n_clusters,
                "status": "completed",
            }
            logger.info(f"Identified {n_clusters} clusters")
        except Exception as e:
            logger.error(f"Clustering failed: {e}", exc_info=True)
            workflow_results["analyses"]["clustering"] = {"error": str(e)}
    
    # Marker genes
    if args.find_markers:
        try:
            logger.info("Finding marker genes...")
            from metainformant.singlecell import find_marker_genes
            
            markers = find_marker_genes(data, cluster_key="cluster")
            
            workflow_results["analyses"]["markers"] = {
                "n_markers": len(markers) if isinstance(markers, (list, dict)) else "unknown",
                "status": "completed",
            }
        except Exception as e:
            logger.error(f"Marker finding failed: {e}", exc_info=True)
            workflow_results["analyses"]["markers"] = {"error": str(e)}
    
    # Save summary
    summary_file = output_dir / "workflow_summary.json"
    dump_json(workflow_results, summary_file, indent=2)
    logger.info(f"Workflow summary saved to {summary_file}")
    
    logger.info("Workflow complete")


def main():
    """Main entry point."""
    args = parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        run_workflow(args)
        return 0
    except Exception as e:
        logger.error(f"Workflow failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
