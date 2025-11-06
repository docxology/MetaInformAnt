#!/usr/bin/env python3
"""Single-cell simulation script.

This script generates synthetic single-cell RNA-seq data including count matrices,
cell type annotations, and trajectory data.

Usage:
    python3 scripts/simulation/simulate_singlecell.py --type counts --n-cells 1000 --n-genes 2000
    python3 scripts/simulation/simulate_singlecell.py --type celltypes --n-cells 500 --n-genes 1500 --n-types 5
    python3 scripts/simulation/simulate_singlecell.py --type trajectory --n-cells 800 --n-genes 2000
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.rna import simulate_counts_negative_binomial

logger = logging.get_logger(__name__)


def simulate_counts(
    output_dir: Path,
    n_cells: int,
    n_genes: int,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate basic single-cell count matrix.
    
    Args:
        output_dir: Output directory for results
        n_cells: Number of cells
        n_genes: Number of genes
        mean_expression: Mean expression level
        dispersion: Dispersion parameter
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_cells, min_val=1, name="n_cells")
    validation.validate_range(n_genes, min_val=1, name="n_genes")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")
    
    logger.info(f"Generating scRNA-seq counts: {n_cells} cells x {n_genes} genes")
    rng = random.Random(seed)
    
    # Transpose: genes x cells (as expected by simulation function)
    counts = simulate_counts_negative_binomial(
        num_genes=n_genes,
        num_samples=n_cells,
        mean_expression=mean_expression,
        dispersion=dispersion,
        rng=rng,
    )
    
    # Create DataFrame (genes x cells)
    df = pd.DataFrame(counts, columns=[f"cell_{i:05d}" for i in range(n_cells)])
    df.index = [f"gene_{i:05d}" for i in range(n_genes)]
    
    csv_file = output_dir / "singlecell_counts.csv"
    df.to_csv(csv_file)
    
    logger.info(f"Single-cell counts saved to {csv_file}")
    
    return {
        "type": "counts",
        "n_cells": n_cells,
        "n_genes": n_genes,
        "output_file": str(csv_file),
    }


def simulate_celltypes(
    output_dir: Path,
    n_cells: int,
    n_genes: int,
    n_types: int,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate single-cell data with cell type annotations.
    
    Args:
        output_dir: Output directory for results
        n_cells: Number of cells
        n_genes: Number of genes
        n_types: Number of cell types
        mean_expression: Mean expression level
        dispersion: Dispersion parameter
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_cells, min_val=1, name="n_cells")
    validation.validate_range(n_genes, min_val=1, name="n_genes")
    validation.validate_range(n_types, min_val=1, max_val=n_cells, name="n_types")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")
    
    logger.info(f"Generating cell type data: {n_cells} cells, {n_types} types")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    # Generate marker genes per cell type
    n_markers = 50
    marker_genes = set()
    type_markers = {}
    
    for type_idx in range(n_types):
        markers = rng.sample(range(n_genes), min(n_markers, n_genes))
        type_markers[f"type_{type_idx:02d}"] = markers
        marker_genes.update(markers)
    
    # Generate counts with cell type-specific expression
    counts = []
    cell_metadata = []
    
    cells_per_type = n_cells // n_types
    
    for type_idx in range(n_types):
        type_name = f"type_{type_idx:02d}"
        markers = type_markers[type_name]
        
        for cell_idx in range(cells_per_type):
            cell_id = f"cell_{type_idx}_{cell_idx:05d}"
            
            # Generate expression
            cell_counts = []
            for gene_idx in range(n_genes):
                if gene_idx in markers:
                    # Higher expression for marker genes
                    mean = mean_expression * 3.0
                else:
                    mean = mean_expression
                
                count = simulate_counts_negative_binomial(
                    1, 1, mean_expression=mean, dispersion=dispersion, rng=rng
                )[0][0]
                cell_counts.append(count)
            
            counts.append(cell_counts)
            cell_metadata.append({
                "cell_id": cell_id,
                "cell_type": type_name,
            })
    
    # Add remaining cells
    for cell_idx in range(n_cells - len(counts)):
        type_idx = rng.randint(0, n_types - 1)
        type_name = f"type_{type_idx:02d}"
        markers = type_markers[type_name]
        
        cell_id = f"cell_{type_idx}_{cell_idx:05d}"
        cell_counts = []
        for gene_idx in range(n_genes):
            if gene_idx in markers:
                mean = mean_expression * 3.0
            else:
                mean = mean_expression
            
            count = simulate_counts_negative_binomial(
                1, 1, mean_expression=mean, dispersion=dispersion, rng=rng
            )[0][0]
            cell_counts.append(count)
        
        counts.append(cell_counts)
        cell_metadata.append({
            "cell_id": cell_id,
            "cell_type": type_name,
        })
    
    # Create count matrix (cells x genes, then transpose)
    count_matrix = np.array(counts).T  # genes x cells
    df = pd.DataFrame(count_matrix, columns=[m["cell_id"] for m in cell_metadata])
    df.index = [f"gene_{i:05d}" for i in range(n_genes)]
    
    csv_file = output_dir / "celltype_counts.csv"
    df.to_csv(csv_file)
    
    # Save metadata
    metadata_df = pd.DataFrame(cell_metadata)
    metadata_file = output_dir / "cell_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    logger.info(f"Cell type data saved to {csv_file}")
    
    return {
        "type": "celltypes",
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_types": n_types,
        "output_file": str(csv_file),
        "metadata_file": str(metadata_file),
    }


def simulate_trajectory(
    output_dir: Path,
    n_cells: int,
    n_genes: int,
    n_states: int,
    mean_expression: float,
    dispersion: float,
    seed: int,
) -> dict:
    """Simulate trajectory/pseudotime data.
    
    Args:
        output_dir: Output directory for results
        n_cells: Number of cells
        n_genes: Number of genes
        n_states: Number of trajectory states
        mean_expression: Mean expression level
        dispersion: Dispersion parameter
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with simulation results and metadata
        
    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_cells, min_val=1, name="n_cells")
    validation.validate_range(n_genes, min_val=1, name="n_genes")
    validation.validate_range(n_states, min_val=1, max_val=n_cells, name="n_states")
    validation.validate_range(mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(dispersion, min_val=0.0, name="dispersion")
    
    logger.info(f"Generating trajectory data: {n_cells} cells, {n_states} states")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    # Generate pseudotime
    pseudotimes = np.random.uniform(0, 1, n_cells)
    pseudotimes.sort()
    
    # Assign cells to states based on pseudotime
    state_assignments = (pseudotimes * n_states).astype(int)
    state_assignments = np.clip(state_assignments, 0, n_states - 1)
    
    # Generate trajectory-specific genes
    n_trajectory_genes = 100
    trajectory_genes = rng.sample(range(n_genes), n_trajectory_genes)
    
    # Generate counts with trajectory-dependent expression
    counts = []
    cell_metadata = []
    
    for cell_idx, (pseudotime, state) in enumerate(zip(pseudotimes, state_assignments)):
        cell_id = f"cell_{cell_idx:05d}"
        
        cell_counts = []
        for gene_idx in range(n_genes):
            if gene_idx in trajectory_genes:
                # Expression changes along trajectory
                mean = mean_expression * (1 + pseudotime * 2)
            else:
                mean = mean_expression
            
            count = simulate_counts_negative_binomial(
                1, 1, mean_expression=mean, dispersion=dispersion, rng=rng
            )[0][0]
            cell_counts.append(count)
        
        counts.append(cell_counts)
        cell_metadata.append({
            "cell_id": cell_id,
            "pseudotime": float(pseudotime),
            "state": int(state),
        })
    
    # Create count matrix
    count_matrix = np.array(counts).T  # genes x cells
    df = pd.DataFrame(count_matrix, columns=[m["cell_id"] for m in cell_metadata])
    df.index = [f"gene_{i:05d}" for i in range(n_genes)]
    
    csv_file = output_dir / "trajectory_counts.csv"
    df.to_csv(csv_file)
    
    # Save metadata
    metadata_df = pd.DataFrame(cell_metadata)
    metadata_file = output_dir / "trajectory_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    logger.info(f"Trajectory data saved to {csv_file}")
    
    return {
        "type": "trajectory",
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_states": n_states,
        "output_file": str(csv_file),
        "metadata_file": str(metadata_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Single-cell simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate basic counts
  %(prog)s --type counts --n-cells 1000 --n-genes 2000

  # Simulate with cell types
  %(prog)s --type celltypes --n-cells 500 --n-genes 1500 --n-types 5

  # Simulate trajectory
  %(prog)s --type trajectory --n-cells 800 --n-genes 2000 --n-states 10
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["counts", "celltypes", "trajectory"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/singlecell"),
        help="Output directory (default: output/simulation/singlecell)",
    )
    parser.add_argument("--n-cells", type=int, default=1000, help="Number of cells")
    parser.add_argument("--n-genes", type=int, default=2000, help="Number of genes")
    parser.add_argument("--mean-expression", type=float, default=5.0, help="Mean expression level")
    parser.add_argument("--dispersion", type=float, default=0.5, help="Dispersion parameter")
    parser.add_argument("--n-types", type=int, default=5, help="Number of cell types (celltypes type)")
    parser.add_argument("--n-states", type=int, default=10, help="Number of trajectory states (trajectory type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    # Validate output directory
    output_dir = paths.ensure_directory(args.output)
    
    # Validate common parameters
    validation.validate_range(args.n_cells, min_val=1, name="n_cells")
    validation.validate_range(args.n_genes, min_val=1, name="n_genes")
    validation.validate_range(args.mean_expression, min_val=0.0, name="mean_expression")
    validation.validate_range(args.dispersion, min_val=0.0, name="dispersion")
    if hasattr(args, "n_types"):
        validation.validate_range(args.n_types, min_val=1, max_val=args.n_cells, name="n_types")
    if hasattr(args, "n_states"):
        validation.validate_range(args.n_states, min_val=1, max_val=args.n_cells, name="n_states")
    
    try:
        if args.type == "counts":
            results = simulate_counts(
                output_dir,
                args.n_cells,
                args.n_genes,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )
        elif args.type == "celltypes":
            results = simulate_celltypes(
                output_dir,
                args.n_cells,
                args.n_genes,
                args.n_types,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )
        elif args.type == "trajectory":
            results = simulate_trajectory(
                output_dir,
                args.n_cells,
                args.n_genes,
                args.n_states,
                args.mean_expression,
                args.dispersion,
                args.seed,
            )
        
        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

