#!/usr/bin/env python3
"""Single-cell RNA-seq analysis example.

This example demonstrates single-cell gene expression analysis using METAINFORMANT's singlecell toolkit.

Usage:
    python examples/singlecell/example_scrna.py

Output:
    output/examples/singlecell/scrna_analysis.json
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate single-cell analysis."""
    # Setup output directory
    output_dir = Path("output/examples/singlecell")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Single-Cell Example ===")

    # Simulated single-cell expression data
    np.random.seed(42)
    n_cells, n_genes = 100, 500

    # Simulate gene expression matrix
    expression_matrix = np.random.negative_binomial(1, 0.8, (n_cells, n_genes))

    # Basic QC metrics
    cells_passing_qc = np.sum(expression_matrix.sum(axis=1) > 100)  # Cells with >100 UMIs
    genes_expressed = np.sum(expression_matrix > 0, axis=0)  # Genes expressed per cell

    results = {
        "cells_total": n_cells,
        "genes_total": n_genes,
        "cells_passing_qc": int(cells_passing_qc),
        "qc_pass_rate": float(cells_passing_qc / n_cells),
        "mean_umis_per_cell": float(np.mean(expression_matrix.sum(axis=1))),
        "genes_expressed_stats": {
            "mean": float(np.mean(genes_expressed)),
            "median": float(np.median(genes_expressed)),
            "max": int(np.max(genes_expressed))
        }
    }

    print(f"✓ Analyzed {n_cells} cells, {n_genes} genes")
    print(f"Cells passing QC: {results['cells_passing_qc']} ({results['qc_pass_rate']:.1%})")

    # Save results
    results_file = output_dir / "scrna_analysis.json"
    io.dump_json({
        "singlecell_analysis": results
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Single-Cell Example Complete ===")

if __name__ == "__main__":
    main()
