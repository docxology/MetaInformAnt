#!/usr/bin/env python3
"""Generate missing PCA scree plot from saved PCA results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from metainformant.gwas.visualization.structure_plots import plot_pca_scree


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate PCA scree plot from saved PCA results")
    parser.add_argument("--results", type=Path, default=Path("output/gwas/pbarbatus/results/pca_results.json"))
    parser.add_argument("--output", type=Path, default=Path("output/gwas/pbarbatus/plots/pca_scree_plot.png"))
    parser.add_argument("--title", type=str, default="PCA Scree Plot\nP. barbatus (n=150 samples, 50,000 variants)")
    args = parser.parse_args()

    with open(args.results) as f:
        pca_data = json.load(f)

    variance_explained = pca_data["explained_variance_ratio"]

    result = plot_pca_scree(
        variance_explained=variance_explained,
        output_path=args.output,
        title=args.title,
    )

    if result["status"] == "success":
        print(f"✓ PCA scree plot generated: {args.output}")
        print(f"  Total variance: {result['total_variance']:.1f}%")
        print(f"  80% variance: PC1-PC{result['pcs_for_80']}")
        print(f"  95% variance: PC1-PC{result['pcs_for_95']}")
    else:
        print(f"✗ Failed: {result.get('error')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
