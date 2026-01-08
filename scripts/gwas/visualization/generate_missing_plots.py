#!/usr/bin/env python3
"""Generate missing PCA scree plot from saved PCA results."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent))

import json
import numpy as np
from visualizations import plot_pca_scree

# Load PCA results
pca_file = Path('output/gwas/pbarbatus/results/pca_results.json')
with open(pca_file) as f:
    pca_data = json.load(f)

variance_explained = np.array(pca_data['explained_variance_ratio'])

# Generate scree plot
output_path = Path('output/gwas/pbarbatus/plots/pca_scree_plot.png')
result = plot_pca_scree(
    variance_explained=variance_explained,
    output_path=output_path,
    title='PCA Scree Plot\nP. barbatus (n=150 samples, 50,000 variants)',
)

if result['status'] == 'success':
    print(f"✓ PCA scree plot generated: {output_path}")
    print(f"  Total variance: {result['total_variance']:.1f}%")
    print(f"  80% variance: PC1-PC{result['pcs_for_80']}")
    print(f"  95% variance: PC1-PC{result['pcs_for_95']}")
else:
    print(f"✗ Failed: {result.get('error')}")





