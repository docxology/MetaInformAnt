#!/usr/bin/env python
"""eQTL Analysis Pipeline Demo for Apis mellifera.

Demonstrates the full eQTL analysis workflow using synthetic data
modeled after the available RNA-seq quantification data.

Usage:
    uv run python scripts/eqtl/run_eqtl_demo.py
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.eqtl.synthetic import create_synthetic_data
from metainformant.gwas.finemapping.eqtl import (
    cis_eqtl_scan,
    eqtl_effect_sizes,
    eqtl_summary_stats,
)
from metainformant.gwas.visualization.eqtl_visualization import (
    plot_eqtl_boxplot,
    plot_eqtl_summary,
    plot_eqtl_volcano,
)

# Output directories
OUTPUT_DIR = Path("output/eqtl/amellifera")
RESULTS_DIR = OUTPUT_DIR / "results"
PLOTS_DIR = OUTPUT_DIR / "plots"
LOGS_DIR = OUTPUT_DIR / "logs"

for directory in (RESULTS_DIR, PLOTS_DIR, LOGS_DIR):
    directory.mkdir(parents=True, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("output/eqtl/amellifera/logs/eqtl_analysis.log"),
    ],
)
logger = logging.getLogger(__name__)


def run_eqtl_analysis():
    """Run the full eQTL analysis pipeline."""
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("Starting eQTL Analysis Pipeline")
    logger.info(f"Start time: {start_time}")
    logger.info("=" * 60)

    # Step 1: Create/load data
    logger.info("\n[Step 1] Preparing data...")
    expr, geno, gene_pos, var_pos = create_synthetic_data(n_genes=100, n_variants=500, n_samples=50)

    # Save input data
    expr.to_csv(RESULTS_DIR / "expression_matrix.tsv", sep="\t")
    geno.to_csv(RESULTS_DIR / "genotype_matrix.tsv", sep="\t")
    gene_pos.to_csv(RESULTS_DIR / "gene_positions.tsv", sep="\t", index=False)
    var_pos.to_csv(RESULTS_DIR / "variant_positions.tsv", sep="\t", index=False)
    logger.info(f"Saved input data to {RESULTS_DIR}")

    # Step 2: Run cis-eQTL scan
    logger.info("\n[Step 2] Running cis-eQTL scan...")
    cis_results = cis_eqtl_scan(
        expression_matrix=expr,
        genotype_matrix=geno,
        gene_positions=gene_pos,
        variant_positions=var_pos,
        cis_window=1_000_000,
        maf_threshold=0.05,
    )

    cis_results.to_csv(RESULTS_DIR / "cis_eqtl_results.tsv", sep="\t", index=False)
    logger.info(f"cis-eQTL results: {len(cis_results)} tests")

    # Step 3: Compute effect sizes
    logger.info("\n[Step 3] Computing effect sizes...")
    if len(cis_results) > 0:
        top_hits = cis_results.nsmallest(50, "pvalue")
        effect_sizes = eqtl_effect_sizes(expr, geno, top_hits)
        effect_sizes.to_csv(RESULTS_DIR / "effect_sizes.tsv", sep="\t", index=False)
        logger.info(f"Computed effect sizes for {len(effect_sizes)} top hits")
    else:
        effect_sizes = pd.DataFrame()

    # Step 4: Generate summary statistics
    logger.info("\n[Step 4] Generating summary statistics...")
    summary = eqtl_summary_stats(cis_results, fdr_threshold=0.05)

    with open(RESULTS_DIR / "summary_stats.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Summary: {summary}")

    # Step 5: Create visualizations
    logger.info("\n[Step 5] Creating visualizations...")

    if len(cis_results) > 0:
        # Volcano plot
        try:
            plot_eqtl_volcano(
                cis_results,
                fdr_threshold=0.05,
                output_path=PLOTS_DIR / "volcano_plot.png",
            )
            logger.info("Created volcano plot")
        except Exception as e:
            logger.warning(f"Could not create volcano plot: {e}")

        # Summary plot
        try:
            plot_eqtl_summary(
                summary,
                output_path=PLOTS_DIR / "summary_plot.png",
            )
            logger.info("Created summary plot")
        except Exception as e:
            logger.warning(f"Could not create summary plot: {e}")

        # Boxplot for top hit
        if len(cis_results) > 0:
            try:
                top_hit = cis_results.nsmallest(1, "pvalue").iloc[0]
                gene_id = top_hit["gene_id"]
                var_id = top_hit["variant_id"]

                plot_eqtl_boxplot(
                    expression=expr.loc[gene_id].values,
                    genotypes=geno.loc[var_id].values,
                    gene_id=gene_id,
                    variant_id=var_id,
                    output_path=PLOTS_DIR / "top_eqtl_boxplot.png",
                )
                logger.info(f"Created boxplot for top eQTL: {gene_id} ~ {var_id}")
            except Exception as e:
                logger.warning(f"Could not create boxplot: {e}")

    # Step 6: Final summary
    end_time = datetime.now()
    elapsed = (end_time - start_time).total_seconds()

    run_summary = {
        "start_time": start_time.isoformat(),
        "end_time": end_time.isoformat(),
        "elapsed_seconds": elapsed,
        "n_genes": len(gene_pos),
        "n_variants": len(var_pos),
        "n_samples": len(expr.columns),
        "n_tests": len(cis_results),
        "n_significant": summary.get("n_eqtls", 0),
        "n_egenes": summary.get("n_egenes", 0),
    }

    with open(RESULTS_DIR / "run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    logger.info("\n" + "=" * 60)
    logger.info("eQTL Analysis Complete!")
    logger.info(f"Elapsed time: {elapsed:.2f} seconds")
    logger.info("=" * 60)
    logger.info("\nOutput locations:")
    logger.info(f"  Results: {RESULTS_DIR}")
    logger.info(f"  Plots:   {PLOTS_DIR}")
    logger.info(f"  Logs:    {LOGS_DIR}")

    return run_summary


if __name__ == "__main__":
    run_eqtl_analysis()
