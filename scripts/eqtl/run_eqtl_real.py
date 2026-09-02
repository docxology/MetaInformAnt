#!/usr/bin/env python
"""Real A. mellifera eQTL Analysis Pipeline.

Uses actual quantified RNA-seq expression data from 2,300+ samples.
Since we don't have matched genotype data, we generate synthetic
genotypes but demonstrate the full analytical workflow.

Usage:
    uv run python scripts/eqtl/run_eqtl_real.py
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

from metainformant.eqtl.synthetic import (
    create_synthetic_genotypes,
    load_real_expression_data,
    parse_gene_positions,
)
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

# Paths
QUANT_DIR = Path("output/amalgkit/apis_mellifera_all/work/quant")
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
        logging.FileHandler("output/eqtl/amellifera/logs/eqtl_real_analysis.log"),
    ],
)
logger = logging.getLogger(__name__)


def run_real_eqtl_analysis():
    """Run eQTL analysis with real A. mellifera expression data."""
    start_time = datetime.now()
    logger.info("=" * 70)
    logger.info("Real A. mellifera eQTL Analysis Pipeline")
    logger.info(f"Start time: {start_time}")
    logger.info("=" * 70)

    # Step 1: Load real expression data
    logger.info("\n[Step 1] Loading REAL expression data...")
    expr_matrix, sample_ids = load_real_expression_data(
        QUANT_DIR, max_samples=200, min_tpm=1.0
    )

    # Step 2: Parse gene positions
    logger.info("\n[Step 2] Parsing gene annotations...")
    gene_positions = parse_gene_positions(list(expr_matrix.index))
    logger.info(f"Parsed positions for {len(gene_positions)} genes")

    # Save gene info with expression stats
    gene_stats = pd.DataFrame(
        {
            "gene_id": expr_matrix.index,
            "mean_tpm": expr_matrix.mean(axis=1),
            "std_tpm": expr_matrix.std(axis=1),
            "n_samples_expressed": (expr_matrix > 0).sum(axis=1),
        }
    )
    gene_stats.to_csv(RESULTS_DIR / "gene_expression_stats.tsv", sep="\t", index=False)

    # Step 3: Create synthetic genotypes (real WGS not available)
    logger.info("\n[Step 3] Creating synthetic genotypes...")
    geno_matrix, var_positions = create_synthetic_genotypes(
        sample_ids, gene_positions, variants_per_gene=3
    )
    logger.info(f"Created {len(var_positions)} variants")

    # Save inputs
    expr_matrix.to_csv(RESULTS_DIR / "real_expression_matrix.tsv", sep="\t")
    geno_matrix.to_csv(RESULTS_DIR / "synthetic_genotype_matrix.tsv", sep="\t")
    gene_positions.to_csv(
        RESULTS_DIR / "gene_positions_real.tsv", sep="\t", index=False
    )
    var_positions.to_csv(RESULTS_DIR / "variant_positions.tsv", sep="\t", index=False)

    # Step 4: Run cis-eQTL scan
    logger.info("\n[Step 4] Running cis-eQTL scan...")
    cis_results = cis_eqtl_scan(
        expression_matrix=expr_matrix,
        genotype_matrix=geno_matrix,
        gene_positions=gene_positions,
        variant_positions=var_positions,
        cis_window=500_000,  # 500kb window
        maf_threshold=0.05,
    )

    cis_results.to_csv(RESULTS_DIR / "cis_eqtl_results_real.tsv", sep="\t", index=False)
    logger.info(f"cis-eQTL results: {len(cis_results)} tests performed")

    # Step 5: Annotate results with gene names
    logger.info("\n[Step 5] Annotating results...")
    if len(cis_results) > 0:
        # Add gene expression stats
        cis_results_annotated = cis_results.merge(
            gene_stats[["gene_id", "mean_tpm"]], on="gene_id", how="left"
        )
        cis_results_annotated.to_csv(
            RESULTS_DIR / "cis_eqtl_annotated.tsv", sep="\t", index=False
        )

        # Top hits
        top_hits = cis_results.nsmallest(100, "pvalue")
        top_hits.to_csv(RESULTS_DIR / "top_100_eqtls.tsv", sep="\t", index=False)
        effect_sizes = eqtl_effect_sizes(expr_matrix, geno_matrix, top_hits)
        effect_sizes.to_csv(RESULTS_DIR / "top_effect_sizes.tsv", sep="\t", index=False)

    # Step 6: Summary statistics
    logger.info("\n[Step 6] Computing summary statistics...")
    summary = eqtl_summary_stats(cis_results, fdr_threshold=0.05)
    summary["species"] = "Apis mellifera"
    summary["n_samples"] = len(sample_ids)
    summary["n_genes_tested"] = len(gene_positions)
    summary["data_source"] = "Real RNA-seq from ENA/NCBI"

    with open(RESULTS_DIR / "summary_stats_real.json", "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Summary: {json.dumps(summary, indent=2)}")

    # Step 7: Visualizations
    logger.info("\n[Step 7] Creating visualizations...")

    if len(cis_results) > 0:
        # Volcano plot
        try:
            plot_eqtl_volcano(
                cis_results,
                fdr_threshold=0.05,
                output_path=PLOTS_DIR / "volcano_plot_real.png",
                title="A. mellifera eQTL Volcano Plot (Real Expression Data)",
            )
            logger.info("Created volcano plot")
        except Exception as e:
            logger.warning(f"Could not create volcano plot: {e}")

        # Summary plot
        try:
            plot_eqtl_summary(summary, output_path=PLOTS_DIR / "summary_plot_real.png")
            logger.info("Created summary plot")
        except Exception as e:
            logger.warning(f"Could not create summary plot: {e}")

        # Boxplot for top eQTL
        if len(cis_results) > 0:
            try:
                top_hit = cis_results.nsmallest(1, "pvalue").iloc[0]
                gene_id = top_hit["gene_id"]
                var_id = top_hit["variant_id"]

                plot_eqtl_boxplot(
                    expression=expr_matrix.loc[gene_id].values,
                    genotypes=geno_matrix.loc[var_id].values,
                    gene_id=gene_id[:30],  # Truncate for display
                    variant_id=var_id[:30],
                    output_path=PLOTS_DIR / "top_eqtl_boxplot_real.png",
                )
                logger.info("Created boxplot for top eQTL")
            except Exception as e:
                logger.warning(f"Could not create boxplot: {e}")

    # Final summary
    end_time = datetime.now()
    elapsed = (end_time - start_time).total_seconds()

    run_summary = {
        "analysis_type": "Real A. mellifera eQTL Analysis",
        "start_time": start_time.isoformat(),
        "end_time": end_time.isoformat(),
        "elapsed_seconds": elapsed,
        "expression_source": str(QUANT_DIR),
        "n_samples": len(sample_ids),
        "n_genes": len(gene_positions),
        "n_variants": len(var_positions),
        "n_tests": len(cis_results),
        "n_significant": summary.get("n_eqtls", 0),
        "n_egenes": summary.get("n_egenes", 0),
        "note": "Expression data is REAL from A. mellifera RNA-seq; genotypes are synthetic",
    }

    with open(RESULTS_DIR / "run_summary_real.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    logger.info("\n" + "=" * 70)
    logger.info("Real A. mellifera eQTL Analysis Complete!")
    logger.info(f"Elapsed time: {elapsed:.2f} seconds")
    logger.info("=" * 70)
    logger.info("\nOutput locations:")
    logger.info(f"  Results: {RESULTS_DIR.absolute()}")
    logger.info(f"  Plots:   {PLOTS_DIR.absolute()}")
    logger.info(f"  Logs:    {LOGS_DIR.absolute()}")
    logger.info("\nKey files:")
    logger.info("  - real_expression_matrix.tsv  (REAL data)")
    logger.info("  - cis_eqtl_annotated.tsv      (annotated results)")
    logger.info("  - top_100_eqtls.tsv           (significant hits)")
    logger.info("  - volcano_plot_real.png       (visualization)")

    return run_summary


if __name__ == "__main__":
    run_real_eqtl_analysis()
