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

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

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

# Output directories
OUTPUT_DIR = Path("output/eqtl/amellifera")
RESULTS_DIR = OUTPUT_DIR / "results"
PLOTS_DIR = OUTPUT_DIR / "plots"
LOGS_DIR = OUTPUT_DIR / "logs"


def create_synthetic_data(
    n_genes: int = 100, n_variants: int = 500, n_samples: int = 50
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create synthetic expression and genotype data for demo.
    
    Simulates realistic eQTL patterns including:
    - Cis-eQTLs with strong effects
    - Background noise variants
    """
    np.random.seed(42)
    
    logger.info(f"Generating synthetic data: {n_genes} genes, {n_variants} variants, {n_samples} samples")
    
    # Sample IDs
    sample_ids = [f"SRR{10000000 + i}" for i in range(n_samples)]
    
    # Gene positions (spread across chromosomes 1-5)
    gene_ids = [f"LOC{100000 + i}" for i in range(n_genes)]
    chromosomes = [str((i % 5) + 1) for i in range(n_genes)]
    tss_positions = [1_000_000 + (i * 50_000) for i in range(n_genes)]
    
    gene_positions = pd.DataFrame({
        "gene_id": gene_ids,
        "chrom": chromosomes,
        "tss_position": tss_positions,
    })
    
    # Variant positions (5 per gene in cis window)
    variant_ids = []
    var_chromosomes = []
    var_positions = []
    
    for gene_idx in range(n_genes):
        for var_offset in range(5):
            variant_ids.append(f"rs{gene_idx * 5 + var_offset}")
            var_chromosomes.append(chromosomes[gene_idx])
            # Position within cis window
            var_positions.append(tss_positions[gene_idx] + (var_offset - 2) * 10_000)
    
    variant_positions = pd.DataFrame({
        "variant_id": variant_ids,
        "chrom": var_chromosomes,
        "position": var_positions,
    })
    
    # Genotypes: dosages 0, 1, 2 with MAF ~ 0.3
    genotypes = np.random.choice([0, 1, 2], size=(len(variant_ids), n_samples), p=[0.5, 0.35, 0.15])
    genotype_matrix = pd.DataFrame(genotypes, index=variant_ids, columns=sample_ids)
    
    # Expression: baseline + eQTL effect for ~30% of gene-variant pairs
    expression = np.zeros((n_genes, n_samples))
    
    for gene_idx in range(n_genes):
        # Baseline expression (log-scale, ~5-10)
        baseline = np.random.uniform(5, 10)
        noise = np.random.normal(0, 0.5, n_samples)
        
        # Add eQTL effect for first variant of each gene (30% chance)
        if np.random.random() < 0.3:
            var_idx = gene_idx * 5  # First variant for this gene
            eqtl_effect = np.random.uniform(0.5, 2.0)
            effect = genotypes[var_idx] * eqtl_effect
        else:
            effect = 0
        
        expression[gene_idx] = baseline + effect + noise
    
    expression_matrix = pd.DataFrame(expression, index=gene_ids, columns=sample_ids)
    
    logger.info(f"Created expression matrix: {expression_matrix.shape}")
    logger.info(f"Created genotype matrix: {genotype_matrix.shape}")
    
    return expression_matrix, genotype_matrix, gene_positions, variant_positions


def run_eqtl_analysis():
    """Run the full eQTL analysis pipeline."""
    start_time = datetime.now()
    logger.info("=" * 60)
    logger.info("Starting eQTL Analysis Pipeline")
    logger.info(f"Start time: {start_time}")
    logger.info("=" * 60)
    
    # Step 1: Create/load data
    logger.info("\n[Step 1] Preparing data...")
    expr, geno, gene_pos, var_pos = create_synthetic_data(
        n_genes=100, n_variants=500, n_samples=50
    )
    
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
