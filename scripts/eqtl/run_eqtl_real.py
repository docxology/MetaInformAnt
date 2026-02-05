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
        logging.FileHandler("output/eqtl/amellifera/logs/eqtl_real_analysis.log"),
    ],
)
logger = logging.getLogger(__name__)

# Paths
QUANT_DIR = Path("output/amalgkit/apis_mellifera_all/work/quant")
OUTPUT_DIR = Path("output/eqtl/amellifera")
RESULTS_DIR = OUTPUT_DIR / "results"
PLOTS_DIR = OUTPUT_DIR / "plots"
LOGS_DIR = OUTPUT_DIR / "logs"


def load_real_expression_data(
    quant_dir: Path, max_samples: int = 100, min_tpm: float = 1.0, max_genes: int = 1000
) -> tuple[pd.DataFrame, list[str]]:
    """Load real kallisto expression data from abundance.tsv files.
    
    Args:
        quant_dir: Directory containing sample subdirectories.
        max_samples: Maximum samples to load (for speed).
        min_tpm: Minimum mean TPM to include a gene.
        max_genes: Maximum genes to include (top by mean expression).
    
    Returns:
        Expression matrix (genes x samples) and list of sample IDs.
    """
    logger.info(f"Loading real expression data from {quant_dir}")
    
    sample_dirs = sorted([d for d in quant_dir.iterdir() if d.is_dir()])
    logger.info(f"Found {len(sample_dirs)} sample directories")
    
    # Limit samples for speed
    sample_dirs = sample_dirs[:max_samples]
    logger.info(f"Loading {len(sample_dirs)} samples")
    
    expression_data = {}
    gene_info = None
    
    for sample_dir in sample_dirs:
        abundance_file = sample_dir / "abundance.tsv"
        if not abundance_file.exists():
            continue
        
        sample_id = sample_dir.name
        df = pd.read_csv(abundance_file, sep="\t")
        
        if gene_info is None:
            gene_info = df[["target_id", "length"]].copy()
        
        expression_data[sample_id] = df.set_index("target_id")["tpm"]
    
    if not expression_data:
        raise ValueError("No expression data loaded")
    
    # Build expression matrix
    expr_matrix = pd.DataFrame(expression_data)
    logger.info(f"Loaded expression matrix: {expr_matrix.shape}")
    
    # Filter low-expression genes
    mean_tpm = expr_matrix.mean(axis=1)
    expressed_genes = mean_tpm[mean_tpm >= min_tpm].index
    expr_matrix = expr_matrix.loc[expressed_genes]
    logger.info(f"After filtering (TPM >= {min_tpm}): {len(expr_matrix)} genes")
    
    # Keep only top genes by expression for speed
    if len(expr_matrix) > max_genes:
        top_genes = mean_tpm.loc[expr_matrix.index].nlargest(max_genes).index
        expr_matrix = expr_matrix.loc[top_genes]
        logger.info(f"Selected top {max_genes} genes by expression")
    
    return expr_matrix, list(expression_data.keys())


def parse_gene_positions(gene_ids: list[str]) -> pd.DataFrame:
    """Parse gene positions from kallisto target IDs.
    
    Target format: lcl|NC_037638.1_mrna_XM_623972.6_1
    """
    positions = []
    
    for gid in gene_ids:
        parts = gid.split("_")
        if len(parts) >= 3:
            # Extract chromosome from NC_XXXXXX.1
            chrom_part = parts[0].replace("lcl|", "")
            if chrom_part.startswith("NC_"):
                chrom = chrom_part
            else:
                chrom = "unknown"
            
            # Use index as position (real positions would come from GFF)
            # Spread across genome (~250Mb / n_genes)
            idx = gene_ids.index(gid)
            position = 1_000_000 + (idx * 10_000)
            
            positions.append({
                "gene_id": gid,
                "chrom": chrom,
                "tss_position": position,
            })
    
    return pd.DataFrame(positions)


def create_synthetic_genotypes(
    sample_ids: list[str], gene_positions: pd.DataFrame, variants_per_gene: int = 1
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create synthetic genotype data for demonstration.
    
    Note: In real analysis, this would come from matched WGS/genotyping data.
    """
    logger.info("Creating synthetic genotypes (real genotypes not available)")
    np.random.seed(42)
    
    variant_ids = []
    var_chroms = []
    var_positions = []
    
    for _, row in gene_positions.iterrows():
        for i in range(variants_per_gene):
            var_id = f"var_{row['gene_id'][:20]}_{i}"
            variant_ids.append(var_id)
            var_chroms.append(row["chrom"])
            var_positions.append(int(row["tss_position"]) + (i - 1) * 5000)
    
    variant_positions = pd.DataFrame({
        "variant_id": variant_ids,
        "chrom": var_chroms,
        "position": var_positions,
    })
    
    # Generate dosages with MAF ~ 0.25
    n_variants = len(variant_ids)
    n_samples = len(sample_ids)
    genotypes = np.random.choice([0, 1, 2], size=(n_variants, n_samples), p=[0.56, 0.32, 0.12])
    
    genotype_matrix = pd.DataFrame(genotypes, index=variant_ids, columns=sample_ids)
    
    return genotype_matrix, variant_positions


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
    gene_stats = pd.DataFrame({
        "gene_id": expr_matrix.index,
        "mean_tpm": expr_matrix.mean(axis=1),
        "std_tpm": expr_matrix.std(axis=1),
        "n_samples_expressed": (expr_matrix > 0).sum(axis=1),
    })
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
    gene_positions.to_csv(RESULTS_DIR / "gene_positions_real.tsv", sep="\t", index=False)
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
            gene_stats[["gene_id", "mean_tpm"]],
            on="gene_id",
            how="left"
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
                logger.info(f"Created boxplot for top eQTL")
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
