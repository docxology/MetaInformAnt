"""Comprehensive visualization generation for GWAS results.

This module provides a unified interface to generate all available
visualization types from GWAS analysis results.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)


def generate_all_plots(
    association_results: Path,
    output_dir: Path,
    *,
    pca_file: Path | None = None,
    kinship_file: Path | None = None,
    vcf_file: Path | None = None,
    significance_threshold: float = 5e-8,
) -> dict[str, Any]:
    """Generate all available GWAS visualization plots.
    
    Creates a comprehensive set of plots from GWAS results including:
    - Genome-wide views (Manhattan, circular, ideogram)
    - Statistical diagnostics (QQ, lambda GC, volcano)
    - Population structure (PCA, kinship)
    - Variant properties (MAF, density, Ts/Tv)
    - Effect sizes (forest plot, direction)
    
    Args:
        association_results: Path to association results TSV
        output_dir: Directory for output plots
        pca_file: Path to PCA components file
        kinship_file: Path to kinship matrix file
        vcf_file: Path to VCF file (for variant-level plots)
        significance_threshold: P-value threshold for significance
    
    Returns:
        Dictionary with plot generation results and paths
    """
    logger.info(f"generate_all_plots: Creating comprehensive visualization suite")
    
    ensure_directory(output_dir)
    plots_generated = {}
    errors = []
    
    # Genome-wide visualizations
    logger.info("generate_all_plots: Creating genome-wide visualizations")
    try:
        from .visualization_genome import (
            manhattan_plot,
            circular_manhattan_plot,
            chromosome_ideogram,
        )
        
        result = manhattan_plot(
            association_results,
            output_dir / "manhattan_genome.png",
            significance_threshold=significance_threshold,
        )
        plots_generated["manhattan_genome"] = result
        
        result = circular_manhattan_plot(
            association_results,
            output_dir / "circular_manhattan.png",
            significance_threshold=significance_threshold,
        )
        plots_generated["circular_manhattan"] = result
        
        result = chromosome_ideogram(
            association_results,
            output_dir / "chromosome_ideogram.png",
            significance_threshold=significance_threshold,
        )
        plots_generated["chromosome_ideogram"] = result
        
    except Exception as e:
        logger.error(f"generate_all_plots: Genome-wide visualization error: {e}")
        errors.append({"category": "genome", "error": str(e)})
    
    # Statistical diagnostics
    logger.info("generate_all_plots: Creating statistical diagnostic plots")
    try:
        from .visualization_statistical import (
            qq_plot,
            qq_plot_stratified,
            lambda_gc_plot,
            volcano_plot,
        )
        
        result = qq_plot(
            association_results,
            output_dir / "qq_plot.png",
            show_ci=True,
            show_lambda_gc=True,
        )
        plots_generated["qq_plot"] = result
        
        result = qq_plot_stratified(
            association_results,
            output_dir / "qq_plot_stratified.png",
        )
        plots_generated["qq_plot_stratified"] = result
        
        result = lambda_gc_plot(
            association_results,
            output_dir / "lambda_gc_by_chrom.png",
        )
        plots_generated["lambda_gc"] = result
        
        result = volcano_plot(
            association_results,
            output_dir / "volcano_plot.png",
            significance_threshold=significance_threshold,
        )
        plots_generated["volcano_plot"] = result
        
    except Exception as e:
        logger.error(f"generate_all_plots: Statistical visualization error: {e}")
        errors.append({"category": "statistical", "error": str(e)})
    
    # Population structure visualizations
    if pca_file and pca_file.exists():
        logger.info("generate_all_plots: Creating population structure plots")
        try:
            from .visualization_population import (
                pca_plot,
                pca_scree_plot,
                kinship_heatmap,
                population_tree,
            )
            
            result = pca_plot(
                pca_file,
                output_dir / "pca_plot.png",
                pc1=1,
                pc2=2,
            )
            plots_generated["pca_plot"] = result
            
            result = pca_scree_plot(
                pca_file,
                output_dir / "pca_scree.png",
                max_pcs=20,
            )
            plots_generated["pca_scree"] = result
            
            if kinship_file and kinship_file.exists():
                result = kinship_heatmap(
                    kinship_file,
                    output_dir / "kinship_heatmap.png",
                    max_samples=200,
                )
                plots_generated["kinship_heatmap"] = result
                
                result = population_tree(
                    kinship_file,
                    output_dir / "population_tree.png",
                )
                plots_generated["population_tree"] = result
            
        except Exception as e:
            logger.error(f"generate_all_plots: Population structure error: {e}")
            errors.append({"category": "population", "error": str(e)})
    
    # Variant property visualizations
    logger.info("generate_all_plots: Creating variant property plots")
    try:
        from .visualization_variants import (
            maf_distribution,
            variant_density_plot,
            transition_transversion_plot,
        )
        
        result = maf_distribution(
            association_results,
            output_dir / "maf_distribution.png",
            bins=50,
        )
        plots_generated["maf_distribution"] = result
        
        result = variant_density_plot(
            association_results,
            output_dir / "variant_density.png",
            window_size=1_000_000,
        )
        plots_generated["variant_density"] = result
        
        result = transition_transversion_plot(
            association_results,
            output_dir / "ts_tv_ratio.png",
        )
        plots_generated["ts_tv_ratio"] = result
        
    except Exception as e:
        logger.error(f"generate_all_plots: Variant property error: {e}")
        errors.append({"category": "variants", "error": str(e)})
    
    # Effect size visualizations
    logger.info("generate_all_plots: Creating effect size plots")
    try:
        from .visualization_effects import (
            effect_size_forest_plot,
            effect_direction_plot,
        )
        
        result = effect_size_forest_plot(
            association_results,
            output_dir / "forest_plot.png",
            top_n=20,
            significance_threshold=significance_threshold,
        )
        plots_generated["forest_plot"] = result
        
        result = effect_direction_plot(
            association_results,
            output_dir / "effect_direction.png",
            significance_threshold=significance_threshold,
        )
        plots_generated["effect_direction"] = result
        
    except Exception as e:
        logger.error(f"generate_all_plots: Effect size error: {e}")
        errors.append({"category": "effects", "error": str(e)})
    
    # Summary
    successful_plots = [k for k, v in plots_generated.items() 
                       if v.get("status") == "success"]
    skipped_plots = [k for k, v in plots_generated.items() 
                    if v.get("status") == "skipped"]
    
    logger.info(f"generate_all_plots: Generated {len(successful_plots)} plots, "
               f"{len(skipped_plots)} skipped, {len(errors)} errors")
    
    return {
        "status": "completed",
        "num_plots_generated": len(successful_plots),
        "num_plots_skipped": len(skipped_plots),
        "num_errors": len(errors),
        "plots": plots_generated,
        "errors": errors,
        "successful_plots": successful_plots,
        "skipped_plots": skipped_plots,
    }





