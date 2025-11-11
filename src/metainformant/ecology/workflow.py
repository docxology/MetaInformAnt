"""Ecology analysis workflow orchestration."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from metainformant.core import io, paths
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def run_community_analysis_workflow(
    abundance_file: Path,
    output_dir: Path,
    calculate_diversity: bool = True,
    calculate_beta: bool = True,
    environmental_file: Path | None = None,
) -> Dict[str, Any]:
    """Run complete community ecology analysis workflow.
    
    Args:
        abundance_file: Path to species abundance table (sites x species)
        output_dir: Output directory for results
        calculate_diversity: If True, calculate alpha diversity metrics
        calculate_beta: If True, calculate beta diversity metrics
        environmental_file: Optional environmental variables file
        
    Returns:
        Dictionary with workflow results
    """
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting community analysis workflow: {abundance_file}")
    
    # Load abundance data
    abundance_df = pd.read_csv(abundance_file, index_col=0)
    logger.info(f"Loaded abundance data: {abundance_df.shape[0]} sites, {abundance_df.shape[1]} species")
    
    results = {
        "input_file": str(abundance_file),
        "n_sites": abundance_df.shape[0],
        "n_species": abundance_df.shape[1],
    }
    
    # Calculate diversity metrics
    if calculate_diversity:
        from metainformant.ecology import (
            chao1_estimator,
            pielou_evenness,
            shannon_diversity,
            simpson_diversity,
            species_richness,
        )
        
        diversity_results = []
        
        for site_id, abundances in abundance_df.iterrows():
            abund_list = abundances.values.tolist()
            
            diversity_results.append({
                "site_id": site_id,
                "species_richness": species_richness(abund_list),
                "shannon_diversity": shannon_diversity(abund_list),
                "simpson_diversity": simpson_diversity(abund_list),
                "pielou_evenness": pielou_evenness(abund_list),
                "chao1_estimate": chao1_estimator([int(x) for x in abund_list]),
            })
        
        diversity_df = pd.DataFrame(diversity_results)
        diversity_output = output_dir / "diversity_metrics.tsv"
        diversity_df.to_csv(diversity_output, sep="\t", index=False)
        results["diversity_file"] = str(diversity_output)
        logger.info(f"Calculated diversity metrics, saved to {diversity_output}")
    
    # Calculate beta diversity
    if calculate_beta and abundance_df.shape[0] > 1:
        from metainformant.ecology import beta_diversity_partitioning, bray_curtis_dissimilarity
        
        site_abundances = [abundances.values.tolist() for _, abundances in abundance_df.iterrows()]
        
        beta_results = beta_diversity_partitioning(site_abundances)
        beta_output = output_dir / "beta_diversity.json"
        io.dump_json(beta_results, beta_output)
        results["beta_diversity_file"] = str(beta_output)
        logger.info(f"Calculated beta diversity, saved to {beta_output}")
    
    # Environmental gradient analysis
    if environmental_file:
        from metainformant.ecology import analyze_environmental_gradient
        
        env_df = pd.read_csv(environmental_file, index_col=0)
        ordination = analyze_environmental_gradient(abundance_df, env_df)
        
        ordination_output = output_dir / "ordination_results.json"
        # Convert DataFrames to dict for JSON serialization
        ordination_dict = {
            "scores": ordination["scores"].to_dict(),
            "species_scores": ordination["species_scores"].to_dict(),
            "environmental_loadings": ordination["environmental_loadings"].to_dict(),
            "explained_variance": ordination["explained_variance"],
        }
        io.dump_json(ordination_dict, ordination_output)
        results["ordination_file"] = str(ordination_output)
        logger.info(f"Performed ordination analysis, saved to {ordination_output}")
    
    results["status"] = "completed"
    return results

