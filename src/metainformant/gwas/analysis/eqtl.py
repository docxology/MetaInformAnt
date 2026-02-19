"""
eQTL Analysis Wrapper.

This module provides high-level interfaces for running expression Quantitative
Trait Locus (eQTL) analysis, linking genotype data with RNA-seq expression levels.
"""

import pandas as pd
from typing import Dict, List, Optional, Tuple, Union

from metainformant.core.utils.logging import get_logger
from metainformant.gwas.analysis.association import association_test_linear

logger = get_logger(__name__)


def run_eqtl_analysis(
    genotype_matrix: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    covariates: Optional[pd.DataFrame] = None,
    cis_window: int = 1000000
) -> pd.DataFrame:
    """
    Run basic linear eQTL analysis for all transcripts.
    
    Note: A full all-vs-all eQTL scan is computationally expensive.
    This is a simplified wrapper for verification of the data pipeline.
    
    Args:
        genotype_matrix: (Samples x Variants) - 0/1/2 dosage.
        expression_matrix: (Samples x Transcripts) - Normalized expression.
        covariates: (Samples x Covariates) - Optional.
        cis_window: Base pairs matching window (not used in simple linear test, placeholder).
        
    Returns:
        pd.DataFrame: Summary statistics (variant, transcript, beta, p_value).
    """
    results = []
    
    # Align samples
    common_samples = genotype_matrix.index.intersection(expression_matrix.index)
    
    if covariates is not None:
        common_samples = common_samples.intersection(covariates.index)
        
    if len(common_samples) == 0:
        logger.error("No common samples between genotype and expression matrices.")
        return pd.DataFrame()
        
    logger.info(f"Running eQTL analysis on {len(common_samples)} samples.")
    
    # Subset data
    G = genotype_matrix.loc[common_samples]
    E = expression_matrix.loc[common_samples]
    C = covariates.loc[common_samples] if covariates is not None else None
    
    # For this verification implementation, we run a naive loop
    # In production, use TensorQTL or MatrixeQTL approaches
    
    # Limit to first few transcripts for speed/demo if matrix is huge
    max_transcripts = 10 
    transcripts = E.columns[:max_transcripts]
    
    for transcript_id in transcripts:
        phenotype = E[transcript_id]
        
        # Run association test (Variant vs Phenotype)
        # This returns specific stats for each variant against THIS transcript
        try:
            assoc_stats = association_test_linear(
                genotypes=G,
                phenotypes=phenotype,
                covariates=C
            )
            
            # Add transcript ID
            assoc_stats["transcript"] = transcript_id
            results.append(assoc_stats)
            
        except Exception as e:
            logger.warning(f"Failed eQTL for {transcript_id}: {e}")
            
    if not results:
        return pd.DataFrame()
        
    final_df = pd.concat(results, ignore_index=True)
    return final_df
