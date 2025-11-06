"""RNA-Protein integration utilities for translation efficiency and abundance prediction."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def calculate_translation_efficiency(
    rna_expression: pd.DataFrame,
    protein_abundance: pd.DataFrame,
    method: str = "correlation",
) -> pd.DataFrame:
    """Calculate translation efficiency from RNA and protein data.
    
    Translation efficiency measures how efficiently mRNA is translated into protein.
    Common metrics include protein/mRNA ratio or correlation-based measures.
    
    Args:
        rna_expression: DataFrame with samples as rows and genes as columns (mRNA levels)
        protein_abundance: DataFrame with samples as rows and proteins as columns
        method: Method for calculating efficiency ('ratio', 'correlation', 'slope')
        
    Returns:
        DataFrame with translation efficiency metrics:
        - 'gene_id': Gene/protein identifier
        - 'efficiency': Translation efficiency value
        - 'method': Method used
    """
    if rna_expression.empty or protein_abundance.empty:
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])
    
    # Align samples
    common_samples = set(rna_expression.index) & set(protein_abundance.index)
    if not common_samples:
        logger.warning("No common samples between RNA and protein data")
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])
    
    rna_aligned = rna_expression.loc[list(common_samples)]
    protein_aligned = protein_abundance.loc[list(common_samples)]
    
    # Find common genes/proteins
    common_genes = set(rna_aligned.columns) & set(protein_aligned.columns)
    
    if not common_genes:
        logger.warning("No common genes between RNA and protein data")
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])
    
    efficiencies = []
    
    for gene_id in common_genes:
        rna_values = rna_aligned[gene_id].values
        protein_values = protein_aligned[gene_id].values
        
        # Filter out zeros and missing values
        valid_mask = (rna_values > 0) & (protein_values > 0) & np.isfinite(rna_values) & np.isfinite(protein_values)
        
        if np.sum(valid_mask) < 3:  # Need at least 3 samples
            continue
        
        rna_valid = rna_values[valid_mask]
        protein_valid = protein_values[valid_mask]
        
        if method == "ratio":
            # Mean protein/mRNA ratio
            efficiency = np.mean(protein_valid / rna_valid)
        elif method == "correlation":
            # Correlation between RNA and protein
            try:
                efficiency = np.corrcoef(rna_valid, protein_valid)[0, 1]
                if np.isnan(efficiency):
                    continue
            except Exception:
                continue
        elif method == "slope":
            # Linear regression slope
            try:
                slope = np.polyfit(rna_valid, protein_valid, 1)[0]
                efficiency = slope
            except Exception:
                continue
        else:
            logger.warning(f"Unknown method: {method}, using ratio")
            efficiency = np.mean(protein_valid / rna_valid)
        
        efficiencies.append({
            "gene_id": gene_id,
            "efficiency": float(efficiency),
            "method": method,
        })
    
    if not efficiencies:
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])
    
    return pd.DataFrame(efficiencies)


def predict_protein_abundance_from_rna(
    rna_expression: pd.DataFrame,
    training_rna: pd.DataFrame | None = None,
    training_protein: pd.DataFrame | None = None,
    method: str = "linear",
) -> pd.DataFrame:
    """Predict protein abundance from RNA expression levels.
    
    Uses a trained model (if provided) or simple linear relationship to predict
    protein levels from mRNA expression.
    
    Args:
        rna_expression: RNA expression data to predict from (samples x genes)
        training_rna: Optional training RNA data for model fitting
        training_protein: Optional training protein data for model fitting
        method: Prediction method ('linear', 'log_linear', 'power_law')
        
    Returns:
        DataFrame with predicted protein abundance (samples x proteins)
    """
    if rna_expression.empty:
        return pd.DataFrame()
    
    # If training data provided, fit model
    if training_rna is not None and training_protein is not None:
        # Align training data
        common_samples = set(training_rna.index) & set(training_protein.index)
        common_genes = set(training_rna.columns) & set(training_protein.columns)
        
        if not common_samples or not common_genes:
            logger.warning("Training data alignment failed, using simple prediction")
            return _simple_protein_prediction(rna_expression, method)
        
        # Fit models per gene
        predictions = {}
        
        for gene_id in common_genes:
            rna_train = training_rna.loc[list(common_samples), gene_id].values
            protein_train = training_protein.loc[list(common_samples), gene_id].values
            
            # Filter valid values
            valid_mask = (rna_train > 0) & (protein_train > 0) & np.isfinite(rna_train) & np.isfinite(protein_train)
            
            if np.sum(valid_mask) < 3:
                continue
            
            rna_valid = rna_train[valid_mask]
            protein_valid = protein_train[valid_mask]
            
            # Fit model
            if method == "linear":
                coeffs = np.polyfit(rna_valid, protein_valid, 1)
                # Predict for new data
                if gene_id in rna_expression.columns:
                    predictions[gene_id] = np.polyval(coeffs, rna_expression[gene_id].values)
            elif method == "log_linear":
                log_rna = np.log(rna_valid + 1)
                log_protein = np.log(protein_valid + 1)
                coeffs = np.polyfit(log_rna, log_protein, 1)
                if gene_id in rna_expression.columns:
                    log_rna_new = np.log(rna_expression[gene_id].values + 1)
                    predictions[gene_id] = np.exp(np.polyval(coeffs, log_rna_new)) - 1
            else:  # power_law
                log_rna = np.log(rna_valid + 1)
                log_protein = np.log(protein_valid + 1)
                coeffs = np.polyfit(log_rna, log_protein, 1)
                if gene_id in rna_expression.columns:
                    log_rna_new = np.log(rna_expression[gene_id].values + 1)
                    predictions[gene_id] = np.exp(np.polyval(coeffs, log_rna_new)) - 1
        
        if predictions:
            return pd.DataFrame(predictions, index=rna_expression.index)
    
    # Fallback to simple prediction
    return _simple_protein_prediction(rna_expression, method)


def _simple_protein_prediction(rna_expression: pd.DataFrame, method: str) -> pd.DataFrame:
    """Simple protein prediction without training data."""
    if method == "linear":
        # Assume linear relationship: protein ≈ 0.1 * RNA
        return rna_expression * 0.1
    elif method == "log_linear":
        # Log-linear: log(protein) ≈ log(RNA) - 2
        return np.exp(np.log(rna_expression + 1) - 2) - 1
    else:  # power_law
        # Power law: protein ≈ RNA^0.8 * 0.05
        return (rna_expression + 1) ** 0.8 * 0.05


def ribosome_profiling_integration(
    rna_expression: pd.DataFrame,
    ribo_profiling: pd.DataFrame,
    gene_annotations: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Integrate ribosome profiling data with RNA expression.
    
    Ribosome profiling (Ribo-seq) measures actively translated mRNAs.
    This function integrates Ribo-seq data with RNA-seq to identify
    translationally regulated genes.
    
    Args:
        rna_expression: RNA-seq expression data (samples x genes)
        ribo_profiling: Ribo-seq profiling data (samples x genes)
        gene_annotations: Optional gene annotations with 'gene_id', 'cds_length' columns
        
    Returns:
        DataFrame with integrated metrics:
        - 'gene_id': Gene identifier
        - 'rna_level': RNA expression level
        - 'ribo_level': Ribosome profiling level
        - 'translation_rate': Ribo/RNA ratio
        - 'translationally_regulated': Boolean indicating translational regulation
    """
    if rna_expression.empty or ribo_profiling.empty:
        return pd.DataFrame(columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"])
    
    # Align samples and genes
    common_samples = set(rna_expression.index) & set(ribo_profiling.index)
    common_genes = set(rna_expression.columns) & set(ribo_profiling.columns)
    
    if not common_samples or not common_genes:
        logger.warning("No common samples or genes between RNA and Ribo-seq data")
        return pd.DataFrame(columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"])
    
    rna_aligned = rna_expression.loc[list(common_samples), list(common_genes)]
    ribo_aligned = ribo_profiling.loc[list(common_samples), list(common_genes)]
    
    # Calculate mean levels and translation rates
    results = []
    
    for gene_id in common_genes:
        rna_mean = rna_aligned[gene_id].mean()
        ribo_mean = ribo_aligned[gene_id].mean()
        
        # Translation rate (Ribo/RNA ratio)
        if rna_mean > 0:
            translation_rate = ribo_mean / rna_mean
        else:
            translation_rate = 0.0
        
        # Identify translationally regulated genes (high Ribo, low RNA or vice versa)
        # Threshold: translation rate > 2x median or < 0.5x median
        translationally_regulated = False  # Would need to calculate median across all genes
        
        results.append({
            "gene_id": gene_id,
            "rna_level": float(rna_mean),
            "ribo_level": float(ribo_mean),
            "translation_rate": float(translation_rate),
            "translationally_regulated": translationally_regulated,
        })
    
    if not results:
        return pd.DataFrame(columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"])
    
    result_df = pd.DataFrame(results)
    
    # Calculate median translation rate for regulation detection
    median_rate = result_df["translation_rate"].median()
    result_df["translationally_regulated"] = (
        (result_df["translation_rate"] > 2 * median_rate) |
        (result_df["translation_rate"] < 0.5 * median_rate)
    )
    
    return result_df

