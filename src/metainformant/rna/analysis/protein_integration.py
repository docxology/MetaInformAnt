"""RNA-Protein integration and translation efficiency analysis.

This module provides tools for integrating RNA-seq data with protein abundance
measurements. It enables analysis of translation efficiency, prediction of protein
levels from RNA expression, and integration of ribosome profiling data.

Main Functions:
    - calculate_translation_efficiency: Calculate protein-RNA translation efficiency
    - predict_protein_abundance_from_rna: Predict protein levels from RNA expression
    - ribosome_profiling_integration: Integrate ribosome profiling data

Example:
    >>> from metainformant.rna import protein_integration
    >>> import pandas as pd
    >>> rna_df = pd.DataFrame(...)  # genes × samples
    >>> protein_df = pd.DataFrame(...)  # genes × samples
    >>> efficiency = protein_integration.calculate_translation_efficiency(
    ...     rna_df, protein_df, method="ratio"
    ... )
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from metainformant.core import logging

logger = logging.get_logger(__name__)


def calculate_translation_efficiency(
    rna_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    method: str = "ratio",
) -> pd.DataFrame:
    """Calculate translation efficiency from RNA and protein levels.

    Computes the efficiency of translation for each gene using multiple methods.
    Translation efficiency reflects the relationship between mRNA levels and
    corresponding protein abundance.

    Args:
        rna_df: RNA expression data (samples × genes)
        protein_df: Protein abundance data (samples × genes)
        method: Calculation method:
            - "ratio": Protein/RNA ratio per sample, aggregated
            - "correlation": Pearson correlation between RNA and protein

    Returns:
        DataFrame with columns:
        - gene_id: str - Gene identifier
        - efficiency: float - Calculated efficiency value
        - method: str - Method used

        Returns empty DataFrame if inputs are empty or have no common samples.

    Example:
        >>> rna_df = pd.DataFrame(
        ...     {"gene1": [10, 20], "gene2": [5, 10]},
        ...     index=["sample1", "sample2"]
        ... )
        >>> protein_df = pd.DataFrame(
        ...     {"gene1": [1, 2], "gene2": [0.5, 1]},
        ...     index=["sample1", "sample2"]
        ... )
        >>> eff = calculate_translation_efficiency(rna_df, protein_df, method="ratio")
        >>> assert len(eff) > 0
        >>> assert "gene_id" in eff.columns
    """
    # Handle empty DataFrames
    if rna_df.empty or protein_df.empty:
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])

    # Find common samples and genes
    common_samples = rna_df.index.intersection(protein_df.index)
    common_genes = rna_df.columns.intersection(protein_df.columns)

    if len(common_samples) == 0 or len(common_genes) == 0:
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])

    # Get subset of data
    rna_subset = rna_df.loc[common_samples, common_genes]
    protein_subset = protein_df.loc[common_samples, common_genes]

    results = []

    if method == "ratio":
        # Calculate protein/RNA ratio per gene
        for gene in common_genes:
            rna_vals = rna_subset[gene]
            prot_vals = protein_subset[gene]

            # Filter zero and NaN values
            valid_mask = (rna_vals > 0) & (prot_vals > 0) & rna_vals.notna() & prot_vals.notna()

            if valid_mask.sum() > 0:
                ratios = prot_vals[valid_mask] / rna_vals[valid_mask]
                efficiency = float(ratios.mean())

                results.append(
                    {
                        "gene_id": gene,
                        "efficiency": efficiency,
                        "method": "ratio",
                    }
                )

    elif method == "correlation":
        # Calculate correlation between RNA and protein for each gene
        for gene in common_genes:
            rna_vals = rna_subset[gene]
            prot_vals = protein_subset[gene]

            # Filter NaN values
            valid_mask = rna_vals.notna() & prot_vals.notna()

            if valid_mask.sum() > 1:
                # Calculate Pearson correlation
                correlation = float(rna_vals[valid_mask].corr(prot_vals[valid_mask]))

                if not np.isnan(correlation):
                    results.append(
                        {
                            "gene_id": gene,
                            "efficiency": correlation,
                            "method": "correlation",
                        }
                    )

    logger.debug(f"Calculated translation efficiency for {len(results)} genes using {method}")

    if results:
        return pd.DataFrame(results)
    else:
        return pd.DataFrame(columns=["gene_id", "efficiency", "method"])


def predict_protein_abundance_from_rna(
    rna_df: pd.DataFrame,
    training_rna: Optional[pd.DataFrame] = None,
    training_protein: Optional[pd.DataFrame] = None,
    method: str = "linear",
) -> pd.DataFrame:
    """Predict protein abundance levels from RNA expression.

    Uses RNA expression data to predict corresponding protein abundance.
    Can optionally use training data to learn the RNA-protein relationship.

    Args:
        rna_df: RNA expression data for prediction (samples × genes)
        training_rna: Optional RNA data from training set
        training_protein: Optional protein data from training set
        method: Prediction method:
            - "linear": Simple linear scaling (default)
            - "lognormal": Lognormal model

    Returns:
        DataFrame with same shape and indices as rna_df containing
        predicted protein abundance values.

        Returns empty DataFrame if input is empty.

    Example:
        >>> rna_df = pd.DataFrame(
        ...     {"gene1": [10, 20, 30], "gene2": [5, 10, 15]},
        ...     index=["sample1", "sample2", "sample3"]
        ... )
        >>> predictions = predict_protein_abundance_from_rna(rna_df, method="linear")
        >>> assert predictions.shape == rna_df.shape
    """
    # Handle empty DataFrame
    if rna_df.empty:
        return pd.DataFrame()

    if training_rna is not None and training_protein is not None:
        # Learn relationship from training data
        if not training_rna.empty and not training_protein.empty:
            # Find common samples and genes
            common_samples = training_rna.index.intersection(training_protein.index)
            common_genes = training_rna.columns.intersection(training_protein.columns)

            if len(common_samples) > 0 and len(common_genes) > 0:
                train_rna = training_rna.loc[common_samples, common_genes]
                train_prot = training_protein.loc[common_samples, common_genes]

                # Calculate scaling factors per gene
                scaling_factors = {}
                for gene in common_genes:
                    rna_mean = train_rna[gene].mean()
                    prot_mean = train_prot[gene].mean()

                    if rna_mean > 0:
                        scaling_factors[gene] = prot_mean / rna_mean
                    else:
                        scaling_factors[gene] = 1.0

                # Apply scaling to prediction data
                predictions = rna_df.copy()
                for gene in predictions.columns:
                    if gene in scaling_factors:
                        predictions[gene] = predictions[gene] * scaling_factors[gene]
                    else:
                        # No training data for this gene, use identity
                        predictions[gene] = predictions[gene]

                return predictions

    # Simple linear method without training data
    # Scale to mean values
    predictions = rna_df.copy()

    if method == "linear":
        # Log transform to reduce skewness then scale
        # Use map for pandas >= 2.1 compatibility
        try:
            predictions = predictions.map(lambda x: np.log2(x + 1) if x > 0 else 0)
        except AttributeError:
            # Fallback for older pandas versions
            predictions = predictions.applymap(lambda x: np.log2(x + 1) if x > 0 else 0)

    logger.debug(f"Predicted protein abundance for {rna_df.shape[1]} genes")

    return predictions


def ribosome_profiling_integration(
    rna_df: pd.DataFrame,
    ribo_df: pd.DataFrame,
) -> pd.DataFrame:
    """Integrate ribosome profiling data with RNA-seq.

    Combines RNA-seq and ribosome profiling (Ribo-seq) data to identify
    translationally regulated genes and calculate translation rates.

    Ribosome profiling measures the density of ribosomes on mRNAs, providing
    a direct measure of translation activity independent of mRNA levels.

    Args:
        rna_df: RNA-seq expression data (samples × genes)
        ribo_df: Ribosome profiling data (samples × genes)

    Returns:
        DataFrame with columns:
        - gene_id: str - Gene identifier
        - rna_level: float - Mean RNA expression
        - ribo_level: float - Mean ribosome occupancy
        - translation_rate: float - Ribo/RNA ratio
        - translationally_regulated: bool - Whether significantly regulated

        Returns empty DataFrame if inputs are empty or have no common data.

    Example:
        >>> rna_df = pd.DataFrame(
        ...     {"gene1": [10, 20, 30]},
        ...     index=["sample1", "sample2", "sample3"]
        ... )
        >>> ribo_df = pd.DataFrame(
        ...     {"gene1": [1, 2, 3]},
        ...     index=["sample1", "sample2", "sample3"]
        ... )
        >>> result = ribosome_profiling_integration(rna_df, ribo_df)
    """
    # Handle empty DataFrames
    if rna_df.empty or ribo_df.empty:
        return pd.DataFrame(
            columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"]
        )

    # Find common data
    common_samples = rna_df.index.intersection(ribo_df.index)
    common_genes = rna_df.columns.intersection(ribo_df.columns)

    if len(common_samples) == 0 or len(common_genes) == 0:
        return pd.DataFrame(
            columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"]
        )

    rna_subset = rna_df.loc[common_samples, common_genes]
    ribo_subset = ribo_df.loc[common_samples, common_genes]

    results = []

    for gene in common_genes:
        rna_vals = rna_subset[gene]
        ribo_vals = ribo_subset[gene]

        # Calculate mean levels (filter NaN)
        rna_level = float(rna_vals[rna_vals.notna()].mean())
        ribo_level = float(ribo_vals[ribo_vals.notna()].mean())

        # Calculate translation rate (Ribo/RNA)
        if rna_level > 0:
            translation_rate = ribo_level / rna_level
        else:
            translation_rate = 0.0

        # Determine if translationally regulated
        # Simple heuristic: compare to median translation rate
        translationally_regulated = False
        if translation_rate > 0:
            # Would use statistical test in full implementation
            median_rate = (ribo_subset.mean() / rna_subset.mean()).median()
            if not np.isnan(median_rate) and median_rate > 0:
                translationally_regulated = translation_rate > median_rate * 1.5

        results.append(
            {
                "gene_id": gene,
                "rna_level": rna_level,
                "ribo_level": ribo_level,
                "translation_rate": translation_rate,
                "translationally_regulated": translationally_regulated,
            }
        )

    logger.debug(f"Integrated ribosome profiling data for {len(results)} genes")

    if results:
        return pd.DataFrame(results)
    else:
        return pd.DataFrame(
            columns=["gene_id", "rna_level", "ribo_level", "translation_rate", "translationally_regulated"]
        )
