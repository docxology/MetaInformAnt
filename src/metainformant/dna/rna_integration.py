"""DNA-RNA integration utilities for variant-to-transcript mapping and eQTL analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def map_variants_to_transcripts(
    variants: pd.DataFrame,
    transcript_annotations: pd.DataFrame,
    window_size: int = 5000,
) -> pd.DataFrame:
    """Map genetic variants to nearby transcripts.
    
    Identifies transcripts that overlap or are near genetic variants,
    useful for understanding variant effects on gene expression.
    
    Args:
        variants: DataFrame with 'chrom', 'pos' columns (and optionally 'ref', 'alt')
        transcript_annotations: DataFrame with 'chrom', 'start', 'end', 'transcript_id' columns
        window_size: Window size around transcript boundaries to consider (default: 5kb)
        
    Returns:
        DataFrame with variant-transcript mappings:
        - 'variant_chrom', 'variant_pos': Variant coordinates
        - 'transcript_id': Transcript identifier
        - 'distance': Distance from variant to transcript (0 if overlapping)
        - 'overlaps': Boolean indicating if variant overlaps transcript
    """
    if variants.empty or transcript_annotations.empty:
        return pd.DataFrame(columns=["variant_chrom", "variant_pos", "transcript_id", "distance", "overlaps"])
    
    mappings = []
    
    for _, variant in variants.iterrows():
        var_chrom = variant["chrom"]
        var_pos = variant["pos"]
        
        # Find transcripts on same chromosome
        chrom_transcripts = transcript_annotations[transcript_annotations["chrom"] == var_chrom]
        
        for _, transcript in chrom_transcripts.iterrows():
            trans_start = transcript["start"]
            trans_end = transcript["end"]
            trans_id = transcript["transcript_id"]
            
            # Check if variant overlaps transcript
            if trans_start <= var_pos <= trans_end:
                distance = 0
                overlaps = True
            else:
                # Calculate distance to transcript
                if var_pos < trans_start:
                    distance = trans_start - var_pos
                else:
                    distance = var_pos - trans_end
                overlaps = False
            
            # Include if within window
            if distance <= window_size:
                mappings.append({
                    "variant_chrom": var_chrom,
                    "variant_pos": var_pos,
                    "transcript_id": trans_id,
                    "distance": distance,
                    "overlaps": overlaps,
                })
    
    if not mappings:
        return pd.DataFrame(columns=["variant_chrom", "variant_pos", "transcript_id", "distance", "overlaps"])
    
    return pd.DataFrame(mappings)


def predict_variant_effect_on_expression(
    variant: Dict[str, Any],
    transcript_sequence: str,
    variant_position_in_transcript: int | None = None,
) -> Dict[str, Any]:
    """Predict effect of variant on transcript expression.
    
    Simplified prediction based on variant location (promoter, coding, UTR, etc.).
    More sophisticated predictions would require additional tools.
    
    Args:
        variant: Dictionary with variant information ('chrom', 'pos', 'ref', 'alt')
        transcript_sequence: Transcript sequence (for coding variants)
        variant_position_in_transcript: Position of variant within transcript (if known)
        
    Returns:
        Dictionary with predicted effects:
        - 'effect_type': Type of predicted effect (promoter, coding, UTR, intergenic)
        - 'severity': Predicted severity (high, medium, low, none)
        - 'notes': Additional notes about prediction
    """
    # Simplified prediction based on variant location
    # In practice, would use tools like VEP, SnpEff, etc.
    
    effect_type = "intergenic"
    severity = "none"
    notes = "Simplified prediction - use specialized tools for accurate annotation"
    
    if variant_position_in_transcript is not None:
        if variant_position_in_transcript < 1000:  # Upstream/promoter region
            effect_type = "promoter"
            severity = "high"
        elif variant_position_in_transcript < len(transcript_sequence) - 500:  # Coding region
            effect_type = "coding"
            severity = "high"
        else:  # 3' UTR
            effect_type = "utr_3"
            severity = "medium"
    
    return {
        "effect_type": effect_type,
        "severity": severity,
        "notes": notes,
    }


def eqtl_analysis(
    genotypes: pd.DataFrame,
    expression: pd.DataFrame,
    variant_positions: pd.DataFrame,
    gene_positions: pd.DataFrame,
    window_size: int = 1000000,
    min_samples: int = 10,
) -> pd.DataFrame:
    """Perform expression quantitative trait locus (eQTL) analysis.
    
    Identifies genetic variants associated with gene expression levels.
    
    Args:
        genotypes: DataFrame with samples as rows and variants as columns (0, 1, 2 encoding)
        expression: DataFrame with samples as rows and genes as columns
        variant_positions: DataFrame with 'variant_id', 'chrom', 'pos' columns
        gene_positions: DataFrame with 'gene_id', 'chrom', 'start', 'end' columns
        window_size: Window size around genes to search for variants (default: 1Mb)
        min_samples: Minimum number of samples required for association test
        
    Returns:
        DataFrame with eQTL associations:
        - 'gene_id': Gene identifier
        - 'variant_id': Variant identifier
        - 'distance': Distance from variant to gene
        - 'correlation': Correlation between genotype and expression
        - 'p_value': Statistical significance (if scipy available)
    """
    if genotypes.empty or expression.empty:
        return pd.DataFrame(columns=["gene_id", "variant_id", "distance", "correlation", "p_value"])
    
    # Ensure samples are aligned
    common_samples = set(genotypes.index) & set(expression.index)
    if len(common_samples) < min_samples:
        logger.warning(f"Only {len(common_samples)} common samples, need at least {min_samples}")
        return pd.DataFrame(columns=["gene_id", "variant_id", "distance", "correlation", "p_value"])
    
    genotypes_aligned = genotypes.loc[list(common_samples)]
    expression_aligned = expression.loc[list(common_samples)]
    
    associations = []
    
    # For each gene, find nearby variants
    for _, gene in gene_positions.iterrows():
        gene_id = gene["gene_id"]
        gene_chrom = gene["chrom"]
        gene_start = gene["start"]
        gene_end = gene["end"]
        gene_center = (gene_start + gene_end) // 2
        
        if gene_id not in expression_aligned.columns:
            continue
        
        gene_expr = expression_aligned[gene_id].values
        
        # Find variants within window
        nearby_variants = variant_positions[
            (variant_positions["chrom"] == gene_chrom) &
            (variant_positions["pos"] >= gene_center - window_size) &
            (variant_positions["pos"] <= gene_center + window_size)
        ]
        
        if nearby_variants.empty:
            continue
        
        # Test association with each variant
        for _, variant in nearby_variants.iterrows():
            variant_id = variant["variant_id"]
            
            if variant_id not in genotypes_aligned.columns:
                continue
            
            variant_genos = genotypes_aligned[variant_id].values
            
            # Calculate correlation
            try:
                correlation = np.corrcoef(gene_expr, variant_genos)[0, 1]
                
                if np.isnan(correlation):
                    continue
                
                # Calculate p-value if scipy available
                try:
                    from scipy.stats import pearsonr
                    _, p_value = pearsonr(gene_expr, variant_genos)
                except ImportError:
                    p_value = np.nan
                
                distance = abs(variant["pos"] - gene_center)
                
                associations.append({
                    "gene_id": gene_id,
                    "variant_id": variant_id,
                    "distance": int(distance),
                    "correlation": float(correlation),
                    "p_value": float(p_value) if not np.isnan(p_value) else None,
                })
            except Exception as e:
                logger.debug(f"Error calculating correlation for {gene_id}-{variant_id}: {e}")
                continue
    
    if not associations:
        return pd.DataFrame(columns=["gene_id", "variant_id", "distance", "correlation", "p_value"])
    
    result_df = pd.DataFrame(associations)
    result_df = result_df.sort_values("p_value" if "p_value" in result_df.columns else "correlation", ascending=True)
    
    return result_df

