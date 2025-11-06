"""ATAC-seq chromatin accessibility analysis utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from metainformant.core.io import dump_json
from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve, prepare_file_path
from metainformant.epigenome.tracks import read_bedgraph

logger = get_logger(__name__)


def identify_accessible_regions(
    signal: pd.DataFrame,
    chrom: str,
    threshold: float | None = None,
    min_width: int = 150,
) -> pd.DataFrame:
    """Identify accessible chromatin regions from ATAC-seq signal.
    
    ATAC-seq typically produces narrow peaks (~150bp) representing nucleosome-free
    regions. This function identifies these accessible regions.
    
    Args:
        signal: DataFrame with columns ['chrom', 'start', 'end', 'value']
        chrom: Chromosome to analyze
        threshold: Signal threshold. If None, uses 75th percentile
        min_width: Minimum width for accessible region (default: 150bp for ATAC-seq)
        
    Returns:
        DataFrame with accessible regions ['chrom', 'start', 'end', 'peak_id', 'value']
    """
    # Filter to chromosome
    chrom_data = signal[signal["chrom"] == chrom].copy()
    if chrom_data.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "value"])
    
    # Sort by start position
    chrom_data = chrom_data.sort_values("start").reset_index(drop=True)
    
    # Determine threshold if not provided (use 75th percentile for ATAC-seq)
    if threshold is None:
        threshold = chrom_data["value"].quantile(0.75)
        logger.info(f"Using 75th percentile threshold: {threshold:.2f}")
    
    # Find regions above threshold
    above_threshold = chrom_data["value"] >= threshold
    
    if not above_threshold.any():
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "value"])
    
    # Identify accessible regions
    accessible_regions = []
    in_region = False
    region_start = None
    region_end = None
    region_max = 0.0
    
    for idx, row in chrom_data.iterrows():
        if above_threshold.iloc[idx]:
            if not in_region:
                # Start new region
                in_region = True
                region_start = row["start"]
                region_max = row["value"]
            else:
                # Extend region
                region_max = max(region_max, row["value"])
            region_end = row["end"]
        else:
            if in_region:
                # End current region
                region_width = region_end - region_start
                if region_width >= min_width:
                    accessible_regions.append({
                        "chrom": chrom,
                        "start": region_start,
                        "end": region_end,
                        "value": region_max,
                    })
                in_region = False
                region_start = None
                region_end = None
                region_max = 0.0
    
    # Handle region at end
    if in_region and region_start is not None and region_end is not None:
        region_width = region_end - region_start
        if region_width >= min_width:
            accessible_regions.append({
                "chrom": chrom,
                "start": region_start,
                "end": region_end,
                "value": region_max,
            })
    
    if not accessible_regions:
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "value"])
    
    result_df = pd.DataFrame(accessible_regions)
    result_df["peak_id"] = [f"{chrom}_accessible_{i+1}" for i in range(len(result_df))]
    
    return result_df[["chrom", "start", "end", "peak_id", "value"]]


def calculate_accessibility_scores(
    signal: pd.DataFrame,
    regions: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate accessibility scores for genomic regions.
    
    Args:
        signal: ATAC-seq signal DataFrame with 'chrom', 'start', 'end', 'value'
        regions: Genomic regions DataFrame with 'chrom', 'start', 'end'
        
    Returns:
        DataFrame with regions and their accessibility scores
    """
    if signal.empty or regions.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "accessibility_score"])
    
    scores = []
    for _, region in regions.iterrows():
        chrom = region["chrom"]
        start = region["start"]
        end = region["end"]
        
        # Find overlapping signal
        chrom_signal = signal[signal["chrom"] == chrom]
        overlapping = chrom_signal[
            ((chrom_signal["start"] < end) & (chrom_signal["end"] > start))
        ]
        
        if overlapping.empty:
            score = 0.0
        else:
            # Calculate weighted average accessibility
            total_signal = 0.0
            total_length = 0
            
            for _, sig in overlapping.iterrows():
                overlap_start = max(start, sig["start"])
                overlap_end = min(end, sig["end"])
                overlap_length = overlap_end - overlap_start
                
                total_signal += sig["value"] * overlap_length
                total_length += overlap_length
            
            score = total_signal / total_length if total_length > 0 else 0.0
        
        scores.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "accessibility_score": score,
        })
    
    return pd.DataFrame(scores)


def compare_accessibility(
    signal1: pd.DataFrame,
    signal2: pd.DataFrame,
    regions: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Compare accessibility between two conditions.
    
    Args:
        signal1: ATAC-seq signal from condition 1
        signal2: ATAC-seq signal from condition 2
        regions: Optional regions to compare. If None, uses all overlapping regions
        
    Returns:
        Dictionary with comparison statistics
    """
    if signal1.empty or signal2.empty:
        return {
            "n_regions": 0,
            "mean_fold_change": 0.0,
            "n_differential": 0,
        }
    
    # Calculate accessibility scores for each condition
    if regions is None:
        # Use all chromosomes present in both signals
        chroms1 = set(signal1["chrom"].unique())
        chroms2 = set(signal2["chrom"].unique())
        common_chroms = chroms1 & chroms2
        
        # Create regions from signal (simplified - use signal bins as regions)
        regions = []
        for chrom in common_chroms:
            chrom_sig1 = signal1[signal1["chrom"] == chrom].sort_values("start")
            chrom_sig2 = signal2[signal2["chrom"] == chrom].sort_values("start")
            
            # Use union of regions
            all_starts = sorted(set(chrom_sig1["start"].tolist() + chrom_sig2["start"].tolist()))
            all_ends = sorted(set(chrom_sig1["end"].tolist() + chrom_sig2["end"].tolist()))
            
            # Create regions (simplified)
            for i in range(min(len(all_starts), len(all_ends))):
                regions.append({
                    "chrom": chrom,
                    "start": all_starts[i] if i < len(all_starts) else all_starts[-1],
                    "end": all_ends[i] if i < len(all_ends) else all_ends[-1],
                })
        
        regions = pd.DataFrame(regions)
    
    scores1 = calculate_accessibility_scores(signal1, regions)
    scores2 = calculate_accessibility_scores(signal2, regions)
    
    # Merge scores
    merged = pd.merge(
        scores1,
        scores2,
        on=["chrom", "start", "end"],
        suffixes=("_cond1", "_cond2"),
    )
    
    # Calculate fold changes (avoid division by zero)
    merged["fold_change"] = np.where(
        merged["accessibility_score_cond2"] > 0,
        merged["accessibility_score_cond1"] / merged["accessibility_score_cond2"],
        np.where(merged["accessibility_score_cond1"] > 0, np.inf, 1.0),
    )
    
    # Log2 fold change
    merged["log2_fold_change"] = np.log2(np.maximum(merged["fold_change"], 1e-10))
    
    # Identify differential regions (fold change > 2 or < 0.5)
    merged["differential"] = (merged["fold_change"] > 2.0) | (merged["fold_change"] < 0.5)
    
    return {
        "n_regions": len(merged),
        "mean_fold_change": float(merged["fold_change"].mean()),
        "mean_log2_fold_change": float(merged["log2_fold_change"].mean()),
        "n_differential": int(merged["differential"].sum()),
        "differential_fraction": float(merged["differential"].mean()),
    }

