"""ChIP-seq peak calling and analysis utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from metainformant.core.io import dump_json
from metainformant.core.logging import get_logger
from metainformant.core.paths import expand_and_resolve, prepare_file_path

logger = get_logger(__name__)


def call_peaks_simple(
    signal: pd.DataFrame,
    chrom: str,
    threshold: float | None = None,
    min_width: int = 100,
    min_gap: int = 200,
) -> pd.DataFrame:
    """Call peaks from ChIP-seq signal using simple threshold-based method.
    
    Identifies regions where signal exceeds threshold and merges nearby peaks.
    
    Args:
        signal: DataFrame with columns ['chrom', 'start', 'end', 'value']
        chrom: Chromosome to analyze
        threshold: Signal threshold for peak calling. If None, uses median + 2*MAD
        min_width: Minimum peak width in base pairs
        min_gap: Minimum gap between peaks in base pairs
        
    Returns:
        DataFrame with columns ['chrom', 'start', 'end', 'peak_id', 'max_value', 'sum_value']
    """
    # Filter to chromosome
    chrom_data = signal[signal["chrom"] == chrom].copy()
    if chrom_data.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "max_value", "sum_value"])
    
    # Sort by start position
    chrom_data = chrom_data.sort_values("start").reset_index(drop=True)
    
    # Determine threshold if not provided
    if threshold is None:
        values = chrom_data["value"].values
        median = np.median(values)
        mad = np.median(np.abs(values - median))
        threshold = median + 2 * mad
        logger.info(f"Using adaptive threshold: {threshold:.2f} (median + 2*MAD)")
    
    # Find regions above threshold
    above_threshold = chrom_data["value"] >= threshold
    
    if not above_threshold.any():
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "max_value", "sum_value"])
    
    # Identify contiguous regions
    peaks = []
    in_peak = False
    peak_start = None
    peak_end = None
    peak_max = 0.0
    peak_sum = 0.0
    
    for idx, row in chrom_data.iterrows():
        if above_threshold.iloc[idx]:
            if not in_peak:
                # Start new peak
                in_peak = True
                peak_start = row["start"]
                peak_max = row["value"]
                peak_sum = row["value"]
            else:
                # Extend peak
                peak_max = max(peak_max, row["value"])
                peak_sum += row["value"]
            peak_end = row["end"]
        else:
            if in_peak:
                # End current peak
                peak_width = peak_end - peak_start
                if peak_width >= min_width:
                    peaks.append({
                        "chrom": chrom,
                        "start": peak_start,
                        "end": peak_end,
                        "max_value": peak_max,
                        "sum_value": peak_sum,
                    })
                in_peak = False
                peak_start = None
                peak_end = None
                peak_max = 0.0
                peak_sum = 0.0
    
    # Handle peak at end of chromosome
    if in_peak and peak_start is not None and peak_end is not None:
        peak_width = peak_end - peak_start
        if peak_width >= min_width:
            peaks.append({
                "chrom": chrom,
                "start": peak_start,
                "end": peak_end,
                "max_value": peak_max,
                "sum_value": peak_sum,
            })
    
    if not peaks:
        return pd.DataFrame(columns=["chrom", "start", "end", "peak_id", "max_value", "sum_value"])
    
    peaks_df = pd.DataFrame(peaks)
    
    # Merge nearby peaks
    merged_peaks = []
    for idx, peak in peaks_df.iterrows():
        if not merged_peaks:
            merged_peaks.append(peak.to_dict())
        else:
            last_peak = merged_peaks[-1]
            gap = peak["start"] - last_peak["end"]
            if gap <= min_gap:
                # Merge with previous peak
                last_peak["end"] = peak["end"]
                last_peak["max_value"] = max(last_peak["max_value"], peak["max_value"])
                last_peak["sum_value"] += peak["sum_value"]
            else:
                merged_peaks.append(peak.to_dict())
    
    # Add peak IDs
    result_df = pd.DataFrame(merged_peaks)
    result_df["peak_id"] = [f"{chrom}_peak_{i+1}" for i in range(len(result_df))]
    
    return result_df[["chrom", "start", "end", "peak_id", "max_value", "sum_value"]]


def analyze_peak_overlap(
    peaks1: pd.DataFrame,
    peaks2: pd.DataFrame,
) -> dict[str, Any]:
    """Analyze overlap between two sets of peaks.
    
    Args:
        peaks1: First peak set (must have 'chrom', 'start', 'end' columns)
        peaks2: Second peak set (must have 'chrom', 'start', 'end' columns)
        
    Returns:
        Dictionary with overlap statistics
    """
    if peaks1.empty or peaks2.empty:
        return {
            "n_peaks1": 0,
            "n_peaks2": 0,
            "n_overlapping": 0,
            "overlap_fraction": 0.0,
        }
    
    # Calculate overlap for each peak pair
    overlaps = []
    for _, peak1 in peaks1.iterrows():
        chrom1 = peak1["chrom"]
        start1 = peak1["start"]
        end1 = peak1["end"]
        
        # Find overlapping peaks in peaks2
        chrom_peaks2 = peaks2[peaks2["chrom"] == chrom1]
        for _, peak2 in chrom_peaks2.iterrows():
            start2 = peak2["start"]
            end2 = peak2["end"]
            
            # Check overlap
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                peak1_length = end1 - start1
                peak2_length = end2 - start2
                
                overlaps.append({
                    "peak1_start": start1,
                    "peak1_end": end1,
                    "peak2_start": start2,
                    "peak2_end": end2,
                    "overlap_length": overlap_length,
                    "overlap_fraction1": overlap_length / peak1_length if peak1_length > 0 else 0.0,
                    "overlap_fraction2": overlap_length / peak2_length if peak2_length > 0 else 0.0,
                })
    
    n_overlapping = len(set((o["peak1_start"], o["peak1_end"]) for o in overlaps))
    
    return {
        "n_peaks1": len(peaks1),
        "n_peaks2": len(peaks2),
        "n_overlapping": n_overlapping,
        "overlap_fraction": n_overlapping / len(peaks1) if len(peaks1) > 0 else 0.0,
        "overlaps": overlaps[:100],  # Limit to first 100 for output
    }


def peak_enrichment_analysis(
    peaks: pd.DataFrame,
    gene_annotations: pd.DataFrame | None = None,
    window_size: int = 5000,
) -> dict[str, Any]:
    """Analyze peak enrichment near genomic features.
    
    Args:
        peaks: Peak DataFrame with 'chrom', 'start', 'end' columns
        gene_annotations: Optional gene annotations with 'chrom', 'start', 'end', 'gene_id' columns
        window_size: Window size around features for enrichment analysis
        
    Returns:
        Dictionary with enrichment statistics
    """
    if peaks.empty:
        return {
            "n_peaks": 0,
            "total_peak_length": 0,
            "genes_with_peaks": 0,
        }
    
    total_peak_length = (peaks["end"] - peaks["start"]).sum()
    
    result = {
        "n_peaks": len(peaks),
        "total_peak_length": int(total_peak_length),
    }
    
    if gene_annotations is not None and not gene_annotations.empty:
        # Count peaks near genes
        genes_with_peaks = 0
        for _, gene in gene_annotations.iterrows():
            chrom = gene["chrom"]
            gene_start = gene["start"]
            gene_end = gene["end"]
            
            # Find peaks in window around gene
            chrom_peaks = peaks[peaks["chrom"] == chrom]
            nearby_peaks = chrom_peaks[
                ((chrom_peaks["start"] >= gene_start - window_size) & (chrom_peaks["start"] <= gene_end + window_size)) |
                ((chrom_peaks["end"] >= gene_start - window_size) & (chrom_peaks["end"] <= gene_end + window_size)) |
                ((chrom_peaks["start"] <= gene_start) & (chrom_peaks["end"] >= gene_end))
            ]
            
            if not nearby_peaks.empty:
                genes_with_peaks += 1
        
        result["genes_with_peaks"] = genes_with_peaks
        result["n_genes_analyzed"] = len(gene_annotations)
        result["enrichment_fraction"] = genes_with_peaks / len(gene_annotations) if len(gene_annotations) > 0 else 0.0
    
    return result

