"""Epigenome analysis workflow orchestration."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from metainformant.core import io, logging, paths
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def run_methylation_workflow(
    methylation_file: Path,
    output_dir: Path,
    compute_beta: bool = True,
    differential_analysis: bool = False,
    condition2_file: Path | None = None,
) -> Dict[str, Any]:
    """Run complete DNA methylation analysis workflow.
    
    Args:
        methylation_file: Path to CpG methylation table
        output_dir: Output directory for results
        compute_beta: If True, compute beta values
        differential_analysis: If True, perform differential methylation analysis
        condition2_file: Optional second condition file for differential analysis
        
    Returns:
        Dictionary with workflow results and output paths
    """
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting methylation workflow: {methylation_file}")
    
    # Load methylation data
    from metainformant.epigenome import load_cpg_table, compute_beta_values, summarize_beta_by_chromosome
    
    cpg_df = load_cpg_table(methylation_file)
    logger.info(f"Loaded {len(cpg_df)} CpG sites")
    
    results = {
        "input_file": str(methylation_file),
        "n_sites": len(cpg_df),
    }
    
    # Compute beta values
    if compute_beta:
        cpg_df = compute_beta_values(cpg_df)
        beta_output = output_dir / "beta_values.tsv"
        cpg_df.to_csv(beta_output, sep="\t", index=False)
        results["beta_file"] = str(beta_output)
        logger.info(f"Computed beta values, saved to {beta_output}")
    
    # Summarize by chromosome
    summary = summarize_beta_by_chromosome(cpg_df)
    summary_output = output_dir / "chromosome_summary.json"
    io.dump_json(summary.to_dict(), summary_output)
    results["summary_file"] = str(summary_output)
    
    # Differential analysis if requested
    if differential_analysis and condition2_file:
        from metainformant.epigenome import differential_methylation
        
        condition2_df = load_cpg_table(condition2_file)
        if compute_beta:
            condition2_df = compute_beta_values(condition2_df)
        
        diff_results = differential_methylation(cpg_df, condition2_df)
        diff_output = output_dir / "differential_methylation.tsv"
        diff_results.to_csv(diff_output, sep="\t", index=False)
        results["differential_file"] = str(diff_output)
        results["n_differential_sites"] = len(diff_results)
        logger.info(f"Identified {len(diff_results)} differentially methylated sites")
    
    results["status"] = "completed"
    return results


def run_chipseq_workflow(
    signal_file: Path,
    output_dir: Path,
    chrom: str,
    threshold: float | None = None,
) -> Dict[str, Any]:
    """Run ChIP-seq peak calling workflow.
    
    Args:
        signal_file: Path to bedGraph signal file
        output_dir: Output directory for results
        chrom: Chromosome to analyze
        threshold: Optional signal threshold for peak calling
        
    Returns:
        Dictionary with workflow results
    """
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting ChIP-seq workflow: {signal_file}")
    
    from metainformant.epigenome import read_bedgraph
    from metainformant.epigenome.chipseq import call_peaks_simple
    
    # Load signal
    signal = read_bedgraph(signal_file)
    logger.info(f"Loaded signal data: {len(signal)} intervals")
    
    # Call peaks
    peaks = call_peaks_simple(signal, chrom, threshold=threshold)
    
    # Save results
    peaks_output = output_dir / f"{chrom}_peaks.tsv"
    peaks.to_csv(peaks_output, sep="\t", index=False)
    
    results = {
        "input_file": str(signal_file),
        "chromosome": chrom,
        "n_peaks": len(peaks),
        "peaks_file": str(peaks_output),
        "status": "completed",
    }
    
    logger.info(f"Called {len(peaks)} peaks, saved to {peaks_output}")
    
    return results


def run_atacseq_workflow(
    signal_file: Path,
    output_dir: Path,
    chrom: str,
    threshold: float | None = None,
) -> Dict[str, Any]:
    """Run ATAC-seq accessibility analysis workflow.
    
    Args:
        signal_file: Path to bedGraph signal file
        output_dir: Output directory for results
        chrom: Chromosome to analyze
        threshold: Optional signal threshold
        
    Returns:
        Dictionary with workflow results
    """
    paths.ensure_directory(output_dir)
    
    logger.info(f"Starting ATAC-seq workflow: {signal_file}")
    
    from metainformant.epigenome import read_bedgraph
    from metainformant.epigenome.atac import identify_accessible_regions
    
    # Load signal
    signal = read_bedgraph(signal_file)
    logger.info(f"Loaded signal data: {len(signal)} intervals")
    
    # Identify accessible regions
    accessible = identify_accessible_regions(signal, chrom, threshold=threshold)
    
    # Save results
    accessible_output = output_dir / f"{chrom}_accessible_regions.tsv"
    accessible.to_csv(accessible_output, sep="\t", index=False)
    
    results = {
        "input_file": str(signal_file),
        "chromosome": chrom,
        "n_accessible_regions": len(accessible),
        "accessible_regions_file": str(accessible_output),
        "status": "completed",
    }
    
    logger.info(f"Identified {len(accessible)} accessible regions, saved to {accessible_output}")
    
    return results

