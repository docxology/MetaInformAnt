"""Regional locus visualization for GWAS.

This module provides detailed views of specific genomic regions,
including LD structure, gene annotations, and recombination rates.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

from ..core.io import ensure_directory

logger = logging.getLogger(__name__)

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


def regional_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    chrom: str,
    start: int,
    end: int,
    lead_snp_pos: int | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Regional association plot around a locus.
    
    Shows association signal in a genomic window, useful for
    fine-mapping and identifying candidate causal variants.
    
    Args:
        results: Association results or path
        output_path: Output path
        chrom: Chromosome
        start: Region start position
        end: Region end position
        lead_snp_pos: Position of lead SNP to highlight
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"regional_plot: Plotting {chrom}:{start}-{end}")
    
    # Load data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Extract regional data
    positions, pvalues = [], []
    for r in results_list:
        try:
            r_chrom = str(r.get("CHROM", ""))
            r_pos = int(r.get("POS", 0))
            r_pval = float(r.get("p_value", 1.0))
            
            if r_chrom == chrom and start <= r_pos <= end and 0 < r_pval <= 1:
                positions.append(r_pos)
                pvalues.append(-math.log10(r_pval))
        except (ValueError, TypeError):
            continue
    
    if not positions:
        return {"status": "failed", "error": f"No variants in region {chrom}:{start}-{end}"}
    
    positions = np.array(positions)
    pvalues = np.array(pvalues)
    
    # Plot
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Scatter plot
    ax.scatter(positions, pvalues, c='#2E86AB', s=20, alpha=0.7,
              edgecolors='black', linewidths=0.5)
    
    # Highlight lead SNP
    if lead_snp_pos:
        lead_idx = np.abs(positions - lead_snp_pos).argmin()
        ax.scatter([positions[lead_idx]], [pvalues[lead_idx]],
                  c='red', s=100, marker='D', edgecolors='darkred',
                  linewidths=2, zorder=10, label='Lead SNP')
    
    # Labels
    ax.set_xlabel(f"Position on {chrom} (bp)", fontsize=12, fontweight="bold")
    ax.set_ylabel("-log₁₀(p-value)", fontsize=12, fontweight="bold")
    
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(f"Regional Plot: {chrom}:{start:,}-{end:,} ({len(positions)} variants)",
                    fontsize=14, fontweight="bold")
    
    ax.grid(True, alpha=0.3, linestyle=':')
    if lead_snp_pos:
        ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"regional_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "region": f"{chrom}:{start}-{end}",
        "num_variants": len(positions),
        "max_neg_log_p": float(pvalues.max()),
    }


def regional_ld_plot(
    vcf_path: Path,
    output_path: str | Path,
    *,
    chrom: str,
    start: int,
    end: int,
    lead_snp_pos: int,
    title: str | None = None,
) -> dict[str, Any]:
    """LD structure around lead SNP.
    
    Shows linkage disequilibrium (r²) between lead SNP and surrounding variants.
    Useful for identifying LD blocks and tagging variants.
    
    Args:
        vcf_path: Path to VCF file
        output_path: Output path
        chrom: Chromosome
        start: Region start
        end: Region end
        lead_snp_pos: Lead SNP position
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"regional_ld_plot: Computing LD around {chrom}:{lead_snp_pos}")
    
    # This requires genotype data parsing and LD calculation
    # For genome-scale data, use PLINK or specialized LD tools
    
    return {
        "status": "skipped",
        "message": "Regional LD calculation requires genotype parsing",
        "recommendation": "Use PLINK --r2 or LDlinkR for LD computation",
        "command_example": f"plink --vcf {vcf_path} --chr {chrom} --from-bp {start} --to-bp {end} --r2 --ld-window-r2 0 --out ld_region",
    }


def gene_annotation_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    gff_path: Path | None = None,
    chrom: str,
    start: int,
    end: int,
    title: str | None = None,
) -> dict[str, Any]:
    """Variants overlaid with gene annotations.
    
    Shows association signal with gene structures (exons, introns).
    Helps identify variants in regulatory or coding regions.
    
    Args:
        results: Association results or path
        output_path: Output path
        gff_path: Path to GFF3 annotation file
        chrom: Chromosome
        start: Region start
        end: Region end
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"gene_annotation_plot: Plotting genes in {chrom}:{start}-{end}")
    
    # Load association data
    if isinstance(results, (Path, str)):
        if PANDAS_AVAILABLE:
            df = pd.read_csv(results, sep="\t")
            results_list = df.to_dict("records")
        else:
            from ..core.io import read_tsv
            data = read_tsv(results)
            header = data[0]
            results_list = [{header[i]: row[i] for i in range(len(header))} for row in data[1:]]
    else:
        results_list = results
    
    # Extract regional data
    positions, pvalues = [], []
    for r in results_list:
        try:
            r_chrom = str(r.get("CHROM", ""))
            r_pos = int(r.get("POS", 0))
            r_pval = float(r.get("p_value", 1.0))
            
            if r_chrom == chrom and start <= r_pos <= end and 0 < r_pval <= 1:
                positions.append(r_pos)
                pvalues.append(-math.log10(r_pval))
        except (ValueError, TypeError):
            continue
    
    if not positions:
        return {"status": "failed", "error": f"No variants in region"}
    
    # Create two-panel plot
    fig, (ax_assoc, ax_genes) = plt.subplots(2, 1, figsize=(14, 8),
                                             gridspec_kw={'height_ratios': [3, 1]},
                                             sharex=True)
    
    # Association panel
    ax_assoc.scatter(positions, pvalues, c='#2E86AB', s=20, alpha=0.7)
    ax_assoc.set_ylabel("-log₁₀(p-value)", fontsize=11, fontweight="bold")
    ax_assoc.grid(True, alpha=0.3, linestyle=':')
    
    # Gene annotation panel (simplified - full implementation requires GFF parsing)
    if gff_path and gff_path.exists():
        ax_genes.text((start + end) / 2, 0.5, "Gene annotations require GFF parsing",
                     ha='center', va='center', fontsize=10, style='italic')
    else:
        ax_genes.text((start + end) / 2, 0.5, "No gene annotations provided",
                     ha='center', va='center', fontsize=10, style='italic', color='gray')
    
    ax_genes.set_xlim(start, end)
    ax_genes.set_ylim(0, 1)
    ax_genes.set_xlabel(f"Position on {chrom} (bp)", fontsize=11, fontweight="bold")
    ax_genes.set_yticks([])
    ax_genes.spines['left'].set_visible(False)
    ax_genes.spines['top'].set_visible(False)
    ax_genes.spines['right'].set_visible(False)
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold")
    else:
        fig.suptitle(f"Regional Association with Genes: {chrom}:{start:,}-{end:,}",
                    fontsize=14, fontweight="bold")
    
    plt.tight_layout()
    
    # Save
    output_path_obj = Path(output_path)
    ensure_directory(output_path_obj.parent)
    plt.savefig(output_path_obj, dpi=300, bbox_inches="tight")
    plt.close()
    
    logger.info(f"gene_annotation_plot: Saved to {output_path_obj}")
    
    return {
        "status": "success",
        "output_path": str(output_path_obj),
        "region": f"{chrom}:{start}-{end}",
        "num_variants": len(positions),
        "note": "Gene annotations require GFF3 parsing (future enhancement)",
    }


def recombination_rate_plot(
    results: list[dict[str, Any]] | Path,
    output_path: str | Path,
    *,
    recomb_map_path: Path | None = None,
    chrom: str,
    start: int,
    end: int,
    title: str | None = None,
) -> dict[str, Any]:
    """Association with recombination hotspots.
    
    Overlays recombination rates with association signal.
    Recombination hotspots can affect LD and fine-mapping.
    
    Args:
        results: Association results or path
        output_path: Output path
        recomb_map_path: Path to recombination map
        chrom: Chromosome
        start: Region start
        end: Region end
        title: Plot title
    
    Returns:
        Plot metadata
    """
    logger.info(f"recombination_rate_plot: {chrom}:{start}-{end}")
    
    return {
        "status": "skipped",
        "message": "Recombination rate plotting requires genetic map data",
        "recommendation": "Provide recombination map or use PLINK genetic map format",
        "note": "Recombination maps available from HapMap, 1000 Genomes for humans; limited for non-model organisms",
    }



