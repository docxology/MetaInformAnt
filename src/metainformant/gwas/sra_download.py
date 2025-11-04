"""Download and process sequencing data from NCBI SRA for variant calling."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Any

from ..core.io import dump_json, ensure_directory

logger = logging.getLogger(__name__)


def check_sra_tools_available() -> bool:
    """Check if SRA Toolkit is available.
    
    Returns:
        True if fastq-dump or fasterq-dump is available
    """
    for tool in ["fasterq-dump", "fastq-dump"]:
        try:
            subprocess.run(
                [tool, "--version"],
                capture_output=True,
                timeout=10,
                check=True,
            )
            logger.info(f"SRA Toolkit available: {tool}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            continue
    
    logger.warning("SRA Toolkit not found (fasterq-dump or fastq-dump)")
    return False


def download_sra_run(
    sra_accession: str,
    dest_dir: str | Path,
    *,
    use_fasterq: bool = True,
    threads: int = 4,
) -> dict[str, Any]:
    """Download SRA run and convert to FASTQ.
    
    Args:
        sra_accession: SRA run accession (e.g., SRR1234567)
        dest_dir: Destination directory for FASTQ files
        use_fasterq: Use fasterq-dump (faster) instead of fastq-dump
        threads: Number of threads for fasterq-dump
    
    Returns:
        Dictionary with download status and file paths
    """
    logger.info(f"Downloading SRA run {sra_accession}")
    
    if not check_sra_tools_available():
        return {
            "status": "failed",
            "error": "SRA Toolkit not available",
            "message": "Install SRA Toolkit: https://github.com/ncbi/sra-tools",
        }
    
    out_dir = ensure_directory(dest_dir)
    
    # Choose tool
    tool = "fasterq-dump" if use_fasterq else "fastq-dump"
    
    # Build command
    if use_fasterq:
        cmd = [
            "fasterq-dump",
            sra_accession,
            "-O", str(out_dir),
            "-e", str(threads),
            "--split-files",  # Split paired-end reads
            "-p",  # Show progress
        ]
    else:
        cmd = [
            "fastq-dump",
            sra_accession,
            "-O", str(out_dir),
            "--split-files",
            "--gzip",
        ]
    
    try:
        logger.info(f"Running: {' '.join(cmd)}")
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
            check=True,
        )
        
        # Find generated FASTQ files
        fastq_files = list(out_dir.glob(f"{sra_accession}*.fastq*"))
        
        if fastq_files:
            logger.info(f"Downloaded {len(fastq_files)} FASTQ file(s)")
            result = {
                "status": "success",
                "sra_accession": sra_accession,
                "fastq_files": [str(f) for f in fastq_files],
                "dest_dir": str(out_dir),
            }
        else:
            result = {
                "status": "failed",
                "error": "No FASTQ files generated",
            }
        
        dump_json(result, out_dir / f"{sra_accession}_download.json", indent=2)
        return result
        
    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Download timeout (>1 hour)",
            "sra_accession": sra_accession,
        }
    except subprocess.CalledProcessError as e:
        return {
            "status": "failed",
            "error": f"SRA download failed: {e}",
            "stderr": e.stderr if hasattr(e, "stderr") else "",
            "sra_accession": sra_accession,
        }


def search_sra_for_organism(
    organism: str,
    strategy: str = "WGS",
    max_results: int = 100,
) -> dict[str, Any]:
    """Search NCBI SRA for sequencing runs of a specific organism.
    
    Args:
        organism: Scientific name (e.g., "Apis mellifera")
        strategy: Sequencing strategy (WGS, RNA-Seq, etc.)
        max_results: Maximum number of results to return
    
    Returns:
        Dictionary with search results and SRA accessions
    
    Note:
        Requires esearch and efetch from NCBI E-utilities, or uses web scraping.
        For programmatic access, users should use NCBI E-utilities directly.
    """
    logger.info(f"Searching SRA for {organism} {strategy} data")
    
    # Check if esearch is available (NCBI E-utilities)
    try:
        subprocess.run(
            ["esearch", "-version"],
            capture_output=True,
            timeout=10,
            check=True,
        )
        has_eutils = True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        has_eutils = False
    
    if not has_eutils:
        return {
            "status": "pending",
            "message": "NCBI E-utilities not installed",
            "note": "Install E-utilities for programmatic SRA search",
            "manual_search_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={organism.replace(' ', '+')}+{strategy}",
        }
    
    # Use esearch to find SRA runs
    try:
        search_term = f'"{organism}"[Organism] AND "{strategy}"[Strategy]'
        cmd = [
            "esearch",
            "-db", "sra",
            "-query", search_term,
        ]
        
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        
        # Parse output for accessions (simplified)
        # Real implementation would use efetch to get full details
        
        return {
            "status": "success",
            "organism": organism,
            "strategy": strategy,
            "message": "Use NCBI SRA Run Selector for detailed results",
            "sra_search_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={organism.replace(' ', '+')}+{strategy}",
        }
        
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
        return {
            "status": "failed",
            "error": str(e),
            "manual_search_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={organism.replace(' ', '+')}+{strategy}",
        }


def download_sra_project(
    bioproject: str,
    dest_dir: str | Path,
    *,
    max_runs: int = 10,
    threads: int = 4,
) -> dict[str, Any]:
    """Download all SRA runs from a BioProject.
    
    Args:
        bioproject: NCBI BioProject accession (e.g., PRJNA292680)
        dest_dir: Destination directory
        max_runs: Maximum number of runs to download
        threads: Threads per download
    
    Returns:
        Dictionary with download status for all runs
    """
    logger.info(f"Downloading SRA project {bioproject}")
    
    out_dir = ensure_directory(dest_dir)
    
    # Note: In practice, users should use SRA Run Selector to get run list
    # This is a simplified implementation
    
    return {
        "status": "pending",
        "bioproject": bioproject,
        "message": "Use NCBI SRA Run Selector to get run list",
        "instructions": [
            f"1. Visit https://www.ncbi.nlm.nih.gov/Traces/study/?acc={bioproject}",
            "2. Select runs of interest",
            "3. Download metadata CSV or use Accession List",
            "4. Use download_sra_run() for each accession",
        ],
        "dest_dir": str(out_dir),
    }





