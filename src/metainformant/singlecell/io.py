"""I/O utilities for single-cell transcriptomics.

This module provides functions for fetching datasets from GEO/SRA and 
handling fastq compression/decompression.
"""
from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import List

from metainformant.core.utils import logging as core_logging

logger = core_logging.get_logger(__name__)

def download_geo_supplementary(accession: str, output_dir: str | Path) -> Path:
    """Download supplementary files for a GEO accession using wget."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Standard GEO FTP pattern
    # ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107451/suppl/
    nnn = accession[:-3] + "nnn"
    url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{accession}/suppl/"
    
    logger.info(f"Checking supplementary files for {accession}...")
    cmd = [
        "wget", "-r", "-np", "-nd", "-q", "--show-progress",
        "-P", str(output_dir),
        url
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return output_dir

def download_sra_fastqs(accession: str, output_dir: str | Path, threads: int = 1) -> List[Path]:
    """Download FASTQ files for an SRA accession using fasterq-dump and compress them."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Downloading {accession} using fasterq-dump...")
    cmd = [
        "fasterq-dump",
        "--outdir", str(output_dir),
        "--threads", str(threads),
        "--progress",
        accession
    ]
    subprocess.run(cmd, check=True)
    
    # Identify downloaded FASTQ files
    fastqs = list(output_dir.glob(f"{accession}*.fastq"))
    logger.info(f"Downloaded {len(fastqs)} FASTQ files for {accession}.")
    
    # Compress files
    for fq in fastqs:
        _compress_fastq(fq)
        
    return list(output_dir.glob(f"{accession}*.fastq.gz"))

def _compress_fastq(path: Path) -> None:
    """Compress a FASTQ file using pigz or gzip."""
    # Check for pigz availability
    try:
        subprocess.run(["pigz", "--version"], capture_output=True, check=True)
        compressor = "pigz"
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.warning("pigz not found, using standard gzip...")
        compressor = "gzip"
        
    logger.info(f"Compressing {path.name} with {compressor}...")
    
    # USE -f TO PREVENT STALLING ON OVERWRITE
    cmd = [compressor, "-f", str(path)]
    subprocess.run(cmd, check=True)

def fetch_atlas_datasets(datasets: List[str], output_base: str | Path) -> List[Path]:
    """Fetch multiple datasets (Davie-2018, Baker-2021, etc.) to the output base."""
    output_base = Path(output_base)
    output_base.mkdir(parents=True, exist_ok=True)
    fetched_paths = []
    
    dataset_map = {
        "Davie-2018": "GSE107451",
        "Baker-2021": "GSE152495",
        "Park-2022": "GSE207799",
        "Dopp-2024": "GSE221239",
        "LeeBenton-2023": "GSE247965",
        "Li-2022": "E-MTAB-10519",
    }
    
    for ds_name in datasets:
        acc = dataset_map.get(ds_name)
        if not acc:
            logger.warning(f"Unknown dataset name: {ds_name}")
            continue
            
        ds_dir = output_base / ds_name
        logger.info(f"Fetching {ds_name} ({acc}) -> {ds_dir}")
        
        if acc.startswith("GSE"):
            download_geo_supplementary(acc, ds_dir)
        elif acc.startswith("E-"):
            logger.info(f"Skipping direct fetch for {acc} (implement custom handler).")
            
        fetched_paths.append(ds_dir)
        
    return fetched_paths

def run_salmon_alevin(fastq_r1: Path, fastq_r2: Path, index_dir: Path, tg_map: Path, output_dir: Path, threads: int = 1) -> None:
    """Run Salmon Alevin for single-cell quantification."""
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running Salmon Alevin quantification into {output_dir}...")
    
    cmd = [
        "salmon", "alevin",
        "-l", "ISR",
        "-1", str(fastq_r1),
        "-2", str(fastq_r2),
        "--chromiumV3",  # Assuming V3, can be parameterized
        "-i", str(index_dir),
        "-p", str(threads),
        "-o", str(output_dir),
        "--tgMap", str(tg_map)
    ]
    logger.info(f"Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_seurat_integration(input_dirs: dict, output_dir: Path, resolution: float = 0.8) -> Path:
    """Stub for Seurat integration (usually calls an R script)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running Seurat integration (resolution={resolution}) into {output_dir}...")
    # In a real scenario, this would invoke Rscript
    rds_path = output_dir / "integrated.rds"
    # Touch file for now to allow pipeline to "progress" or fail informatively
    rds_path.touch()
    return rds_path

def run_differential_expression(input_rds: Path, output_csv: Path, contrast: str) -> Path:
    """Stub for differential expression."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running Differential Expression with contrast {contrast}...")
    output_csv.touch()
    return output_csv

def run_scenic_grnboost2(expression_loom: Path, tf_txt: Path, output_adjacencies: Path, threads: int = 1) -> Path:
    """Stub for pySCENIC GRNBoost2."""
    output_adjacencies.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running pySCENIC GRNBoost2...")
    output_adjacencies.touch()
    return output_adjacencies

def run_scenic_cistarget(expression_loom: Path, adjacencies: Path, motif_database: Path, motif_annotations: Path, output_regulons: Path) -> Path:
    """Stub for pySCENIC cisTarget."""
    output_regulons.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running pySCENIC cisTarget...")
    output_regulons.touch()
    return output_regulons

def run_scenic_aucell(expression_loom: Path, regulons: Path, output_loom: Path) -> Path:
    """Stub for pySCENIC AUCell."""
    output_loom.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running pySCENIC AUCell...")
    output_loom.touch()
    return output_loom

def filter_regulons_post_scenic(regulons_csv: Path, filtered_output: Path, raw_expression_csv: Path) -> Path:
    """Stub for filtering regulons."""
    filtered_output.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Filtering regulons...")
    filtered_output.touch()
    return filtered_output

def run_pseudotime_trajectory(input_rds: Path, output_rds: Path, root_selector: str) -> Path:
    """Stub for Monocle3 trajectory."""
    output_rds.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running pseudotime trajectory...")
    output_rds.touch()
    return output_rds

def run_go_enrichment(input_deseq_csv: Path, output_dir: Path, organism: str) -> Path:
    """Stub for clusterProfiler GO enrichment."""
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running GO enrichment for {organism}...")
    return output_dir

def run_shinycell_export(input_rds: Path, output_dir: Path) -> Path:
    """Stub for ShinyCell export."""
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Exporting to ShinyCell...")
    return output_dir
