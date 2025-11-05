"""Sample processing pipeline: SRA→FASTQ→Quant→Delete.

This module provides orchestration functions for processing individual samples
through the complete pipeline: convert SRA to FASTQ (if needed), quantify,
and delete FASTQ files to free disk space.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited
from .getfastq import convert_sra_to_fastq, delete_sample_fastqs
from .quant import quantify_sample

logger = logging.getLogger(__name__)


def process_sample_pipeline(
    sample_id: str,
    config_path: Path,
    status: str,
    *,
    log_dir: Path | None = None,
) -> tuple[bool, str]:
    """Process a complete sample: SRA→FASTQ (if needed) → Quant → Delete.
    
    This is the main orchestration function that handles the complete pipeline
    for a single sample. It determines what stage the sample is at and processes
    it accordingly.
    
    Args:
        sample_id: SRA accession ID (e.g., "SRR1234567")
        config_path: Path to species workflow config file
        status: Current status - "sra" (needs SRA→FASTQ conversion) or "fastq" (ready for quant)
        log_dir: Optional directory for log files
        
    Returns:
        Tuple of (success: bool, message: str)
    """
    try:
        from ..workflow import load_workflow_config
        cfg = load_workflow_config(config_path)
        
        # Get directories from config
        getfastq_params = dict(cfg.per_step.get("getfastq", {}))
        fastq_dir = Path(getfastq_params.get("out_dir", cfg.work_dir / "fastq"))
        quant_params = dict(cfg.per_step.get("quant", {}))
        quant_dir = Path(quant_params.get("out_dir", cfg.work_dir / "quant"))
        
        # Get metadata file
        metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"
        
        if not metadata_file.exists():
            return False, f"No metadata file found"
        
        # Read metadata to find this sample
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        sample_rows = [row for row in rows if row.get("run") == sample_id]
        
        if not sample_rows:
            return False, f"Sample {sample_id} not found in metadata"
        
        # Step 1: Convert SRA to FASTQ if needed
        if status == "sra":
            # Find SRA file
            sample_dir = fastq_dir / "getfastq" / sample_id
            if not sample_dir.exists():
                sample_dir = fastq_dir / sample_id
            
            if not sample_dir.exists():
                return False, f"Sample directory not found: {sample_dir}"
            
            sra_files = list(sample_dir.glob("*.sra"))
            if not sra_files:
                return False, f"No SRA file found for {sample_id}"
            
            sra_file = sra_files[0]
            threads = cfg.threads or getfastq_params.get("threads", 4)
            
            logger.info(f"Converting SRA to FASTQ for {sample_id}...")
            success, message, fastq_files = convert_sra_to_fastq(
                sample_id,
                sra_file,
                sample_dir,
                threads=threads,
                log_dir=log_dir,
            )
            
            if not success:
                return False, f"SRA conversion failed: {message}"
            
            # Verify FASTQ files were created
            if not fastq_files:
                return False, "SRA conversion completed but no FASTQ files found"
        
        # Step 2: Quantify
        logger.info(f"Quantifying {sample_id}...")
        quant_success, quant_message, abundance_file = quantify_sample(
            sample_id,
            sample_rows,
            quant_params,
            log_dir=log_dir,
            step_name=f"quant_{sample_id}",
        )
        
        if not quant_success:
            # Still delete FASTQ to free space even if quant failed
            logger.warning(f"Quantification failed, but deleting FASTQ files to free space")
            delete_sample_fastqs(sample_id, fastq_dir)
            return False, f"Quantification failed: {quant_message}"
        
        # Step 3: Delete FASTQ files
        logger.info(f"Deleting FASTQ files for {sample_id}...")
        delete_sample_fastqs(sample_id, fastq_dir)
        
        return True, f"Completed pipeline for {sample_id}: converted, quantified, and deleted"
    
    except Exception as e:
        logger.error(f"Error processing sample {sample_id}: {e}", exc_info=True)
        return False, str(e)

