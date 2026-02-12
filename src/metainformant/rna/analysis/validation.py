"""RNA-seq workflow sample validation utilities.

This module provides comprehensive validation functions to track samples
through the complete pipeline: download → extract → quant → merge.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging
from metainformant.core import io
from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig

logger = logging.get_logger(__name__)


def get_sample_pipeline_status(
    sample_id: str, work_dir: Path, fastq_dir: Optional[Path] = None, quant_dir: Optional[Path] = None
) -> Dict[str, Any]:
    """Get current pipeline status for a single sample.

    Args:
        sample_id: Sample identifier (SRA run ID, e.g., SRR123456)
        work_dir: Workflow working directory
        fastq_dir: FASTQ directory (default: work_dir/fastq or work_dir/getfastq)
        quant_dir: Quantification directory (default: work_dir/quant)

    Returns:
        Dictionary with status at each pipeline stage:
        - download: bool (SRA file exists)
        - extraction: bool (FASTQ files exist)
        - quantification: bool (abundance file exists)
        - merge: bool (sample in merged matrix)
        - stage: str (current stage: 'download', 'extraction', 'quantification', 'merge', 'complete')
        - diagnostics: dict (detailed file paths and sizes)
    """
    # Resolve fastq_dir - handle both direct and amalgkit getfastq/ subdirectory structures
    if fastq_dir is None:
        # Check both possible locations in work_dir
        fastq_dir = work_dir / "fastq"
        getfastq_dir = work_dir / "getfastq"
        if getfastq_dir.exists():
            fastq_dir = getfastq_dir
    else:
        fastq_dir = Path(fastq_dir)
        # When fastq_dir is provided, check if it has a getfastq/ subdirectory (amalgkit structure)
        # Prefer amalgkit structure if it exists
        # But don't check if fastq_dir is already a getfastq directory (avoid getfastq/getfastq)
        if fastq_dir.name != "getfastq":
            getfastq_subdir = fastq_dir / "getfastq"
            if getfastq_subdir.exists() and getfastq_subdir.is_dir():
                fastq_dir = getfastq_subdir

    if quant_dir is None:
        quant_dir = work_dir / "quant"
    else:
        quant_dir = Path(quant_dir)

    status = {
        "sample_id": sample_id,
        "download": False,
        "extraction": False,
        "quantification": False,
        "merge": False,
        "stage": "not_started",
        "diagnostics": {},
    }

    # Check download stage: SRA file exists
    # Try both possible locations: fastq_dir/sample_id/ and fastq_dir.parent/sample_id/
    # (in case fastq_dir is already the getfastq subdirectory)
    sra_file = fastq_dir / sample_id / f"{sample_id}.sra"
    if not sra_file.exists() and fastq_dir.name == "getfastq":
        # Also check parent directory structure
        alt_sra_file = fastq_dir.parent / sample_id / f"{sample_id}.sra"
        if alt_sra_file.exists():
            sra_file = alt_sra_file

    if sra_file.exists():
        status["download"] = True
        status["diagnostics"]["sra_file"] = str(sra_file)
        status["diagnostics"]["sra_size"] = sra_file.stat().st_size
        status["stage"] = "download"

    # Check extraction stage: FASTQ files exist
    # Try both possible locations
    sample_fastq_dir = fastq_dir / sample_id
    if not sample_fastq_dir.exists() and fastq_dir.name == "getfastq":
        # Also check parent directory structure
        alt_sample_dir = fastq_dir.parent / sample_id
        if alt_sample_dir.exists():
            sample_fastq_dir = alt_sample_dir

    if sample_fastq_dir.exists():
        fastq_files = list(sample_fastq_dir.glob("*.fastq*"))
        if fastq_files:
            status["extraction"] = True
            status["diagnostics"]["fastq_files"] = [str(f) for f in fastq_files]
            status["diagnostics"]["fastq_count"] = len(fastq_files)
            status["diagnostics"]["fastq_total_size"] = sum(f.stat().st_size for f in fastq_files)
            status["stage"] = "extraction"

    # Check quantification stage: abundance file exists
    # Amalgkit creates quant/{sample_id}/abundance.tsv (kallisto) or quant/{sample_id}/quant.sf (salmon)
    sample_quant_dir = quant_dir / sample_id
    if sample_quant_dir.exists():
        abundance_tsv = sample_quant_dir / "abundance.tsv"
        sample_abundance_tsv = sample_quant_dir / f"{sample_id}_abundance.tsv"
        quant_sf = sample_quant_dir / "quant.sf"

        if (abundance_tsv.exists() and abundance_tsv.stat().st_size > 0) or (
            sample_abundance_tsv.exists() and sample_abundance_tsv.stat().st_size > 0
        ):
            target_abundance = abundance_tsv if abundance_tsv.exists() else sample_abundance_tsv
            status["quantification"] = True
            status["diagnostics"]["abundance_file"] = str(target_abundance)
            status["diagnostics"]["abundance_size"] = target_abundance.stat().st_size
            status["stage"] = "quantification"
        elif quant_sf.exists() and quant_sf.stat().st_size > 0:
            status["quantification"] = True
            status["diagnostics"]["abundance_file"] = str(quant_sf)
            status["diagnostics"]["abundance_size"] = quant_sf.stat().st_size
            status["stage"] = "quantification"

    # Check merge stage: sample appears in merged abundance matrix
    # Optimized: stream file reading with early exit
    merge_dir = work_dir / "merge"
    merged_abundance = merge_dir / "merged_abundance.tsv"
    if merged_abundance.exists():
        try:
            # Stream file reading for large files - check header first, then sample rows
            with open(merged_abundance, "r") as f:
                header = f.readline()
                if sample_id in header:
                    status["merge"] = True
                    status["stage"] = "merge"
                else:
                    # Check first 1000 data lines for sample ID (early exit optimization)
                    # Most samples appear in first few rows, avoid reading entire large file
                    max_lines_to_check = 1000
                    lines_checked = 0
                    for line in f:
                        lines_checked += 1
                        if line.startswith(sample_id + "\t") or line.startswith(sample_id + ","):
                            status["merge"] = True
                            status["stage"] = "merge"
                            break
                        if lines_checked >= max_lines_to_check:
                            # For very large files, stop after checking reasonable number of lines
                            # If sample not found in first 1000 rows, likely not present
                            break
        except Exception as e:
            logger.debug(
                f"Could not check merge status for {sample_id} in {merged_abundance}: {e}\n"
                f"Remediation: Ensure merged abundance file exists and is readable."
            )

    # Determine overall stage
    if status["merge"]:
        status["stage"] = "complete"
    elif status["quantification"]:
        status["stage"] = "quantification"
    elif status["extraction"]:
        status["stage"] = "extraction"
    elif status["download"]:
        status["stage"] = "download"

    return status


def validate_sample_end_to_end(
    sample_id: str, work_dir: Path, fastq_dir: Optional[Path] = None, quant_dir: Optional[Path] = None
) -> Dict[str, Any]:
    """Validate that a sample completed all pipeline stages end-to-end.

    Args:
        sample_id: Sample identifier
        work_dir: Workflow working directory
        fastq_dir: FASTQ directory (default: inferred from work_dir)
        quant_dir: Quantification directory (default: inferred from work_dir)

    Returns:
        Validation result dictionary:
        - valid: bool (all stages complete)
        - sample_id: str
        - stages: dict (status for each stage)
        - missing_stages: list (stages that are incomplete)
        - diagnostics: dict (detailed file information)
    """
    status = get_sample_pipeline_status(sample_id, work_dir, fastq_dir, quant_dir)

    result = {
        "valid": False,
        "sample_id": sample_id,
        "stages": {
            "download": status["download"],
            "extraction": status["extraction"],
            "quantification": status["quantification"],
            "merge": status["merge"],
        },
        "missing_stages": [],
        "current_stage": status["stage"],
        "diagnostics": status["diagnostics"],
    }

    # Check which stages are missing with diagnostic information
    # Handle None values for directories in error messages
    fastq_path_str = str(fastq_dir / sample_id) if fastq_dir else "fastq_dir/sample_id"
    quant_path_str = str(quant_dir / sample_id) if quant_dir else "quant_dir/sample_id"

    if not status["download"] and not status["extraction"]:
        result["missing_stages"].append("download")
        result["diagnostics"]["download_issue"] = (
            f"No SRA or FASTQ files found for {sample_id}.\n"
            f"Expected locations:\n"
            f"  - {fastq_path_str}/{sample_id}.sra\n"
            f"  - {fastq_path_str}/*.fastq*\n"
            f"Remediation: Run 'getfastq' step to download and extract sample."
        )
    if not status["extraction"]:
        result["missing_stages"].append("extraction")
        if "extraction_issue" not in result["diagnostics"]:
            result["diagnostics"]["extraction_issue"] = (
                f"No FASTQ files found for {sample_id}.\n"
                f"Expected location: {fastq_path_str}/*.fastq*\n"
                f"Remediation: Ensure 'getfastq' step completed successfully."
            )
    if not status["quantification"]:
        result["missing_stages"].append("quantification")
        result["diagnostics"]["quantification_issue"] = (
            f"No abundance file found for {sample_id}.\n"
            f"Expected locations:\n"
            f"  - {quant_path_str}/abundance.tsv (kallisto)\n"
            f"  - {quant_path_str}/quant.sf (salmon)\n"
            f"Remediation: Run 'quant' step to quantify sample."
        )
    if not status["merge"]:
        result["missing_stages"].append("merge")
        # Don't add diagnostic for merge as it's optional

    # Sample is valid if it has at least extraction and quantification
    # (merge is optional for per-sample validation)
    result["valid"] = status["extraction"] and status["quantification"]

    return result


def validate_all_samples(config: AmalgkitWorkflowConfig, stage: Optional[str] = None) -> Dict[str, Any]:
    """Validate all samples from metadata through the pipeline.

    Args:
        config: Workflow configuration
        stage: Specific stage to validate ('download', 'extraction', 'quantification', 'merge')
               If None, validates all stages

    Returns:
        Validation summary:
        - total_samples: int
        - validated: int (samples that passed validation)
        - failed: int (samples that failed)
        - missing_stages: dict (count of samples missing each stage)
        - per_sample: dict (detailed status for each sample)
        - summary: dict (stage-specific summaries)
    """
    work_dir = config.work_dir

    # Load metadata to get list of samples
    metadata_file = work_dir / "metadata" / "metadata_selected.tsv"
    if not metadata_file.exists():
        metadata_file = work_dir / "metadata" / "metadata.tsv"

    if not metadata_file.exists():
        error_msg = (
            f"Metadata file not found: {metadata_file}\n"
            f"Expected locations:\n"
            f"  - {work_dir / 'metadata' / 'metadata_selected.tsv'}\n"
            f"  - {work_dir / 'metadata' / 'metadata.tsv'}\n"
            f"Remediation: Run 'metadata' and 'select' steps to generate metadata files."
        )
        logger.warning(error_msg)
        return {
            "total_samples": 0,
            "validated": 0,
            "failed": 0,
            "missing_stages": {},
            "per_sample": {},
            "summary": {},
            "error": "Metadata file not found",
            "error_details": error_msg,
        }

    # Read sample IDs from metadata
    # Optimized: case-insensitive column matching, more robust detection
    sample_ids = []
    try:
        import csv

        with open(metadata_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            # Get column names (case-insensitive matching)
            fieldnames_lower = {name.lower(): name for name in reader.fieldnames or []}

            # Try various column name variations (case-insensitive)
            sample_id_key = None
            for key_variant in ["run", "sra_run", "sample_id", "accession"]:
                if key_variant.lower() in fieldnames_lower:
                    sample_id_key = fieldnames_lower[key_variant.lower()]
                    break

            # Fallback to original case-sensitive matching if needed
            if sample_id_key is None:
                sample_id_key = "run"  # Most common default

            for row in reader:
                # Get sample ID with fallback to multiple column names
                sample_id = (
                    row.get(sample_id_key, "")
                    or row.get("run", "")
                    or row.get("Run", "")
                    or row.get("SRA_Run", "")
                    or row.get("sample_id", "")
                    or row.get("accession", "")
                )
                if sample_id and sample_id.strip():  # Skip empty/whitespace-only IDs
                    sample_ids.append(sample_id.strip())
    except Exception as e:
        error_msg = (
            f"Failed to read metadata file: {metadata_file}\n"
            f"Error: {e}\n"
            f"Remediation: Check that metadata file is valid TSV format with 'run' column."
        )
        logger.error(error_msg)
        return {
            "total_samples": 0,
            "validated": 0,
            "failed": 0,
            "missing_stages": {},
            "per_sample": {},
            "summary": {"error": str(e)},
            "error_details": error_msg,
        }

    if not sample_ids:
        error_msg = (
            f"No samples found in metadata file: {metadata_file}\n"
            f"Remediation: Ensure metadata file contains sample IDs in 'run', 'Run', 'SRA_Run', "
            f"'sample_id', or 'accession' column."
        )
        logger.warning(error_msg)
        return {
            "total_samples": 0,
            "validated": 0,
            "failed": 0,
            "missing_stages": {},
            "per_sample": {},
            "summary": {},
            "error_details": error_msg,
        }

    # Get step-specific directories from config
    fastq_dir = None
    quant_dir = None

    steps = config.extra_config.get("steps", {})
    if "getfastq" in steps:
        fastq_dir_raw = steps["getfastq"].get("out_dir", work_dir / "fastq")
        fastq_dir = Path(fastq_dir_raw)
        # Check if amalgkit getfastq/ subdirectory exists (prefer it if present)
        # Also check if fastq_dir itself is already the getfastq directory
        if fastq_dir.name != "getfastq":
            getfastq_subdir = fastq_dir / "getfastq"
            if getfastq_subdir.exists():
                fastq_dir = getfastq_subdir
        # If fastq_dir is already getfastq, use it as-is
    if "quant" in steps:
        quant_dir = Path(steps["quant"].get("out_dir", work_dir / "quant"))

    # Validate each sample
    per_sample = {}
    missing_stages_count = {"download": 0, "extraction": 0, "quantification": 0, "merge": 0}

    validated_count = 0
    failed_count = 0

    # Progress indication for large sample sets
    total_samples = len(sample_ids)
    log_progress_interval = max(1, total_samples // 10)  # Log every 10% or at least every sample

    for idx, sample_id in enumerate(sample_ids, start=1):
        # Log progress for large sample sets
        if total_samples > 10 and idx % log_progress_interval == 0:
            logger.debug(f"Validation progress: {idx}/{total_samples} samples ({100*idx//total_samples}%)")

        if stage:
            # Validate specific stage only
            status = get_sample_pipeline_status(sample_id, work_dir, fastq_dir, quant_dir)
            stage_status = {
                "download": status["download"],
                "extraction": status["extraction"],
                "quantification": status["quantification"],
                "merge": status["merge"],
            }

            if stage == "download":
                # For download, check if SRA file exists
                is_valid = status["download"]
            elif stage == "extraction":
                # For extraction, check if FASTQ files exist
                is_valid = status["extraction"]
            elif stage == "quantification":
                is_valid = status["quantification"]
            elif stage == "merge":
                is_valid = status["merge"]
            else:
                is_valid = False

            per_sample[sample_id] = {"valid": is_valid, "stages": stage_status, "current_stage": status["stage"]}

            if is_valid:
                validated_count += 1
            else:
                failed_count += 1
                missing_stages_count[stage] += 1
        else:
            # Validate end-to-end
            validation = validate_sample_end_to_end(sample_id, work_dir, fastq_dir, quant_dir)
            per_sample[sample_id] = validation

            if validation["valid"]:
                validated_count += 1
            else:
                failed_count += 1
                for missing in validation["missing_stages"]:
                    missing_stages_count[missing] += 1

    # Generate summary by stage
    stage_summaries = {
        "download": {
            "total": len(sample_ids),
            "complete": sum(1 for s in per_sample.values() if s.get("stages", {}).get("download", False)),
            "missing": missing_stages_count["download"],
        },
        "extraction": {
            "total": len(sample_ids),
            "complete": sum(1 for s in per_sample.values() if s.get("stages", {}).get("extraction", False)),
            "missing": missing_stages_count["extraction"],
        },
        "quantification": {
            "total": len(sample_ids),
            "complete": sum(1 for s in per_sample.values() if s.get("stages", {}).get("quantification", False)),
            "missing": missing_stages_count["quantification"],
        },
        "merge": {
            "total": len(sample_ids),
            "complete": sum(1 for s in per_sample.values() if s.get("stages", {}).get("merge", False)),
            "missing": missing_stages_count["merge"],
        },
    }

    return {
        "total_samples": len(sample_ids),
        "validated": validated_count,
        "failed": failed_count,
        "missing_stages": missing_stages_count,
        "per_sample": per_sample,
        "summary": stage_summaries,
        "validation_stage": stage or "all",
    }


def save_validation_report(validation_result: Dict[str, Any], output_path: Path) -> None:
    """Save validation report to JSON file.

    Args:
        validation_result: Validation result from validate_all_samples()
        output_path: Path to save JSON report
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io.dump_json(validation_result, output_path)
    logger.info(f"Saved validation report to {output_path}")
