"""Step runner for `amalgkit quant` (quantify transcript abundances)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited, write_delimited
from ...core.logging import get_logger
from ..amalgkit import run_amalgkit

logger = get_logger(__name__)


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit quant` (quantification, e.g., Salmon)."""
    from ..amalgkit import quant as _quant
    return _quant(params, work_dir=work_dir, log_dir=log_dir, step_name="quant", check=check)


def quantify_sample(
    sample_id: str,
    metadata_rows: list[dict[str, Any]],
    quant_params: Mapping[str, Any],
    *,
    log_dir: Path | None = None,
    step_name: str | None = None,
) -> tuple[bool, str, Path | None]:
    """Quantify a single sample using amalgkit quant.
    
    Creates a temporary metadata file with just this sample and runs quantification.
    
    Args:
        sample_id: SRA accession ID (e.g., "SRR1234567")
        metadata_rows: List of metadata rows (dicts) for this sample
        quant_params: Parameters for amalgkit quant step
        log_dir: Optional directory for log files
        step_name: Optional step name for logging
        
    Returns:
        Tuple of (success: bool, message: str, abundance_file: Path | None)
        abundance_file is the path to abundance.tsv if successful, None otherwise
    """
    if not metadata_rows:
        return False, f"Sample {sample_id} not found in metadata", None
    
    # Get output directory from quant_params
    quant_dir = Path(quant_params.get("out_dir", "."))
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if already quantified
    abundance_file = quant_dir / sample_id / "abundance.tsv"
    if abundance_file.exists():
        logger.info(f"Sample {sample_id} already quantified")
        return True, f"Sample {sample_id} already quantified", abundance_file
    
    # Create temporary single-sample metadata
    temp_metadata = quant_dir.parent / f"metadata.quant.{sample_id}.tsv"
    if not temp_metadata.parent.exists():
        temp_metadata = Path.cwd() / f"metadata.quant.{sample_id}.tsv"
    
    try:
        write_delimited(metadata_rows, temp_metadata, delimiter="\t")
        
        # Prepare quant params (remove work_dir if present - it's not an amalgkit parameter)
        quant_params_copy = dict(quant_params)
        work_dir_for_quant = None
        fastq_source_dir = None
        fastq_temp_dir = None
        
        if "work_dir" in quant_params_copy:
            work_dir_path = Path(quant_params_copy.pop("work_dir"))
            # Find fastq directory - could be sibling of work_dir or in work_dir
            potential_fastq = work_dir_path.parent / "fastq"
            if potential_fastq.exists():
                fastq_source_dir = potential_fastq
            elif (work_dir_path / "fastq").exists():
                fastq_source_dir = work_dir_path / "fastq"
            
            # If work_dir is a sibling of quant_dir, use parent so amalgkit can find fastq
            # Structure: parent/{work,fastq,quant}
            if work_dir_path.parent == quant_dir.parent:
                # Use parent directory so amalgkit can find fastq relative to it
                work_dir_for_quant = quant_dir.parent
            else:
                work_dir_for_quant = work_dir_path
        elif quant_dir.parent.name == "work":
            work_dir_for_quant = quant_dir.parent
            fastq_source_dir = quant_dir.parent.parent / "fastq"
        else:
            # Try to find work_dir by looking for common parent
            # quant_dir is typically output/amalgkit/<species>/quant
            # work_dir is typically output/amalgkit/<species>/work
            # fastq_dir is typically output/amalgkit/<species>/fastq
            # All are siblings, so use parent directory
            potential_work_dir = quant_dir.parent / "work"
            if potential_work_dir.exists():
                # Use parent so amalgkit can find fastq relative to it
                work_dir_for_quant = quant_dir.parent
                fastq_source_dir = quant_dir.parent / "fastq"
            else:
                work_dir_for_quant = quant_dir.parent
                fastq_source_dir = quant_dir.parent / "fastq"
        
        # amalgkit quant looks for FASTQ files in out_dir/getfastq/<SRR>/ or out_dir/fastq/getfastq/<SRR>/
        # If fastq is a sibling, we need to create a link or copy in quant directory
        # Since symlinks don't work on ext6, create a temporary copy
        sample_fastq_source = None
        if fastq_source_dir:
            # Check both getfastq subdirectory and direct structure
            sample_fastq_source = fastq_source_dir / "getfastq" / sample_id
            if not sample_fastq_source.exists():
                sample_fastq_source = fastq_source_dir / sample_id
        
        # Create temporary copy in quant directory if needed
        if sample_fastq_source and sample_fastq_source.exists():
            fastq_temp_dir = quant_dir / "getfastq" / sample_id
            fastq_temp_dir.parent.mkdir(parents=True, exist_ok=True)
            if not fastq_temp_dir.exists() or not list(fastq_temp_dir.glob("*.fastq*")):
                import shutil
                fastq_temp_dir.mkdir(parents=True, exist_ok=True)
                # Copy only FASTQ files (not SRA files or other directories)
                copied_count = 0
                for fastq_file in sample_fastq_source.glob("*.fastq*"):
                    if fastq_file.is_file():
                        target = fastq_temp_dir / fastq_file.name
                        shutil.copy2(fastq_file, target)
                        copied_count += 1
                        logger.debug(f"Copied {fastq_file.name} to {fastq_temp_dir}")
                if copied_count > 0:
                    logger.info(f"Created temporary FASTQ copy: {fastq_temp_dir} ({copied_count} files)")
                else:
                    logger.warning(f"No FASTQ files found to copy from {sample_fastq_source}")
        
        quant_params_copy["metadata"] = str(temp_metadata.absolute())
        quant_params_copy["out_dir"] = str(quant_dir.absolute())
        
        # Run quantification
        # work_dir is passed to run_amalgkit() to set current working directory
        # amalgkit quant looks for FASTQ files relative to out_dir
        step_label = step_name or f"quant_{sample_id}"
        result = run_amalgkit(
            "quant",
            quant_params_copy,
            work_dir=str(work_dir_for_quant.absolute()) if work_dir_for_quant and work_dir_for_quant.exists() else None,
            log_dir=log_dir,
            step_name=step_label,
            check=False,
        )
        
        # Clean up temporary FASTQ copy
        if fastq_temp_dir and fastq_temp_dir.exists():
            try:
                import shutil
                shutil.rmtree(fastq_temp_dir)
                logger.debug(f"Removed temporary FASTQ copy: {fastq_temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to remove temporary FASTQ copy: {e}")
        
        # Verify quantification succeeded
        if result.returncode == 0 and abundance_file.exists():
            logger.info(f"Quantified {sample_id}")
            return True, f"Quantified {sample_id}", abundance_file
        else:
            logger.warning(f"Quantification failed for {sample_id} (code {result.returncode})")
            return False, f"Quantification failed (code {result.returncode})", None
    except Exception as e:
        logger.error(f"Error quantifying {sample_id}: {e}", exc_info=True)
        return False, str(e), None
    finally:
        # Clean up temp metadata
        try:
            if temp_metadata.exists():
                temp_metadata.unlink()
        except Exception:
            pass
