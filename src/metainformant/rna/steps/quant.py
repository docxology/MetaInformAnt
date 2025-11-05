"""Step runner for `amalgkit quant` (quantify transcript abundances)."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited, write_delimited
from ..amalgkit import run_amalgkit

logger = logging.getLogger(__name__)


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
        
        # Prepare quant params
        quant_params_copy = dict(quant_params)
        quant_params_copy["metadata"] = str(temp_metadata.absolute())
        quant_params_copy["out_dir"] = str(quant_dir.absolute())
        
        # Run quantification
        step_label = step_name or f"quant_{sample_id}"
        result = run_amalgkit(
            "quant",
            quant_params_copy,
            work_dir=None,
            log_dir=log_dir,
            step_name=step_label,
            check=False,
        )
        
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
