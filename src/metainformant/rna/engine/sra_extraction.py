"""SRA file extraction and fallback recovery utilities.

Extracted from workflow.py to reduce its complexity. These functions handle
direct SRA file extraction using fasterq-dump when amalgkit's built-in
extraction fails or skips files.
"""

from __future__ import annotations

import concurrent.futures
import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List

from metainformant.core import logging

if TYPE_CHECKING:
    from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig

logger = logging.get_logger(__name__)


def extract_sra_directly(config: AmalgkitWorkflowConfig, sra_dir: Path, output_dir: Path) -> int:
    """Manually extract SRA files using fasterq-dump when amalgkit fails.

    This acts as a fallback when amalgkit getfastq skips extraction due to
    existing SRA files (redo: no) or other internal checks.

    Args:
        config: Workflow configuration
        sra_dir: Directory containing .sra files
        output_dir: Directory to output .fastq.gz files

    Returns:
        Number of successfully extracted samples
    """
    from metainformant.core.io.download import monitor_subprocess_directory_growth

    # Check if fasterq-dump is available
    fasterq_dump = shutil.which("fasterq-dump")
    if not fasterq_dump:
        logger.error("fasterq-dump not found in PATH - cannot run fallback extraction")
        return 0

    # Check for gzip/pigz
    gzip_cmd = shutil.which("pigz") or shutil.which("gzip")
    if not gzip_cmd:
        logger.error("gzip/pigz not found - cannot compress FASTQ output")
        return 0

    # Find SRA files (root, sra/, or sample subdirs)
    sra_files = list(sra_dir.glob("*.sra"))
    sra_files.extend(list(sra_dir.glob("sra/*.sra")))
    sra_files.extend(list(sra_dir.glob("*/*.sra")))
    # Deduplicate by path
    sra_files = list({str(f): f for f in sra_files}.values())
    if not sra_files:
        logger.warning(f"No SRA files found in {sra_dir} for fallback extraction")
        return 0

    logger.info(f"Attempting fallback extraction for {len(sra_files)} files using {fasterq_dump}...")

    # Configure path for fasterq-dump temp files
    repo_root = Path(__file__).resolve().parent.parent.parent.parent.parent.parent
    tmp_dir = repo_root / ".tmp" / "fasterq-dump"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    def process_sra(sra_file: Path) -> bool:
        sample_id = sra_file.stem
        sample_out_dir = output_dir / sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)

        # Check if already extracted
        if list(sample_out_dir.glob("*.fastq.gz")):
            return True

        try:
            cmd = [
                fasterq_dump,
                "--split-3",
                "--threads",
                str(min(config.threads, 4)),
                "--outdir",
                str(sample_out_dir),
                "--temp",
                str(tmp_dir),
                str(sra_file),
            ]

            logger.info(f"Extracting {sample_id} with fasterq-dump...")

            heartbeat_file = config.work_dir / "heartbeat" / f"fallback_{sample_id}.json"
            heartbeat_file.parent.mkdir(exist_ok=True, parents=True)

            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

            rc, _ = monitor_subprocess_directory_growth(
                process=process,
                watch_dir=sample_out_dir,
                heartbeat_path=heartbeat_file,
                desc=f"Extracting {sample_id}",
                heartbeat_interval=10,
            )

            if rc != 0:
                raise RuntimeError(f"fasterq-dump exited with {rc}")

            # Gzip output files
            fastqs = list(sample_out_dir.glob("*.fastq"))
            if not fastqs:
                logger.warning(f"No FASTQ exported for {sample_id}")
                return False

            for fq in fastqs:
                subprocess.run([gzip_cmd, "-f", str(fq)], check=True)

            return True

        except Exception as e:
            logger.error(f"Failed to extract {sample_id}: {e}")
            return False

    # Run in parallel
    max_workers = max(1, config.threads // 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_sra, sra_files))

    extracted_count = sum(results)
    logger.info(f"Fallback extraction completed: {extracted_count}/{len(sra_files)} samples extracted")

    # Clean up temp dir
    try:
        shutil.rmtree(tmp_dir)
    except (OSError, PermissionError, FileNotFoundError):
        pass

    return extracted_count


def manual_integration_fallback(config: AmalgkitWorkflowConfig) -> bool:
    """Fallback for when 'amalgkit integrate' fails (e.g. path detection bugs).

    Manually exposes FASTQ files from getfastq subdirectories to the location
    expected by 'amalgkit quant'.

    Args:
        config: Workflow configuration with work_dir and step settings

    Returns:
        True if any FASTQ files were successfully linked/copied.
    """
    import re

    logger.info("Attempting Manual Integration Fallback...")

    steps_config = config.extra_config.get("steps", {})
    getfastq_dir_raw = steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq" / "getfastq")
    getfastq_dir = Path(getfastq_dir_raw)

    if not getfastq_dir.exists():
        getfastq_dir = config.work_dir / "fastq" / "getfastq"

    if not getfastq_dir.exists():
        logger.error(f"Cannot perform manual integration: getfastq dir {getfastq_dir} not found")
        return False

    quant_out_dir = Path(steps_config.get("quant", {}).get("out_dir", config.work_dir))
    quant_input_dir = quant_out_dir / "getfastq"
    quant_input_dir.mkdir(parents=True, exist_ok=True)

    found_any = False

    for fastq_file in getfastq_dir.glob("**/*.fastq*"):
        if fastq_file.is_file():
            sample_id = None
            if fastq_file.parent.name.startswith(("SRR", "ERR", "DRR")):
                sample_id = fastq_file.parent.name
            else:
                match = re.search(r"(SRR\d+|ERR\d+|DRR\d+)", fastq_file.name)
                if match:
                    sample_id = match.group(1)

            if sample_id:
                sample_dest_dir = quant_input_dir / sample_id
                sample_dest_dir.mkdir(parents=True, exist_ok=True)
                dest_path = sample_dest_dir / fastq_file.name

                if not dest_path.exists():
                    try:
                        dest_path.symlink_to(fastq_file.resolve())
                        found_any = True
                        logger.info(f"Linked {fastq_file.name} -> {sample_id}")
                    except Exception:
                        try:
                            shutil.copy2(fastq_file, dest_path)
                            found_any = True
                        except Exception as e:
                            logger.warning(f"Failed to link {fastq_file}: {e}")
                else:
                    found_any = True

    if found_any:
        metadata_tsv = config.work_dir / "metadata" / "metadata.tsv"
        selected_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
        if not metadata_tsv.exists() and selected_metadata.exists():
            shutil.copy2(selected_metadata, metadata_tsv)
            logger.info("Created metadata.tsv from metadata_selected.tsv for quant step")

        logger.info("Manual integration completed: valid FASTQ files exposed for quant.")
        return True

    return False
