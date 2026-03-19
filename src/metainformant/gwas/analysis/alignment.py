"""Read alignment and BAM processing for GWAS analysis.

Provides wrappers for BWA-MEM and samtools to align sequencing reads
to a reference genome and process the resulting alignments for downstream
variant calling.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def align_reads_bwa(
    reference_fasta: str | Path,
    fastq_files: List[str | Path],
    output_sam: str | Path,
    threads: int = 4,
    sample_id: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """Align sequencing reads to a reference using BWA-MEM.

    Args:
        reference_fasta: Reference genome FASTA.
        fastq_files: List of FASTQ files (1 for single-end, 2 for paired-end).
        output_sam: Output SAM file path.
        threads: Number of threads to use.
        sample_id: Optional ID for the Read Group (@RG).

    Returns:
        CompletedProcess with execution results.
    """
    output_sam = Path(output_sam)
    output_sam.parent.mkdir(parents=True, exist_ok=True)

    # Check if BWA is available
    try:
        subprocess.run(["bwa"], capture_output=True, check=False)
    except FileNotFoundError:
        raise FileNotFoundError("bwa not found on PATH")

    # Build BWA-MEM command
    cmd = [
        "bwa",
        "mem",
        "-t",
        str(threads),
    ]

    if sample_id:
        cmd.extend(["-R", f"@RG\\tID:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tSM:{sample_id}"])

    cmd.append(str(reference_fasta))
    cmd.extend([str(f) for f in fastq_files])

    logger.info(f"Running BWA-MEM alignment: {' '.join(cmd[:5])}...")

    try:
        with open(output_sam, "w") as sam_out:
            result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True)
        result.check_returncode()
        
        logger.info(f"Alignment completed: {output_sam}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"BWA-MEM alignment failed: {e.stderr}")
        raise


def sort_index_bam(
    input_sam_or_bam: str | Path, output_bam: str | Path, threads: int = 4
) -> Path:
    """Sort and index an alignment file using samtools.

    Args:
        input_sam_or_bam: Input SAM or uncompressed BAM.
        output_bam: Path for the final sorted BAM.
        threads: Number of threads for sorting.

    Returns:
        Path to the final indexed BAM.
    """
    input_path = Path(input_sam_or_bam)
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. View and Sort
        sort_cmd = [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-o",
            str(output_bam),
            str(input_path),
        ]
        logger.info(f"Sorting alignments: {output_bam}")
        subprocess.run(sort_cmd, check=True, capture_output=True)

        # 2. Index
        index_cmd = ["samtools", "index", str(output_bam)]
        logger.info(f"Indexing BAM: {output_bam}")
        subprocess.run(index_cmd, check=True, capture_output=True)

        return output_bam
    except subprocess.CalledProcessError as e:
        logger.error(f"Samtools processing failed: {e.stderr.decode()}")
        raise
    except FileNotFoundError:
        raise FileNotFoundError("samtools not found on PATH")


def check_alignment_tools() -> bool:
    """Check if required alignment tools (bwa, samtools) are available.

    Returns:
        True if all tools are found.
    """
    tools = ["bwa", "samtools"]
    for tool in tools:
        try:
            subprocess.run([tool], capture_output=True, check=False)
        except FileNotFoundError:
            logger.error(f"Required tool not found: {tool}")
            return False
    return True
