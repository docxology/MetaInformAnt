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

    This function provides a Python wrapper around the BWA-MEM aligner,
    handling read alignment for both single-end and paired-end sequencing
    data. It is typically used as part of a GWAS pipeline for variant
    calling from short-read sequencing data.

    Args:
        reference_fasta: Path to reference genome FASTA file. Must be indexed
            (i.e., reference.fasta.fai must exist). If not indexed, BWA will
            create indices automatically on first use.
        fastq_files: List of FASTQ files. For single-end data, provide one
            file path. For paired-end data, provide two files in order:
            [forward_reads.fastq, reverse_reads.fastq]. Gzip-compressed
            files (.fastq.gz) are supported as BWA uses stdin/stdout
            processing which integrates with the shell's pipe handling.
        output_sam: Output SAM file path. Will be created with parent
            directories if they don't exist.
        threads: Number of CPU threads for parallel alignment. Default is 4.
            Higher thread counts improve speed but consume more memory.
        sample_id: Optional sample identifier for Read Group (@RG) tag in
            SAM header. This is critical for downstream tools like GATK
            that require proper read grouping for duplicate marking and
            base quality score recalibration.

    Returns:
        subprocess.CompletedProcess: Object containing the return code,
        stdout, and stderr from the BWA-MEM command. A return code of 0
        indicates successful alignment.

    Raises:
        FileNotFoundError: If 'bwa' is not found in PATH or reference file
            doesn't exist.
        subprocess.CalledProcessError: If BWA-MEM alignment fails (non-zero
            return code). The error message contains stderr output.

    Example:
        >>> # Single-end alignment
        >>> result = align_reads_bwa(
        ...     reference_fasta="ref/genome.fa",
        ...     fastq_files=["reads/sample1.fastq.gz"],
        ...     output_sam="alignments/sample1.sam",
        ...     threads=8,
        ...     sample_id="sample1"
        ... )

        >>> # Paired-end alignment
        >>> result = align_reads_bwa(
        ...     reference_fasta="ref/genome.fa",
        ...     fastq_files=["reads/sample1_R1.fastq.gz", "reads/sample1_R2.fastq.gz"],
        ...     output_sam="alignments/sample1_pe.sam",
        ...     threads=8,
        ...     sample_id="sample1"
        ... )

    Note:
        - BWA-MEM is recommended for reads >70bp and is the default for
          Illumina data. For older datasets with shorter reads, consider
          using BWA-backtrack (bwa aln).
        - The reference genome should be indexed with 'bwa index' before
          running alignment for faster performance.
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
        cmd.extend(
            [
                "-R",
                f"@RG\\tID:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tSM:{sample_id}",
            ]
        )

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


def align_and_sort_bwa(
    reference_fasta: str | Path,
    fastq_files: List[str | Path],
    output_bam: str | Path,
    bwa_threads: int = 8,
    sort_threads: int = 2,
    sort_mem_per_thread: str = "1G",
    sample_id: Optional[str] = None,
) -> Path:
    """Align sequencing reads using BWA-MEM and pipe directly to samtools sort in memory.

    This function provides a dramatic IO optimization by linking the stdout of bwa directly
    to the stdin of samtools sort. By completely bypassing intermediate .sam disk writes,
    this configuration eliminates major thermal bottlenecks on NVMe/SSD drives and provides
    linear CPU-bound scaling across large bioinformatic cohorts.

    Args:
        reference_fasta: Path to reference genome FASTA file.
        fastq_files: List of input FASTQ files.
        output_bam: Output sorted BAM file path.
        bwa_threads: Number of threads for the alignment heavy-lifting.
        sort_threads: Number of threads dedicated strictly to sorting.
        sort_mem_per_thread: Memory boundary per sort thread.
        sample_id: Sample identifier for Read Group tags.

    Returns:
        Path object pointing to the produced, sorted, and indexed BAM.
    """
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    try:
        for tool in ["bwa", "samtools"]:
            subprocess.run([tool], capture_output=True, check=False)
    except FileNotFoundError:
        raise FileNotFoundError("bwa/samtools bioinformatic tools not found on PATH")

    bwa_cmd = ["bwa", "mem", "-t", str(bwa_threads)]
    if sample_id:
        bwa_cmd.extend(["-R", f"@RG\\tID:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tSM:{sample_id}"])
    bwa_cmd.append(str(reference_fasta))
    bwa_cmd.extend([str(f) for f in fastq_files])

    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(sort_threads),
        "-m",
        sort_mem_per_thread,
        "-T",
        f"/tmp/gwas_sort_{output_bam.stem}",
        "-o",
        str(output_bam),
    ]

    logger.info(
        f"Running piped BWA-MEM | samtools sort: output={output_bam.name} (bwa={bwa_threads}t, sort={sort_threads}t)"
    )

    try:
        bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=False)
        sort_proc = subprocess.Popen(
            sort_cmd, stdin=bwa_proc.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=False
        )

        if bwa_proc.stdout:
            bwa_proc.stdout.close()

        _, sort_stderr = sort_proc.communicate()
        bwa_proc.wait()

        if sort_proc.returncode != 0 or bwa_proc.returncode != 0:
            errs = []
            if bwa_proc.returncode != 0:
                errs.append(f"BWA code {bwa_proc.returncode}: (stderr dispatched to devnull)")
            if sort_proc.returncode != 0:
                errs.append(f"Samtools code {sort_proc.returncode}: {sort_stderr.decode('utf-8', errors='replace')}")
            raise subprocess.CalledProcessError(
                sort_proc.returncode or bwa_proc.returncode, "Piped Alignment", stderr="\\n".join(errs).encode()
            )

        logger.info("Indexing generated BAM...")
        index_cmd = ["samtools", "index", str(output_bam)]
        subprocess.run(index_cmd, check=True)

        logger.info(f"Piped alignment + indexing fully completed: {output_bam.name}")
        return output_bam

    except subprocess.CalledProcessError as e:
        err_msg = e.stderr.decode("utf-8", errors="replace") if e.stderr else str(e)
        logger.error(f"Piped alignment execution failed with exit code {e.returncode}:\\n{err_msg}")
        if output_bam.exists():
            output_bam.unlink()
        raise
    except Exception as e:
        logger.error(f"Piped alignment runtime execution failed: {e}")
        if output_bam.exists():
            output_bam.unlink()
        raise


def sort_index_bam(
    input_sam_or_bam: str | Path,
    output_bam: str | Path,
    threads: int = 4,
    sort_mem_per_thread: str = "1G",
) -> Path:
    """Sort and index an alignment file using samtools.

    Args:
        input_sam_or_bam: Input SAM or uncompressed BAM.
        output_bam: Path for the final sorted BAM.
        threads: Number of threads for sorting.
        sort_mem_per_thread: Memory per sort thread (e.g. '1G', '768M'). Defaults
            to '1G', capping total sort RAM at threads × sort_mem_per_thread.
            Samtools default is 768M; with 16 threads that is 12 GB which can
            exhaust RAM on a laptop. Use '512M' for safety on machines with
            16 GB or less.

    Returns:
        Path to the final indexed BAM.
    """
    input_path = Path(input_sam_or_bam)
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. View and Sort (with explicit per-thread memory cap)
        sort_cmd = [
            "samtools",
            "sort",
            "-@",
            str(threads),
            "-m",
            sort_mem_per_thread,
            "-o",
            str(output_bam),
            str(input_path),
        ]
        logger.info(f"Sorting alignments: {output_bam} (-@ {threads} -m {sort_mem_per_thread})")
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

    This utility function verifies that the essential bioinformatics
    command-line tools for read alignment are installed and accessible
    in the system PATH. It is typically used in pre-flight validation
    before running GWAS pipelines.

    Returns:
        bool: True if all required tools (bwa, samtools) are found in
            PATH. False if any tool is missing.

    Example:
        >>> if check_alignment_tools():
        ...     print("All alignment tools available")
        ... else:
        ...     print("Missing required tools - install bwa and samtools")

    Note:
        - bwa: Required for read alignment to reference genome
        - samtools: Required for SAM/BAM file manipulation and indexing
        - Installation: `conda install -c bioconda bwa samtools`
    """
    tools = ["bwa", "samtools"]
    for tool in tools:
        try:
            subprocess.run([tool], capture_output=True, check=False)
        except FileNotFoundError:
            logger.error(f"Required tool not found: {tool}")
            return False
    return True
