"""External-tool wrappers for transcriptome SNP calling (HISAT2 + bcftools).

Extracted from ``scripts/eqtl/rna_snp_pipeline.py`` so the pipeline logic is
reusable and testable. Each function wraps a real external tool and checks
for pre-existing outputs (idempotent resume). Tools required on PATH:
``hisat2``, ``hisat2-build``, ``samtools``, ``bcftools``, ``curl``.

Example:
    >>> from metainformant.eqtl.variant_calling import check_tools
    >>> isinstance(check_tools(), bool)
    True
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.retrieval.ena_downloader import ENADownloader

logger = logging.getLogger(__name__)

# NCBI genome accession map per species
GENOME_ACCESSIONS = {
    "amellifera": "GCF_003254395.2",
    "acromyrmex_echinatior": "GCF_000204515.1",
}

# Minimum tools required
REQUIRED_TOOLS = ["hisat2", "hisat2-build", "samtools", "bcftools", "curl"]


def check_tools() -> bool:
    """Verify all required external tools are on PATH.

    Returns:
        True when every tool in REQUIRED_TOOLS is resolvable, else False.
    """
    missing = []
    for tool in REQUIRED_TOOLS:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        logger.error(f"Missing required tools: {', '.join(missing)}")
        logger.error("Install with: sudo apt install hisat2 samtools bcftools")
        return False
    return True


def find_completed_samples(
    amalgkit_output: Path,
    species: str,
    max_samples: int | None = None,
) -> list[str]:
    """Find sample SRR IDs that have completed Amalgkit quantification.

    Args:
        amalgkit_output: Amalgkit data root (AMALGKIT_DATA_ROOT).
        species: Species name (must match amalgkit output dir).
        max_samples: Optional cap on the number of samples returned.

    Returns:
        Sorted, deduplicated list of sample directory names.
    """
    quant_dirs = []
    # Search in multiple possible locations
    for pattern in [
        amalgkit_output / species / "quant" / "quant",
        amalgkit_output / species / "work" / "quant",
    ]:
        if pattern.exists():
            for d in sorted(pattern.iterdir()):
                if d.is_dir() and find_quantification_file(d, d.name) is not None:
                    quant_dirs.append(d.name)

    # Deduplicate
    samples = sorted(set(quant_dirs))
    if max_samples and len(samples) > max_samples:
        samples = samples[:max_samples]
    logger.info(f"Found {len(samples)} completed samples for {species}")
    return samples


def find_reference_genome(amalgkit_output: Path, species: str) -> Path | None:
    """Find the reference genome FASTA for a species.

    Args:
        amalgkit_output: Amalgkit data root (AMALGKIT_DATA_ROOT).
        species: Species name (must match amalgkit output dir).

    Returns:
        Path to the genomic FASTA (plain or .fna.gz), or None when absent.
    """
    genome_dir = amalgkit_output / species / "genome"
    if not genome_dir.exists():
        return None
    # Look for the genomic FASTA (not cds or rna)
    for f in genome_dir.iterdir():
        if "genomic.fna" in f.name and "cds" not in f.name and "rna" not in f.name:
            return f
    # Fallback to any fna.gz
    for f in genome_dir.iterdir():
        if f.name.endswith(".fna.gz"):
            return f
    return None


def decompress_if_needed(gz_path: Path) -> Path:
    """Decompress a .gz file if the uncompressed version doesn't exist.

    Args:
        gz_path: Path to a possibly-gzipped file.

    Returns:
        Path to the usable (decompressed or original) file.

    Raises:
        RuntimeError: If gunzip fails.
    """
    if not gz_path.name.endswith(".gz"):
        return gz_path
    decompressed = gz_path.with_suffix("")  # Remove .gz
    if decompressed.exists() and decompressed.stat().st_size > 0:
        return decompressed
    logger.info(f"Decompressing {gz_path.name}...")
    result = subprocess.run(
        ["gunzip", "-k", str(gz_path)], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"gunzip failed for {gz_path}: {result.stderr[:300]}")
    return decompressed


def build_hisat2_index(genome_fasta: Path, output_dir: Path, threads: int = 4) -> Path:
    """Build HISAT2 index from reference genome (idempotent).

    Args:
        genome_fasta: Reference genome FASTA (possibly gzipped).
        output_dir: Directory for the index (marker file ``.index_built``).
        threads: Threads for hisat2-build.

    Returns:
        Index prefix path (output_dir / "genome").

    Raises:
        RuntimeError: If hisat2-build fails.
    """
    index_prefix = output_dir / "genome"
    marker = output_dir / ".index_built"

    if marker.exists():
        logger.info("HISAT2 index already built, skipping")
        return index_prefix

    output_dir.mkdir(parents=True, exist_ok=True)
    fasta = decompress_if_needed(genome_fasta)

    logger.info(f"Building HISAT2 index from {fasta.name}...")
    cmd = ["hisat2-build", "-p", str(threads), str(fasta), str(index_prefix)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"hisat2-build failed:\n{result.stderr}")
        raise RuntimeError("hisat2-build failed")

    marker.touch()
    logger.info("HISAT2 index built successfully")
    return index_prefix


def download_fastq(
    srr_id: str,
    output_dir: Path,
    species: str = "",
    amalgkit_output: Path | None = None,
) -> list[Path]:
    """Download FASTQ files for a sample from ENA, or reuse local copies.

    Order of preference: (1) FASTQs already in output_dir, (2) local
    amalgkit FASTQs (symlinked), (3) fresh ENA download.

    Args:
        srr_id: Sample run accession.
        output_dir: Destination directory.
        species: Optional species name to enable local reuse.
        amalgkit_output: Amalgkit data root for local FASTQ reuse.

    Returns:
        List of FASTQ paths, or an empty list when download fails.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if FASTQs already exist in our output
    existing = sorted(output_dir.glob(f"{srr_id}*.fastq.gz"))
    if existing:
        logger.info(f"FASTQs already present for {srr_id}: {len(existing)} files")
        return existing

    # Check if FASTQs are available locally from amalgkit
    if species and amalgkit_output is not None:
        amalgkit_fastq_dir = amalgkit_output / species / "fastq"
        local_fqs = sorted(amalgkit_fastq_dir.glob(f"{srr_id}*.fastq.gz"))
        if local_fqs:
            logger.info(
                f"Reusing {len(local_fqs)} local amalgkit FASTQ(s) for {srr_id}"
            )
            linked = []
            for fq in local_fqs:
                dest = output_dir / fq.name
                if not dest.exists():
                    os.symlink(fq.resolve(), dest)
                linked.append(dest)
            return linked

    logger.info(f"Downloading FASTQs for {srr_id} via ENA...")
    downloader = ENADownloader(timeout=1800, retries=3)
    success, message, fastqs = downloader.download_run(srr_id, output_dir)
    if not success:
        logger.warning(f"Download failed for {srr_id}: {message[:200]}")
        return []
    logger.info(f"{message} for {srr_id}")
    return fastqs


def align_reads(
    fastqs: list[Path],
    index_prefix: Path,
    output_bam: Path,
    threads: int = 4,
) -> bool:
    """Align RNA-seq reads with HISAT2 and sort/index with samtools.

    Args:
        fastqs: One (single-end) or two (paired-end) FASTQ paths.
        index_prefix: HISAT2 index prefix.
        output_bam: Destination BAM path (sorted + indexed).
        threads: Threads for alignment and sorting.

    Returns:
        True on success (or when a non-empty BAM already exists).
    """
    if output_bam.exists() and output_bam.stat().st_size > 0:
        logger.info(f"BAM already exists: {output_bam.name}")
        return True

    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Build hisat2 command
    hisat2_cmd = ["hisat2", "--dta", "-p", str(threads), "-x", str(index_prefix)]

    if len(fastqs) == 2:
        hisat2_cmd += ["-1", str(fastqs[0]), "-2", str(fastqs[1])]
    else:
        hisat2_cmd += ["-U", str(fastqs[0])]

    # Pipe through samtools sort
    logger.info(f"Aligning {len(fastqs)} FASTQ file(s)...")
    hisat2_proc = subprocess.Popen(
        hisat2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam)]
    sort_proc = subprocess.Popen(
        sort_cmd, stdin=hisat2_proc.stdout, stderr=subprocess.PIPE
    )

    hisat2_proc.stdout.close()
    sort_stderr = sort_proc.communicate()[1].decode()
    hisat2_stderr = hisat2_proc.stderr.read().decode()
    hisat2_proc.wait()

    if hisat2_proc.returncode != 0:
        logger.error(f"HISAT2 failed:\n{hisat2_stderr[:500]}")
        return False

    if sort_proc.returncode != 0:
        logger.error(f"samtools sort failed:\n{sort_stderr[:500]}")
        return False

    # Index the BAM
    subprocess.run(["samtools", "index", str(output_bam)], check=True)

    # Log alignment stats from HISAT2 stderr
    for line in hisat2_stderr.strip().split("\n"):
        if "overall alignment rate" in line or "aligned" in line.lower():
            logger.info(f"  {line.strip()}")

    logger.info(f"Alignment complete: {output_bam}")
    return True


def call_variants(bam_path: Path, ref_fasta: Path, output_vcf: Path) -> bool:
    """Call variants using bcftools mpileup + call (idempotent).

    Args:
        bam_path: Input BAM (must be indexed).
        ref_fasta: Reference genome FASTA (possibly gzipped).
        output_vcf: Destination compressed VCF.

    Returns:
        True on success (or when a non-empty VCF already exists).
    """
    if output_vcf.exists() and output_vcf.stat().st_size > 0:
        logger.info(f"VCF already exists: {output_vcf.name}")
        return True

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    fasta = decompress_if_needed(ref_fasta)

    # Index reference if needed
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        logger.info("Indexing reference FASTA...")
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)

    logger.info(f"Calling variants from {bam_path.name}...")

    # bcftools mpileup | bcftools call
    mpileup_cmd = [
        "bcftools",
        "mpileup",
        "-f",
        str(fasta),
        "-Q",
        "20",  # min base quality
        "-q",
        "20",  # min mapping quality
        "--max-depth",
        "10000",
        str(bam_path),
    ]
    call_cmd = [
        "bcftools",
        "call",
        "-mv",  # multiallelic caller, output variants only
        "--ploidy",
        "2",
        "-Oz",  # compressed VCF output
        "-o",
        str(output_vcf),
    ]

    mpileup_proc = subprocess.Popen(
        mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    call_proc = subprocess.Popen(
        call_cmd, stdin=mpileup_proc.stdout, stderr=subprocess.PIPE
    )

    mpileup_proc.stdout.close()
    call_stderr = call_proc.communicate()[1].decode()
    mpileup_proc.stderr.read().decode()
    mpileup_proc.wait()

    if call_proc.returncode != 0:
        logger.error(f"bcftools call failed:\n{call_stderr[:500]}")
        return False

    # Index VCF
    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)

    # Count variants
    count_result = subprocess.run(
        ["bcftools", "stats", str(output_vcf)],
        capture_output=True,
        text=True,
    )
    for line in count_result.stdout.split("\n"):
        if line.startswith("SN") and "number of records" in line:
            logger.info(f"  {line.strip()}")

    return True


def filter_variants(input_vcf: Path, output_vcf: Path) -> bool:
    """Apply quality filters to VCF (QUAL<30 || DP<10 -> LowQual; idempotent).

    Args:
        input_vcf: Input compressed VCF.
        output_vcf: Destination compressed, filtered VCF.

    Returns:
        True on success (or when a non-empty output already exists).
    """
    if output_vcf.exists() and output_vcf.stat().st_size > 0:
        return True

    cmd = [
        "bcftools",
        "filter",
        "-s",
        "LowQual",
        "-e",
        "QUAL<30 || DP<10",
        "-Oz",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Filter failed: {result.stderr[:300]}")
        return False

    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
    return True


def merge_vcfs(vcf_files: list[Path], output_vcf: Path) -> bool:
    """Merge per-sample VCFs into a population VCF.

    Args:
        vcf_files: Per-sample compressed VCFs (bcftools-indexed).
        output_vcf: Destination merged compressed VCF.

    Returns:
        True on success; a single input VCF is copied when only one exists.
    """
    if len(vcf_files) < 2:
        logger.warning("Need at least 2 VCFs to merge")
        if vcf_files:
            shutil.copy2(vcf_files[0], output_vcf)
            subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
        return bool(vcf_files)

    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bcftools",
        "merge",
        "-Oz",
        "-o",
        str(output_vcf),
    ] + [str(v) for v in vcf_files]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Merge failed: {result.stderr[:300]}")
        return False

    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
    logger.info(f"Merged {len(vcf_files)} VCFs -> {output_vcf}")
    return True
