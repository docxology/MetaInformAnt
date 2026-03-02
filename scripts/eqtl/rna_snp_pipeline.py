#!/usr/bin/env python3
"""Transcriptome SNP Calling Pipeline.

Re-downloads FASTQs for completed Amalgkit samples, aligns to a reference
genome with HISAT2, calls SNPs with bcftools, and outputs per-biosample VCF
files plus population genetics summaries.

Usage:
    uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
    uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --samples SRR21601882,SRR21601883
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
AMALGKIT_OUTPUT = PROJECT_ROOT / "output" / "amalgkit"
EQTL_OUTPUT = PROJECT_ROOT / "output" / "eqtl"
DOWNLOAD_SCRIPT = PROJECT_ROOT / "scripts" / "rna" / "download_ena.py"

# NCBI genome accession map per species
GENOME_ACCESSIONS = {
    "amellifera": "GCF_003254395.2",
    "acromyrmex_echinatior": "GCF_000204515.1",
}

# Minimum tools required
REQUIRED_TOOLS = ["hisat2", "hisat2-build", "samtools", "bcftools", "wget"]


def check_tools() -> bool:
    """Verify all required tools are on PATH."""
    missing = []
    for tool in REQUIRED_TOOLS:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        logger.error(f"Missing required tools: {', '.join(missing)}")
        logger.error("Install with: sudo apt install hisat2 samtools bcftools")
        return False
    return True


def find_completed_samples(species: str, max_samples: int | None = None) -> list[str]:
    """Find sample SRR IDs that have completed Amalgkit quantification."""
    quant_dirs = []
    # Search in multiple possible locations
    for pattern in [
        AMALGKIT_OUTPUT / species / "quant" / "quant",
        AMALGKIT_OUTPUT / species / "work" / "quant",
    ]:
        if pattern.exists():
            for d in sorted(pattern.iterdir()):
                if d.is_dir() and (d / "abundance.tsv").exists():
                    quant_dirs.append(d.name)

    # Deduplicate
    samples = sorted(set(quant_dirs))
    if max_samples and len(samples) > max_samples:
        samples = samples[:max_samples]
    logger.info(f"Found {len(samples)} completed samples for {species}")
    return samples


def find_reference_genome(species: str) -> Path | None:
    """Find the reference genome FASTA for a species."""
    genome_dir = AMALGKIT_OUTPUT / species / "genome"
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
    """Decompress a .gz file if the uncompressed version doesn't exist."""
    if not gz_path.name.endswith(".gz"):
        return gz_path
    decompressed = gz_path.with_suffix("")  # Remove .gz
    if decompressed.exists() and decompressed.stat().st_size > 0:
        return decompressed
    logger.info(f"Decompressing {gz_path.name}...")
    subprocess.run(["gunzip", "-k", str(gz_path)], check=True)
    return decompressed


def build_hisat2_index(genome_fasta: Path, output_dir: Path) -> Path:
    """Build HISAT2 index from reference genome."""
    index_prefix = output_dir / "genome"
    marker = output_dir / ".index_built"

    if marker.exists():
        logger.info("HISAT2 index already built, skipping")
        return index_prefix

    output_dir.mkdir(parents=True, exist_ok=True)
    fasta = decompress_if_needed(genome_fasta)

    logger.info(f"Building HISAT2 index from {fasta.name}...")
    cmd = ["hisat2-build", "-p", "4", str(fasta), str(index_prefix)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"hisat2-build failed:\n{result.stderr}")
        raise RuntimeError("hisat2-build failed")

    marker.touch()
    logger.info("HISAT2 index built successfully")
    return index_prefix


def download_fastq(srr_id: str, output_dir: Path, species: str = "") -> list[Path]:
    """Download FASTQ files for a sample from ENA, or reuse local amalgkit copies."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if FASTQs already exist in our output
    existing = sorted(output_dir.glob(f"{srr_id}*.fastq.gz"))
    if existing:
        logger.info(f"FASTQs already present for {srr_id}: {len(existing)} files")
        return existing

    # Check if FASTQs are available locally from amalgkit
    if species:
        amalgkit_fastq_dir = AMALGKIT_OUTPUT / species / "fastq"
        local_fqs = sorted(amalgkit_fastq_dir.glob(f"{srr_id}*.fastq.gz"))
        if local_fqs:
            logger.info(f"Reusing {len(local_fqs)} local amalgkit FASTQ(s) for {srr_id}")
            linked = []
            for fq in local_fqs:
                dest = output_dir / fq.name
                if not dest.exists():
                    os.symlink(fq.resolve(), dest)
                linked.append(dest)
            return linked

    logger.info(f"Downloading FASTQs for {srr_id}...")
    cmd = ["python3", str(DOWNLOAD_SCRIPT), srr_id, str(output_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

    if result.returncode != 0:
        logger.warning(f"Download failed for {srr_id}: {result.stderr[:200]}")
        return []

    fastqs = sorted(output_dir.glob(f"{srr_id}*.fastq.gz"))
    logger.info(f"Downloaded {len(fastqs)} FASTQ files for {srr_id}")
    return fastqs


def align_reads(
    fastqs: list[Path],
    index_prefix: Path,
    output_bam: Path,
    threads: int = 4,
) -> bool:
    """Align RNA-seq reads with HISAT2 and sort with samtools."""
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
    hisat2_proc = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam)]
    sort_proc = subprocess.Popen(sort_cmd, stdin=hisat2_proc.stdout, stderr=subprocess.PIPE)

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
    """Call variants using bcftools mpileup + call."""
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
        "bcftools", "mpileup",
        "-f", str(fasta),
        "-Q", "20",        # min base quality
        "-q", "20",        # min mapping quality
        "--max-depth", "10000",
        str(bam_path),
    ]
    call_cmd = [
        "bcftools", "call",
        "-mv",              # multiallelic caller, output variants only
        "--ploidy", "2",
        "-Oz",              # compressed VCF output
        "-o", str(output_vcf),
    ]

    mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stderr=subprocess.PIPE)

    mpileup_proc.stdout.close()
    call_stderr = call_proc.communicate()[1].decode()
    mpileup_stderr = mpileup_proc.stderr.read().decode()
    mpileup_proc.wait()

    if call_proc.returncode != 0:
        logger.error(f"bcftools call failed:\n{call_stderr[:500]}")
        return False

    # Index VCF
    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)

    # Count variants
    count_result = subprocess.run(
        ["bcftools", "stats", str(output_vcf)],
        capture_output=True, text=True,
    )
    for line in count_result.stdout.split("\n"):
        if line.startswith("SN") and "number of records" in line:
            logger.info(f"  {line.strip()}")

    return True


def filter_variants(input_vcf: Path, output_vcf: Path) -> bool:
    """Apply quality filters to VCF."""
    if output_vcf.exists() and output_vcf.stat().st_size > 0:
        return True

    cmd = [
        "bcftools", "filter",
        "-s", "LowQual",
        "-e", "QUAL<30 || DP<10",
        "-Oz", "-o", str(output_vcf),
        str(input_vcf),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Filter failed: {result.stderr[:300]}")
        return False

    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
    return True


def compute_sample_stats(vcf_path: Path, output_json: Path) -> dict:
    """Compute per-sample variant statistics."""
    result = subprocess.run(
        ["bcftools", "stats", str(vcf_path)],
        capture_output=True, text=True,
    )

    stats = {
        "vcf_file": str(vcf_path),
        "n_records": 0, "n_snps": 0, "n_indels": 0,
        "n_pass": 0, "ts_tv_ratio": 0.0,
    }

    for line in result.stdout.split("\n"):
        if not line.startswith("SN"):
            continue
        if "number of records:" in line:
            stats["n_records"] = int(line.split("\t")[-1])
        elif "number of SNPs:" in line:
            stats["n_snps"] = int(line.split("\t")[-1])
        elif "number of indels:" in line:
            stats["n_indels"] = int(line.split("\t")[-1])

    # Ti/Tv from bcftools stats
    for line in result.stdout.split("\n"):
        if line.startswith("TSTV"):
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    stats["ts_tv_ratio"] = float(parts[4])
                except (ValueError, IndexError):
                    pass

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(stats, f, indent=2)

    return stats


def merge_vcfs(vcf_files: list[Path], output_vcf: Path) -> bool:
    """Merge per-sample VCFs into a population VCF."""
    if len(vcf_files) < 2:
        logger.warning("Need at least 2 VCFs to merge")
        if vcf_files:
            shutil.copy2(vcf_files[0], output_vcf)
            subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
        return bool(vcf_files)

    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bcftools", "merge",
        "-Oz", "-o", str(output_vcf),
    ] + [str(v) for v in vcf_files]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Merge failed: {result.stderr[:300]}")
        return False

    subprocess.run(["bcftools", "index", str(output_vcf)], check=True)
    logger.info(f"Merged {len(vcf_files)} VCFs → {output_vcf}")
    return True


def compute_allele_frequencies(merged_vcf: Path, output_tsv: Path) -> None:
    """Extract per-site allele frequencies from merged VCF."""
    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
        str(merged_vcf),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as f:
        f.write("chrom\tpos\tref\talt\taf\n")
        f.write(result.stdout)

    n_lines = result.stdout.count("\n")
    logger.info(f"Wrote allele frequencies for {n_lines} variants → {output_tsv}")


def compute_popgen_summary(merged_vcf: Path, sample_stats: list[dict], output_json: Path) -> dict:
    """Compute population-level genetics summary."""
    result = subprocess.run(
        ["bcftools", "stats", str(merged_vcf)],
        capture_output=True, text=True,
    )

    summary = {
        "n_samples": len(sample_stats),
        "total_variants": 0,
        "total_snps": 0,
        "total_indels": 0,
        "ts_tv_ratio": 0.0,
        "per_sample_stats": sample_stats,
    }

    for line in result.stdout.split("\n"):
        if not line.startswith("SN"):
            continue
        if "number of records:" in line:
            summary["total_variants"] = int(line.split("\t")[-1])
        elif "number of SNPs:" in line:
            summary["total_snps"] = int(line.split("\t")[-1])
        elif "number of indels:" in line:
            summary["total_indels"] = int(line.split("\t")[-1])

    for line in result.stdout.split("\n"):
        if line.startswith("TSTV"):
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    summary["ts_tv_ratio"] = float(parts[4])
                except (ValueError, IndexError):
                    pass

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def run_pipeline(
    species: str,
    sample_ids: list[str] | None = None,
    n_samples: int = 3,
    threads: int = 4,
    cleanup_fastq: bool = True,
) -> dict:
    """Run the full transcriptome SNP calling pipeline.

    Args:
        species: Species name (must match amalgkit output dir).
        sample_ids: Optional explicit list of SRR IDs.
        n_samples: Number of samples to process (if sample_ids not given).
        threads: Threads for alignment and indexing.
        cleanup_fastq: Delete FASTQs after alignment to save disk space.

    Returns:
        Pipeline run summary dict.
    """
    start_time = time.time()
    logger.info("=" * 70)
    logger.info(f"Transcriptome SNP Calling Pipeline — {species}")
    logger.info("=" * 70)

    # Check tools
    if not check_tools():
        sys.exit(1)

    # Find or use provided samples
    if sample_ids:
        samples = sample_ids
    else:
        samples = find_completed_samples(species, max_samples=n_samples)

    if not samples:
        logger.error(f"No completed samples found for {species}")
        sys.exit(1)

    logger.info(f"Processing {len(samples)} samples: {samples}")

    # Find reference genome
    ref_genome = find_reference_genome(species)
    if not ref_genome:
        logger.error(f"No reference genome found for {species}")
        sys.exit(1)
    logger.info(f"Reference genome: {ref_genome}")

    # Output directories
    species_out = EQTL_OUTPUT / species
    index_dir = species_out / "index"
    pop_dir = species_out / "population"
    log_dir = species_out / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Set up file logging
    fh = logging.FileHandler(log_dir / "pipeline.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    # Step 1: Build HISAT2 index
    logger.info("\n[Step 1/6] Building HISAT2 index...")
    index_prefix = build_hisat2_index(ref_genome, index_dir)

    # Process each sample
    filtered_vcfs = []
    all_sample_stats = []

    for i, srr_id in enumerate(samples, 1):
        logger.info(f"\n{'='*50}")
        logger.info(f"Sample {i}/{len(samples)}: {srr_id}")
        logger.info(f"{'='*50}")

        sample_dir = species_out / "samples" / srr_id
        fastq_dir = sample_dir / "fastq"
        bam_path = sample_dir / "aligned.bam"
        raw_vcf = sample_dir / "variants_raw.vcf.gz"
        filt_vcf = sample_dir / "variants.vcf.gz"
        stats_json = sample_dir / "variant_stats.json"

        # Skip if already fully processed
        if filt_vcf.exists() and filt_vcf.stat().st_size > 0 and stats_json.exists():
            logger.info(f"Sample {srr_id} already processed, loading cached stats")
            with open(stats_json) as f:
                stats = json.load(f)
            stats["sample_id"] = srr_id
            all_sample_stats.append(stats)
            filtered_vcfs.append(filt_vcf)
            logger.info(f"  {srr_id}: {stats.get('n_snps', '?')} SNPs (cached)")
            continue

        # Step 2: Download FASTQs (checks local amalgkit copies first)
        logger.info(f"\n[Step 2/6] Downloading FASTQs for {srr_id}...")
        fastqs = download_fastq(srr_id, fastq_dir, species=species)
        if not fastqs:
            logger.warning(f"Skipping {srr_id}: no FASTQs downloaded")
            continue

        # Step 3: Align reads
        logger.info(f"\n[Step 3/6] Aligning {srr_id}...")
        ref_fasta = decompress_if_needed(ref_genome)
        ok = align_reads(fastqs, index_prefix, bam_path, threads=threads)
        if not ok:
            logger.warning(f"Skipping {srr_id}: alignment failed")
            continue

        # Step 4: Call variants
        logger.info(f"\n[Step 4/6] Calling variants for {srr_id}...")
        ok = call_variants(bam_path, ref_genome, raw_vcf)
        if not ok:
            logger.warning(f"Skipping {srr_id}: variant calling failed")
            continue

        # Step 5: Filter variants
        logger.info(f"\n[Step 5/6] Filtering variants for {srr_id}...")
        ok = filter_variants(raw_vcf, filt_vcf)
        if not ok:
            logger.warning(f"Skipping {srr_id}: filtering failed")
            continue

        # Compute per-sample stats
        stats = compute_sample_stats(filt_vcf, stats_json)
        stats["sample_id"] = srr_id
        all_sample_stats.append(stats)
        filtered_vcfs.append(filt_vcf)

        logger.info(f"  {srr_id}: {stats['n_snps']} SNPs, {stats['n_indels']} indels")

        # Cleanup FASTQs to save disk
        if cleanup_fastq and fastq_dir.exists():
            shutil.rmtree(fastq_dir)
            logger.info(f"  Cleaned up FASTQs for {srr_id}")

    # Step 6: Merge and population analysis
    logger.info(f"\n[Step 6/6] Merging {len(filtered_vcfs)} sample VCFs...")
    if filtered_vcfs:
        merged_vcf = pop_dir / "merged.vcf.gz"
        merge_vcfs(filtered_vcfs, merged_vcf)

        # Allele frequencies
        compute_allele_frequencies(merged_vcf, pop_dir / "allele_freqs.tsv")

        # Population summary
        pop_summary = compute_popgen_summary(
            merged_vcf, all_sample_stats, pop_dir / "popgen_summary.json"
        )
    else:
        pop_summary = {"n_samples": 0, "error": "No samples completed successfully"}

    elapsed = time.time() - start_time

    run_summary = {
        "species": species,
        "n_samples_requested": len(samples),
        "n_samples_completed": len(filtered_vcfs),
        "samples": [s["sample_id"] for s in all_sample_stats],
        "total_snps": pop_summary.get("total_snps", 0),
        "total_indels": pop_summary.get("total_indels", 0),
        "ts_tv_ratio": pop_summary.get("ts_tv_ratio", 0.0),
        "elapsed_seconds": round(elapsed, 1),
        "output_dir": str(species_out),
    }

    with open(species_out / "run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)

    logger.info("\n" + "=" * 70)
    logger.info("Pipeline Complete!")
    logger.info(f"  Species: {species}")
    logger.info(f"  Samples: {run_summary['n_samples_completed']}/{run_summary['n_samples_requested']}")
    logger.info(f"  SNPs: {run_summary['total_snps']}")
    logger.info(f"  Ti/Tv: {run_summary['ts_tv_ratio']}")
    logger.info(f"  Time: {elapsed:.0f}s")
    logger.info(f"  Output: {species_out}")
    logger.info("=" * 70)

    return run_summary


def load_config(config_path: str) -> dict:
    """Load pipeline configuration from YAML file."""
    import yaml

    with open(config_path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Transcriptome SNP Calling Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --config config/eqtl/eqtl_amellifera.yaml
  %(prog)s --species amellifera --n-samples 3
  %(prog)s --species amellifera --samples SRR21601882,SRR21601883
        """,
    )
    parser.add_argument("--config", help="Path to YAML config file")
    parser.add_argument("--species", help="Species name (e.g. amellifera)")
    parser.add_argument("--n-samples", type=int, default=3, help="Number of samples to process")
    parser.add_argument("--samples", help="Comma-separated list of SRR IDs to process")
    parser.add_argument("--threads", type=int, default=4, help="Threads for alignment")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep FASTQ files after alignment")
    args = parser.parse_args()

    # Load from config file if provided
    if args.config:
        cfg = load_config(args.config)
        species = cfg.get("species", args.species)
        samples_cfg = cfg.get("samples", {})
        if samples_cfg.get("mode") == "explicit" and "explicit_ids" in samples_cfg:
            sample_ids = samples_cfg["explicit_ids"]
        else:
            sample_ids = None
        n_samples = samples_cfg.get("max_samples", args.n_samples)
        threads = cfg.get("alignment", {}).get("threads", args.threads)
        cleanup = cfg.get("output", {}).get("cleanup_fastq", not args.no_cleanup)
    else:
        if not args.species:
            parser.error("--species or --config is required")
        species = args.species
        sample_ids = args.samples.split(",") if args.samples else None
        n_samples = args.n_samples
        threads = args.threads
        cleanup = not args.no_cleanup

    run_pipeline(
        species=species,
        sample_ids=sample_ids,
        n_samples=n_samples,
        threads=threads,
        cleanup_fastq=cleanup,
    )


if __name__ == "__main__":
    main()

