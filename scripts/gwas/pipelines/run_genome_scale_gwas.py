#!/usr/bin/env python3
"""Run genome-scale GWAS with comprehensive visualization.

This script orchestrates the complete GWAS workflow from SRA download
through variant calling, association testing, and comprehensive visualization.

Usage:
    python scripts/run_genome_scale_gwas.py --config config/gwas/gwas_amellifera.yaml
    python scripts/run_genome_scale_gwas.py --download-only  # Just download data
    python scripts/run_genome_scale_gwas.py --skip-download  # Use existing data
"""

import argparse
import logging
import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core.io import dump_json
from metainformant.gwas import (
    check_sra_tools_available,
    download_sra_run,
    execute_gwas_workflow,
    generate_all_plots,
    load_gwas_config,
)

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


def check_dependencies():
    """Check if required bioinformatics tools are available."""
    import shutil

    tools = {
        "fasterq-dump": check_sra_tools_available(),
        "bwa": shutil.which("bwa") is not None,
        "samtools": shutil.which("samtools") is not None,
        "bcftools": shutil.which("bcftools") is not None,
    }

    logger.info("Checking dependencies:")
    for tool, available in tools.items():
        status = "✓" if available else "✗"
        logger.info(f"  {status} {tool}")

    return tools


def download_sra_data(accessions: list[str], output_dir: Path, threads: int = 8):
    """Download SRA data for specified accessions."""
    logger.info(f"Downloading {len(accessions)} SRA samples to {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for acc in accessions:
        logger.info(f"Downloading {acc}...")
        result = download_sra_run(
            sra_accession=acc,
            dest_dir=str(output_dir),
            threads=threads,
        )
        results.append(result)

        if result.get("status") == "success":
            logger.info(f"  ✓ {acc} downloaded successfully")
        else:
            logger.error(f"  ✗ {acc} failed: {result.get('error', 'Unknown error')}")

    successful = sum(1 for r in results if r.get("status") == "success")
    logger.info(f"Downloaded {successful}/{len(accessions)} samples successfully")

    return results


def align_reads(fastq_dir: Path, reference: Path, output_dir: Path, threads: int = 8):
    """Align FASTQ reads to reference genome using BWA."""
    import subprocess

    logger.info(f"Aligning reads from {fastq_dir} to {reference}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find FASTQ pairs
    fastq_files = sorted(fastq_dir.glob("*_1.fastq"))

    if not fastq_files:
        logger.error("No FASTQ files found")
        return []

    bam_files = []
    for fq1 in fastq_files:
        sample_name = fq1.stem.replace("_1", "")
        fq2 = fq1.parent / f"{sample_name}_2.fastq"

        if not fq2.exists():
            logger.warning(f"Skipping {sample_name}: paired file not found")
            continue

        output_bam = output_dir / f"{sample_name}.sorted.bam"

        if output_bam.exists():
            logger.info(f"  ✓ {sample_name}.bam already exists, skipping")
            bam_files.append(output_bam)
            continue

        logger.info(f"  Aligning {sample_name}...")

        # Run BWA MEM and SAMtools sort
        cmd = f"bwa mem -t {threads} {reference} {fq1} {fq2} | samtools sort -@ {threads} -o {output_bam}"

        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            subprocess.run(f"samtools index {output_bam}", shell=True, check=True)
            logger.info(f"  ✓ {sample_name} aligned successfully")
            bam_files.append(output_bam)
        except subprocess.CalledProcessError as e:
            logger.error(f"  ✗ {sample_name} failed: {e}")

    logger.info(f"Aligned {len(bam_files)} samples")
    return bam_files


def call_variants(bam_files: list[Path], reference: Path, output_vcf: Path, threads: int = 8):
    """Call variants from aligned BAM files using bcftools."""
    import subprocess

    logger.info(f"Calling variants from {len(bam_files)} BAM files")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    if output_vcf.exists():
        logger.info(f"  ✓ VCF already exists: {output_vcf}")
        return {"status": "success", "vcf_file": str(output_vcf)}

    # Join BAM files
    bam_list = " ".join(str(b) for b in bam_files)

    # bcftools mpileup + call
    cmd = f"bcftools mpileup -f {reference} -Ou {bam_list} | bcftools call -mv -Oz -o {output_vcf}"

    try:
        logger.info("  Running bcftools mpileup and call...")
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f"bcftools index {output_vcf}", shell=True, check=True)

        # Get variant count
        result = subprocess.run(f"bcftools view -H {output_vcf} | wc -l", shell=True, capture_output=True, text=True)
        n_variants = int(result.stdout.strip())

        logger.info(f"  ✓ Called {n_variants:,} variants")

        return {
            "status": "success",
            "vcf_file": str(output_vcf),
            "num_variants": n_variants,
        }
    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Variant calling failed: {e}")
        return {"status": "failed", "error": str(e)}


def main():
    parser = argparse.ArgumentParser(description="Run genome-scale GWAS workflow")
    parser.add_argument(
        "--config", type=Path, default=Path("config/gwas/gwas_amellifera.yaml"), help="GWAS configuration file"
    )
    parser.add_argument("--download-only", action="store_true", help="Only download SRA data, don't run analysis")
    parser.add_argument("--skip-download", action="store_true", help="Skip SRA download (use existing data)")
    parser.add_argument("--skip-align", action="store_true", help="Skip alignment (use existing BAM files)")
    parser.add_argument("--skip-calling", action="store_true", help="Skip variant calling (use existing VCF)")
    parser.add_argument(
        "--comprehensive-plots", action="store_true", default=True, help="Generate comprehensive visualization suite"
    )
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument(
        "--sra-accessions",
        nargs="+",
        default=["SRR2096937", "SRR2096938", "SRR2096939"],
        help="SRA accessions to download",
    )
    args = parser.parse_args()

    # Setup logging
    logging.getLogger().setLevel(logging.INFO)

    logger.info("═" * 80)
    logger.info("GENOME-SCALE GWAS WORKFLOW")
    logger.info("═" * 80)

    # Check dependencies
    deps = check_dependencies()

    if not all(deps.values()):
        logger.warning("Some dependencies missing. See docs/gwas/INSTALL_TOOLS.md for installation")
        if not args.skip_download and not deps["fasterq-dump"]:
            logger.error("fasterq-dump required for SRA download")
            return 1
        if not args.skip_align and not all([deps["bwa"], deps["samtools"]]):
            logger.error("BWA and SAMtools required for alignment")
            return 1
        if not args.skip_calling and not deps["bcftools"]:
            logger.error("bcftools required for variant calling")
            return 1

    # Load config
    logger.info(f"\nLoading configuration: {args.config}")
    config = load_gwas_config(args.config)

    # Define directories
    sra_dir = Path("data/raw/sra")
    aligned_dir = Path("data/aligned")
    variants_dir = Path("data/variants/amellifera/real")

    # Step 1: Download SRA data
    if not args.skip_download and not args.skip_align:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 1: Download SRA Data")
        logger.info("─" * 80)
        download_results = download_sra_data(args.sra_accessions, sra_dir, args.threads)

        if args.download_only:
            logger.info("\n✓ Download complete (download-only mode)")
            return 0
    else:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 1: Download SRA Data [SKIPPED]")
        logger.info("─" * 80)

    # Step 2: Align reads
    if not args.skip_align and not args.skip_calling:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 2: Align Reads with BWA")
        logger.info("─" * 80)

        reference = config.work_dir.parent / "genome" / "genomic.fna"
        if not reference.exists():
            logger.error(f"Reference genome not found: {reference}")
            logger.info("Run: python -m metainformant.gwas download_reference first")
            return 1

        bam_files = align_reads(sra_dir, reference, aligned_dir, args.threads)

        if not bam_files:
            logger.error("No BAM files created")
            return 1
    else:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 2: Align Reads [SKIPPED]")
        logger.info("─" * 80)
        bam_files = list(aligned_dir.glob("*.sorted.bam"))

    # Step 3: Call variants
    if not args.skip_calling:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 3: Call Variants with bcftools")
        logger.info("─" * 80)

        output_vcf = variants_dir / "genome_scale_cohort.vcf.gz"
        reference = config.work_dir.parent / "genome" / "genomic.fna"

        call_result = call_variants(bam_files, reference, output_vcf, args.threads)

        if call_result.get("status") != "success":
            logger.error("Variant calling failed")
            return 1

        # Update config to use this VCF
        config.variants = {"vcf_files": [str(output_vcf)]}
    else:
        logger.info("\n" + "─" * 80)
        logger.info("STEP 3: Call Variants [SKIPPED]")
        logger.info("─" * 80)

    # Step 4: Run GWAS workflow
    logger.info("\n" + "─" * 80)
    logger.info("STEP 4: Run GWAS Analysis")
    logger.info("─" * 80)

    # Enable visualization suite (all plot types)
    if args.comprehensive_plots:
        if not config.output:
            config.output = {}
        config.output["comprehensive_plots"] = True
        logger.info("Visualization suite enabled")

    logger.info("Executing full GWAS workflow...")
    workflow_result = execute_gwas_workflow(config)

    # Save final results
    summary_file = config.work_dir / "workflow_summary.json"
    dump_json(workflow_result, summary_file, indent=2)

    logger.info("\n" + "═" * 80)
    logger.info("WORKFLOW COMPLETE")
    logger.info("═" * 80)
    logger.info(f"\nResults saved to: {config.work_dir}")
    logger.info(f"Summary: {summary_file}")

    # Print statistics
    for step in workflow_result.get("steps", []):
        step_name = step.get("step", "unknown")
        result = step.get("result", {})
        status = result.get("status", "unknown")
        logger.info(f"  {step_name}: {status}")

    logger.info("\n✓ Genome-scale GWAS complete!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
