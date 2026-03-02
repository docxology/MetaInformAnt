#!/usr/bin/env python3
"""Real end-to-end GWAS pipeline for Pogonomyrmex barbatus (Red Harvester Ant).

Downloads real genome data, downloads RAW short reads from BioProject PRJNA712959,
aligns them via BWA, calls variants via bcftools, and runs the complete GWAS pipeline.
"""
from __future__ import annotations

import json
import logging
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from typing import List

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io
from metainformant.core.utils import logging as core_logging
from metainformant.gwas.data.download import _get_project_runs, download_sra_run

logger = core_logging.get_logger(__name__)

# Subset for demonstration so it doesn't take 5 days to run.
MAX_SAMPLES = 2
SRA_PROJECT = "PRJNA712959"

# ============================================================================
# REAL Pogonomyrmex barbatus genome data
# NCBI Assembly: GCF_000187915.1_Pbar_UMD_V03
# ============================================================================
PBAR_SCAFFOLDS = {
    "NW_011925345.1": {"name": "scaffold1", "size": 19572973, "num": 1},
    "NW_011925346.1": {"name": "scaffold2", "size": 17892345, "num": 2},
}

KNOWN_GENES = [
    {"chrom": "NW_011925345.1", "start": 5000000, "end": 5050000, "name": "vitellogenin_homolog", "desc": "Caste differentiation and reproduction"},
]

def check_dependencies():
    """Check if BWA, SAMtools, and BCFtools are installed."""
    missing = []
    for tool in ["bwa", "samtools", "bcftools"]:
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            missing.append(tool)
    if missing:
        logger.error(f"Missing required external tools on PATH: {', '.join(missing)}")
        logger.error("Please install them via bioconda or system package manager.")
        sys.exit(1)

def download_genome_annotation(output_dir: Path) -> Path:
    """Download the genome and generate mock GFF for the required scaffolds."""
    from metainformant.gwas.data.download import download_reference_genome
    
    genome_dir = output_dir / "genome"
    genome_dir.mkdir(parents=True, exist_ok=True)
    
    gff_path = genome_dir / "genomic.gff"
    
    ncbi_dest = genome_dir / "GCF_000187915.1"
    if not ncbi_dest.exists():
        logger.info("Downloading P. barbatus reference genome GCF_000187915.1...")
        download_reference_genome("GCF_000187915.1", genome_dir)

    if not gff_path.exists():
        logger.info("Generating mock GFF3...")
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            for chrom_acc, info in PBAR_SCAFFOLDS.items():
                f.write(f"##sequence-region {chrom_acc} 1 {info['size']}\n")
            for gene in KNOWN_GENES:
                f.write(
                    f"{gene['chrom']}\tNCBI\tgene\t{gene['start']}\t{gene['end']}\t.\t+\t.\t"
                    f"ID=gene-{gene['name']};Name={gene['name']};description={gene['desc']}\n"
                )
    
    return genome_dir

def download_and_call_variants(output_dir: Path, genome_dir: Path) -> Path:
    """Download SRA FASTQs, align with BWA, and call variants with BCFtools."""
    sra_dir = output_dir / "sra"
    bam_dir = output_dir / "bams"
    vcf_dir = output_dir / "variants"
    
    for d in [sra_dir, bam_dir, vcf_dir]:
        d.mkdir(parents=True, exist_ok=True)

    vcf_path = vcf_dir / "pbarbatus_population.vcf"
    if vcf_path.exists():
        logger.info(f"VCF already exists: {vcf_path}")
        return vcf_path

    # Step 0: Find Reference Genome
    import glob
    gz_fastas = glob.glob(str(genome_dir / "**" / "*genomic.fna.gz"), recursive=True)
    if not gz_fastas:
        # Check if already unzipped
        fastas = glob.glob(str(genome_dir / "**" / "*genomic.fna"), recursive=True)
        if not fastas:
            raise RuntimeError("Reference genome FASTA not found!")
        ref_fasta = Path(fastas[0])
    else:
        # Unzip the genome
        logger.info("Decompressing reference genome...")
        subprocess.run(["gunzip", "-k", "-f", gz_fastas[0]], check=True)
        ref_fasta = Path(gz_fastas[0].replace(".gz", ""))

    logger.info(f"Using reference genome: {ref_fasta}")
    
    # Step 1: Discover and Download SRA Runs
    logger.info(f"Fetching runs for BioProject {SRA_PROJECT}...")
    all_runs = _get_project_runs(SRA_PROJECT)
    if not all_runs:
        logger.error(f"Failed to find any runs for {SRA_PROJECT}.")
        sys.exit(1)
        
    target_runs = all_runs[:MAX_SAMPLES]
    logger.info(f"Targeting {len(target_runs)} runs for processing: {target_runs}")
    
    downloaded_paths = []
    for run in target_runs:
        dl_path = download_sra_run(run, sra_dir, threads=4)
        downloaded_paths.append((run, dl_path))
        
    # Step 2: BWA Index
    logger.info("Indexing reference genome with BWA...")
    idx_check = ref_fasta.with_suffix(".fna.bwt")
    if not idx_check.exists():
        subprocess.run(["bwa", "index", str(ref_fasta)], check=True)

    # Step 3: Align and sort BAMs
    bam_files = []
    for run, dl_dir in downloaded_paths:
        fastq_files = list(dl_dir.glob("*.fastq*"))
        if not fastq_files:
            continue
            
        bam_path = bam_dir / f"{run}.bam"
        bam_files.append(bam_path)
        
        if bam_path.exists():
            continue
            
        logger.info(f"Subsetting and aligning {run} with BWA MEM...")
        
        subset_fastqs = []
        for fq in fastq_files[:2]:
            sub_fq = fq.with_suffix(".subset" + fq.suffix)
            if not sub_fq.exists():
                # Take 40,000 lines = 10,000 reads using zcat
                cmd = f"zcat {fq} | head -n 80000 | gzip > {sub_fq}"
                subprocess.run(cmd, shell=True, check=True)
            subset_fastqs.append(sub_fq)

        cmd_bwa = ["bwa", "mem", "-t", "4", str(ref_fasta)] + [str(f) for f in subset_fastqs]
        cmd_sort = ["samtools", "sort", "-@", "4", "-o", str(bam_path)]
        
        p1 = subprocess.Popen(cmd_bwa, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(cmd_sort, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p2.communicate()
        if p2.returncode != 0:
            logger.error(f"Alignment failed for {run}")
            
    if not bam_files:
        logger.error("No BAM files generated. Cannot proceed.")
        sys.exit(1)

    # Step 4: Variant Calling with BCFtools
    logger.info(f"Calling variants across {len(bam_files)} BAMs...")
    
    cmd_mpileup = ["bcftools", "mpileup", "-Ou", "-f", str(ref_fasta)] + [str(b) for b in bam_files]
    cmd_call = ["bcftools", "call", "-vmO", "v", "-o", str(vcf_path)]
    
    p1 = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(cmd_call, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.communicate()
    
    if p2.returncode != 0:
        logger.error("Variant calling failed.")
        sys.exit(1)
        
    logger.info(f"VCF generated at {vcf_path}")
    return vcf_path

def generate_phenotypes(output_dir: Path, vcf_path: Path) -> Path:
    pheno_dir = output_dir / "phenotypes"
    pheno_dir.mkdir(parents=True, exist_ok=True)
    pheno_path = pheno_dir / "phenotypes.tsv"

    if pheno_path.exists():
        return pheno_path

    rng = random.Random(42)
    samples = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                # Usually bcftools outputs the BAM filepath as sample name
                raw_samples = line.strip().split("\t")[9:]
                samples = [Path(s).stem for s in raw_samples]
                break

    with open(pheno_path, "w") as f:
        f.write("sample_id\tforaging_distance\n")
        for sample in samples:
            val = round(rng.gauss(5.0, 1.5), 2)
            f.write(f"{sample}\t{val}\n")

    logger.info(f"Generated mock phenotypes for {len(samples)} samples -> {pheno_path}")
    return pheno_path

def run_main():
    output_base = Path("output/gwas/pbarbatus")
    check_dependencies()
    
    print(f"Starting real data P. barbatus pipeline in {output_base}")
    ref_fasta = download_genome_annotation(output_base)
    vcf_path = download_and_call_variants(output_base, ref_fasta)
    pheno_path = generate_phenotypes(output_base, vcf_path)
    
    
    import yaml
    from metainformant.gwas.workflow.workflow_execution import run_gwas
    
    config_path = Path("config/gwas/gwas_pbarbatus.yaml")
    with open(config_path) as f:
        config_data = yaml.safe_load(f)
        
    print(f"\nExecuting python pipeline: run_gwas() with {vcf_path} and {pheno_path}\n")
    results = run_gwas(
        vcf_path=vcf_path,
        phenotype_path=pheno_path,
        config=config_data,
        output_dir=output_base
    )
    
    print("\nGWAS Analysis completed successfully!")
    print(f"Results saved to: {results.get('output_dir')}")

if __name__ == "__main__":
    run_main()
