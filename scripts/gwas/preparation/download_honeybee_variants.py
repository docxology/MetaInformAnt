#!/usr/bin/env python3
"""Download real Apis mellifera variant data from public repositories.

This script demonstrates how to access and download real honeybee genomic data
from NCBI SRA and other public repositories.

Key Data Sources:
    - NCBI BioProject PRJNA292680: Scout/recruit behavioral caste variants
    - NCBI SRA: Raw sequencing data for variant calling
    - Public VCF repositories: Pre-called variant sets

Usage:
    python3 scripts/download_honeybee_variants.py

Workflow:
    1. Search for available datasets
    2. Download SRA runs (if SRA Toolkit available)
    3. OR: Download pre-called VCF files (if available)
    4. Align reads and call variants (if starting from FASTQ)

Requirements:
    - SRA Toolkit (fasterq-dump): https://github.com/ncbi/sra-tools
    - NCBI E-utilities (optional): https://www.ncbi.nlm.nih.gov/books/NBK179288/
    - wget or curl for direct downloads
"""

import json
import logging
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core.io import ensure_directory
from metainformant.gwas.download import download_variant_data
from metainformant.gwas.sra_download import (
    check_sra_tools_available,
    download_sra_project,
    search_sra_for_organism,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def print_section(title: str) -> None:
    """Print formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80 + "\n")


def main() -> None:
    """Download real Apis mellifera variant data."""

    print_section("REAL APIS MELLIFERA VARIANT DATA DOWNLOAD")

    # Setup directories
    data_dir = Path("data/variants/amellifera/real")
    ensure_directory(data_dir)

    print(f"📁 Output directory: {data_dir}")

    # Check available tools
    print_section("CHECKING AVAILABLE TOOLS")

    has_sra = check_sra_tools_available()
    print(f"  SRA Toolkit: {'✅ Available' if has_sra else '❌ Not installed'}")

    if not has_sra:
        print(
            """
  ⚠️  SRA Toolkit not found!
  
  To download raw sequencing data from NCBI SRA, install SRA Toolkit:
  
  Ubuntu/Debian:
    sudo apt-get install sra-toolkit
  
  macOS (Homebrew):
    brew install sra-tools
  
  Or download from: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
  
  After installation, configure with: vdb-config --interactive
        """
        )

    # Search for Apis mellifera data
    print_section("SEARCHING FOR APIS MELLIFERA GENOMIC DATA")

    search_result = search_sra_for_organism(
        organism="Apis mellifera",
        strategy="WGS",
        max_results=100,
    )

    print(f"📊 Search status: {search_result['status']}")
    if search_result.get("sra_search_url"):
        print(f"🔗 Manual search: {search_result['sra_search_url']}")
    if search_result.get("note"):
        print(f"📝 Note: {search_result['note']}")

    # Key datasets for Apis mellifera
    print_section("KEY APIS MELLIFERA DATASETS")

    datasets = [
        {
            "name": "Scout/Recruit Behavioral Caste Variants",
            "bioproject": "PRJNA292680",
            "description": "Genomic variants associated with scout and recruit behavioral castes",
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/292680",
            "runs_example": ["SRR2096937", "SRR2096938", "SRR2096939"],
        },
        {
            "name": "Honey Bee Genome Sequencing",
            "bioproject": "PRJNA13343",
            "description": "Reference genome sequencing project",
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/13343",
        },
        {
            "name": "Apis mellifera Population Genomics",
            "bioproject": "PRJNA392242",
            "description": "Population genomics across multiple subspecies",
            "url": "https://www.ncbi.nlm.nih.gov/bioproject/392242",
        },
    ]

    for i, ds in enumerate(datasets, 1):
        print(f"{i}. {ds['name']}")
        print(f"   BioProject: {ds['bioproject']}")
        print(f"   Description: {ds['description']}")
        print(f"   URL: {ds['url']}")
        if "runs_example" in ds:
            print(f"   Example runs: {', '.join(ds['runs_example'][:3])}")
        print()

    # Demonstrate BioProject download (instructions)
    print_section("DOWNLOADING FROM BIOPROJECT (INSTRUCTIONS)")

    bioproject = "PRJNA292680"
    project_result = download_sra_project(
        bioproject=bioproject,
        dest_dir=data_dir / bioproject,
        max_runs=5,
    )

    print(f"📦 BioProject: {bioproject}")
    print(f"Status: {project_result['status']}")
    print(f"\n{project_result['message']}")

    if "instructions" in project_result:
        print("\nSteps to download:")
        for instruction in project_result["instructions"]:
            print(f"  {instruction}")

    # Public VCF repositories
    print_section("PUBLIC VCF REPOSITORIES")

    print(
        """
  For pre-called variant data, check these resources:
  
  1. European Variation Archive (EVA)
     https://www.ebi.ac.uk/eva/
     Search for "Apis mellifera" to find submitted variant studies
  
  2. NCBI dbSNP
     https://www.ncbi.nlm.nih.gov/snp/
     Limited data for non-model organisms, but worth checking
  
  3. Zenodo / FigShare
     Search for "Apis mellifera VCF" or "honeybee variants"
     Researchers often deposit supplementary data here
  
  4. Research Group Repositories
     - Honey Bee Genome Consortium
     - Individual lab websites (check recent papers)
    """
    )

    # Try downloading from a known public URL (example)
    print_section("EXAMPLE: DOWNLOAD FROM DIRECT URL")

    print(
        """
  If you have a direct URL to a VCF file, you can download it:
  
  Example (hypothetical URL):
    from metainformant.gwas.download import download_variant_data
    
    result = download_variant_data(
        source="custom",
        url="https://example.org/amellifera_variants.vcf.gz",
        dest_dir="data/variants/amellifera/real",
    )
    
  This will download the file using wget or curl.
    """
    )

    # Real workflow example
    print_section("REAL WORKFLOW: SRA TO VARIANTS")

    print(
        """
  Complete workflow to go from SRA data to variants:
  
  1. DOWNLOAD SRA DATA
     fasterq-dump SRR2096937 -O data/raw/sra/ -e 8 --split-files
  
  2. CHECK QUALITY
     fastqc data/raw/sra/SRR2096937_*.fastq -o data/qc/
  
  3. ALIGN TO REFERENCE
     bwa mem -t 8 genome.fna \\
       data/raw/sra/SRR2096937_1.fastq \\
       data/raw/sra/SRR2096937_2.fastq | \\
     samtools sort -@ 8 -o data/aligned/SRR2096937.bam
     samtools index data/aligned/SRR2096937.bam
  
  4. CALL VARIANTS
     bcftools mpileup -f genome.fna data/aligned/SRR2096937.bam | \\
     bcftools call -mv -Oz -o data/variants/SRR2096937.vcf.gz
  
  5. RUN GWAS
     Use the METAINFORMANT GWAS module with the generated VCF
    """
    )

    # Summary and next steps
    print_section("SUMMARY AND NEXT STEPS")

    print(
        """
  ✅ CHECKED: Available tools and data sources
  📚 PROVIDED: Links to key Apis mellifera datasets
  📝 DOCUMENTED: Complete workflow from SRA to variants
  
  NEXT STEPS:
  
  1. Install SRA Toolkit if not already installed
     https://github.com/ncbi/sra-tools
  
  2. Visit NCBI SRA Run Selector for PRJNA292680:
     https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA292680
  
  3. Select runs and download metadata (CSV with accession list)
  
  4. Download SRA runs:
     while read acc; do fasterq-dump $acc -O data/raw/; done < accessions.txt
  
  5. Align reads and call variants (see workflow above)
  
  6. Run GWAS pipeline:
     python3 -c "from metainformant.gwas import execute_gwas_workflow, load_gwas_config; \\
                 config = load_gwas_config('config/gwas/gwas_amellifera.yaml'); \\
                 execute_gwas_workflow(config)"
  
  ALTERNATIVE: Look for pre-called VCF files in:
    - European Variation Archive (EVA)
    - Supplementary data from published papers
    - Zenodo/FigShare repositories
    """
    )

    print("\n" + "=" * 80)
    print("  For questions or issues, see docs/gwas/data_acquisition.md")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
