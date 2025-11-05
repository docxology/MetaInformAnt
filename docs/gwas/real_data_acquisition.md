# Accessing Real Apis mellifera Variant Data

## Overview

This document provides comprehensive guidance on obtaining and processing **real** *Apis mellifera* genomic variant data from public repositories, not synthetic test data.

## Quick Start: Key Data Sources

### 1. NCBI BioProjects with Honeybee Genomic Data

| BioProject | Description | URL | Data Type |
|------------|-------------|-----|-----------|
| **PRJNA292680** | Scout/recruit behavioral caste variants | [Link](https://www.ncbi.nlm.nih.gov/bioproject/292680) | WGS (paired-end) |
| PRJNA13343 | Reference genome sequencing | [Link](https://www.ncbi.nlm.nih.gov/bioproject/13343) | Reference |
| PRJNA392242 | Population genomics | [Link](https://www.ncbi.nlm.nih.gov/bioproject/392242) | WGS |
| PRJNA429464 | Microbiome (16S) | [Link](https://db.cngb.org/data_resources/project/PRJNA429464) | 16S V4 |

### 2. European Variation Archive (EVA)

Search for submitted Apis mellifera studies:
- **URL**: https://www.ebi.ac.uk/eva/
- **Search**: "Apis mellifera" or "honeybee"
- **Format**: Pre-called VCF files

### 3. Research Publications with Supplementary Data

Many recent honeybee GWAS papers include supplementary VCF files:
- Check journal websites (Nature, PLOS, BMC, etc.)
- Look in Supplementary Materials sections
- Often hosted on Zenodo, FigShare, or Dryad

## Method 1: Download from NCBI SRA (Raw Sequencing Data)

### Prerequisites

Install SRA Toolkit:

```bash
# Ubuntu/Debian
sudo apt-get install sra-toolkit

# macOS (Homebrew)
brew install sra-tools

# Or download from: https://github.com/ncbi/sra-tools
```

Configure SRA Toolkit:
```bash
vdb-config --interactive
```

### Step-by-Step Workflow

#### 1. Find SRA Runs

Visit the NCBI SRA Run Selector:
- https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA292680

Select runs of interest and download the metadata CSV or accession list.

**Example Accessions from PRJNA292680:**
- `SRR2096937` - Scout caste sample 1
- `SRR2096938` - Scout caste sample 2
- `SRR2096939` - Recruit caste sample 1
- `SRR2096940` - Recruit caste sample 2
- Additional runs available (~15 total)

#### 2. Download SRA Data

```bash
# Download single run
fasterq-dump SRR2096937 -O data/raw/sra/ -e 8 --split-files -p

# Download multiple runs
for acc in SRR2096937 SRR2096938 SRR2096939; do
    fasterq-dump $acc -O data/raw/sra/ -e 8 --split-files -p
done

# Or use parallel download
cat accessions.txt | parallel -j 4 "fasterq-dump {} -O data/raw/sra/ -e 8 --split-files"
```

**Expected Output:**
- `SRR2096937_1.fastq` (read 1, ~1-2 GB)
- `SRR2096937_2.fastq` (read 2, ~1-2 GB)

#### 3. Quality Control

```bash
# Install FastQC if not available
sudo apt-get install fastqc

# Run quality control
fastqc data/raw/sra/SRR2096937_*.fastq -o data/qc/ -t 4

# View reports
firefox data/qc/SRR2096937_1_fastqc.html
```

#### 4. Align to Reference Genome

```bash
# Index reference genome (once)
bwa index output/gwas/amellifera/genome/genomic.fna

# Align reads
bwa mem -t 8 -R '@RG\tID:SRR2096937\tSM:SRR2096937\tPL:ILLUMINA' \
    output/gwas/amellifera/genome/genomic.fna \
    data/raw/sra/SRR2096937_1.fastq \
    data/raw/sra/SRR2096937_2.fastq | \
samtools sort -@ 8 -o data/aligned/SRR2096937.bam

# Index BAM file
samtools index data/aligned/SRR2096937.bam

# Check alignment stats
samtools flagstat data/aligned/SRR2096937.bam
```

#### 5. Call Variants

**Option A: bcftools**
```bash
# Call variants for single sample
bcftools mpileup -f output/gwas/amellifera/genome/genomic.fna \
    data/aligned/SRR2096937.bam | \
bcftools call -mv -Oz -o data/variants/amellifera/real/SRR2096937.vcf.gz

# Index VCF
tabix -p vcf data/variants/amellifera/real/SRR2096937.vcf.gz
```

**Option B: GATK HaplotypeCaller** (more accurate)
```bash
# Call variants
gatk HaplotypeCaller \
    -R output/gwas/amellifera/genome/genomic.fna \
    -I data/aligned/SRR2096937.bam \
    -O data/variants/amellifera/real/SRR2096937.vcf.gz \
    --native-pair-hmm-threads 8
```

**For Multiple Samples:**
```bash
# Call variants jointly for better accuracy
bcftools mpileup -f output/gwas/amellifera/genome/genomic.fna \
    data/aligned/*.bam | \
bcftools call -mv -Oz -o data/variants/amellifera/real/cohort.vcf.gz
```

#### 6. Run GWAS with Real Data

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Update config to point to real VCF
config = load_gwas_config("config/gwas/gwas_amellifera.yaml")
config.variants["vcf_files"] = ["data/variants/amellifera/real/cohort.vcf.gz"]

# Execute workflow
result = execute_gwas_workflow(config)
```

---

## Method 2: Download Pre-Called VCF Files

### From European Variation Archive (EVA)

1. Visit https://www.ebi.ac.uk/eva/
2. Search for "Apis mellifera"
3. Browse submitted studies
4. Download VCF files directly

### From Research Publications

Example: Check supplementary data from recent papers:

- "Genomic signatures of selection in honeybees" (Nature, 2023)
- "Population genomics of Apis mellifera" (PLOS Genetics, 2022)

**Typical locations:**
- Journal supplementary files
- Zenodo: https://zenodo.org/ (search "Apis mellifera VCF")
- Dryad: https://datadryad.org/
- FigShare: https://figshare.com/

### Download Example

```bash
# Example: Download from direct URL (hypothetical)
wget -O data/variants/amellifera/real/published_variants.vcf.gz \
    "https://zenodo.org/record/XXXXX/files/amellifera_variants.vcf.gz"

# Verify and index
gunzip -c data/variants/amellifera/real/published_variants.vcf.gz | head -20
tabix -p vcf data/variants/amellifera/real/published_variants.vcf.gz
```

---

## Method 3: Use METAINFORMANT Download Functions

### Python API

```python
from metainformant.gwas.download import download_variant_data
from metainformant.gwas.sra_download import download_sra_run

# Download from custom URL
result = download_variant_data(
    source="custom",
    url="https://example.org/amellifera_variants.vcf.gz",
    dest_dir="data/variants/amellifera/real",
)

# Download SRA run
result = download_sra_run(
    sra_accession="SRR2096937",
    dest_dir="data/raw/sra",
    threads=8,
)
```

### Command Line

```bash
# Use the provided download script
python3 scripts/download_real_honeybee_variants.py
```

---

## Complete Workflow Example

### Scenario: GWAS on Scout vs Recruit Behavior

**Data Source:** PRJNA292680 (Scout/recruit behavioral caste variants)

#### 1. Download Samples

```bash
# Create directories
mkdir -p data/raw/sra data/aligned data/variants/amellifera/real

# Download scout samples
for acc in SRR2096937 SRR2096938 SRR2096941 SRR2096944 SRR2096947; do
    fasterq-dump $acc -O data/raw/sra/ -e 8 --split-files
done

# Download recruit samples
for acc in SRR2096939 SRR2096940 SRR2096943 SRR2096946 SRR2096949; do
    fasterq-dump $acc -O data/raw/sra/ -e 8 --split-files
done
```

#### 2. Process All Samples

```bash
# Align all samples
for fastq1 in data/raw/sra/*_1.fastq; do
    base=$(basename $fastq1 _1.fastq)
    fastq2="data/raw/sra/${base}_2.fastq"
    
    echo "Processing $base..."
    bwa mem -t 8 -R "@RG\tID:$base\tSM:$base\tPL:ILLUMINA" \
        output/gwas/amellifera/genome/genomic.fna \
        $fastq1 $fastq2 | \
    samtools sort -@ 8 -o data/aligned/${base}.bam
    
    samtools index data/aligned/${base}.bam
done
```

#### 3. Call Variants (Joint Calling)

```bash
# Call variants for all samples
bcftools mpileup -f output/gwas/amellifera/genome/genomic.fna \
    data/aligned/*.bam | \
bcftools call -mv -Oz -o data/variants/amellifera/real/scout_recruit_cohort.vcf.gz

# Index
tabix -p vcf data/variants/amellifera/real/scout_recruit_cohort.vcf.gz

# Check variant count
bcftools stats data/variants/amellifera/real/scout_recruit_cohort.vcf.gz | grep "number of SNPs:"
```

#### 4. Create Phenotype File

```bash
# Create phenotype file (behavior: 0 = recruit, 1 = scout)
cat > data/phenotypes/amellifera/scout_recruit_behavior.tsv << 'EOF'
sample_id	behavior	colony_size	age_days
SRR2096937	1	45000	21
SRR2096938	1	42000	19
SRR2096939	0	48000	18
SRR2096940	0	45000	20
SRR2096941	1	44000	22
SRR2096943	0	46000	19
SRR2096944	1	43000	21
SRR2096946	0	47000	18
SRR2096947	1	44500	20
SRR2096949	0	46500	19
EOF
```

#### 5. Run GWAS

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load and update config
config = load_gwas_config("config/gwas/gwas_amellifera.yaml")
config.variants["vcf_files"] = ["data/variants/amellifera/real/scout_recruit_cohort.vcf.gz"]
config.samples["phenotype_file"] = "data/phenotypes/amellifera/scout_recruit_behavior.tsv"
config.association["trait"] = "behavior"
config.association["model"] = "logistic"  # Binary trait

# Execute
result = execute_gwas_workflow(config)

# View results
import pandas as pd
results = pd.read_csv("output/gwas/amellifera/work/results/association_results.tsv", sep="\t")
significant = results[results['bonferroni'] < 0.05]
print(f"Found {len(significant)} significant variants")
```

---

## Troubleshooting

### SRA Download Issues

**Error:** "Failed to call external services"
- **Solution:** Configure SRA Toolkit: `vdb-config --interactive`
- Enable remote access and prefetch options

**Error:** Slow downloads
- **Solution:** Use `prefetch` first: `prefetch SRR2096937` then `fasterq-dump SRR2096937`

### Alignment Issues

**Error:** Reference not indexed
- **Solution:** Run `bwa index genome.fna` and `samtools faidx genome.fna`

**Low alignment rate (<80%)**
- Check that reference genome matches the species/assembly
- Verify FASTQ quality with FastQC

### Variant Calling Issues

**Too few variants**
- Increase sequencing depth (download more samples)
- Adjust bcftools call parameters: `-v` for variants only

**Too many low-quality variants**
- Filter by quality: `bcftools view -i 'QUAL>=30' input.vcf.gz`
- Use GATK VQSR for variant quality score recalibration

---

## Data Size Estimates

| Data Type | Per Sample | 10 Samples | 50 Samples |
|-----------|------------|------------|------------|
| Raw FASTQ | 2-4 GB | 20-40 GB | 100-200 GB |
| Aligned BAM | 1-2 GB | 10-20 GB | 50-100 GB |
| VCF (cohort) | - | 100-500 MB | 500 MB-2 GB |

**Storage Requirements:**
- Minimum: 50 GB for 10 samples (with cleanup)
- Recommended: 500 GB for 50+ sample GWAS
- Large-scale (100+): 1-2 TB

---

## Additional Resources

### Tools

- **SRA Toolkit**: https://github.com/ncbi/sra-tools
- **BWA**: http://bio-bwa.sourceforge.net/
- **SAMtools/BCFtools**: http://www.htslib.org/
- **GATK**: https://gatk.broadinstitute.org/
- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

### Databases

- **NCBI SRA**: https://www.ncbi.nlm.nih.gov/sra
- **EVA**: https://www.ebi.ac.uk/eva/
- **BeeBiome**: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-025-06229-7
- **Hymenoptera Genome Database**: https://hymenopteragenome.org/

### Publications

- Wallberg et al. (2014). "A worldwide survey of genome sequence variation provides insight into the evolutionary history of the honeybee Apis mellifera." *Nature Genetics* 46:1081-1088.
- Harpur et al. (2014). "Population genomics of the honey bee reveals strong signatures of positive selection on worker traits." *PNAS* 111(7):2614-2619.
- Recent papers on BioRxiv/medRxiv: Search "Apis mellifera genomics"

---

## Summary

**Real Data Options:**

1. ✅ **SRA Download** (Most flexible, requires processing)
   - Download raw FASTQ with `fasterq-dump`
   - Align with BWA, call variants with bcftools/GATK
   - Full control over quality and parameters

2. ✅ **Pre-Called VCF** (Fastest, limited availability)
   - Check EVA, Zenodo, publication supplements
   - Direct download with wget/curl
   - Ready for immediate analysis

3. ✅ **Custom URL** (Moderate, requires known source)
   - Use `download_variant_data(source="custom", url=...)`
   - Download from research group servers
   - Verify data quality

**Recommended Workflow:**
- For learning: Use synthetic data (fast)
- For pilot study: Download 5-10 samples from PRJNA292680
- For publication: Download full cohort (50+ samples) or use pre-called VCF

All methods integrate seamlessly with METAINFORMANT GWAS pipeline!






