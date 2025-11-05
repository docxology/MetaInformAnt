# Installing Bioinformatics Tools for Genome-Scale GWAS

## Required Tools

To process genome-scale data from SRA to variants, you need:

1. **BWA** (Burrows-Wheeler Aligner) - Read alignment
2. **SAMtools** - BAM file manipulation
3. **bcftools** - Variant calling and VCF processing

## Installation Instructions

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install -y bwa samtools bcftools
```

### macOS (Homebrew)

```bash
brew install bwa samtools bcftools
```

### Conda (Cross-platform)

```bash
conda install -c bioconda bwa samtools bcftools
```

### Verify Installation

```bash
bwa 2>&1 | head -3
samtools --version
bcftools --version
```

Expected output should show version information for each tool.

## Tool Versions Tested

- BWA: v0.7.17 or higher
- SAMtools: v1.15 or higher
- bcftools: v1.15 or higher

## Disk Space Requirements

- Raw FASTQ: ~10-15 GB (3-5 samples)
- Aligned BAM: ~4-8 GB
- Variant VCF: ~500 MB - 1 GB
- **Total**: ~20-30 GB for complete workflow
- **Recommended**: 50-100 GB free space

## Memory Requirements

- Alignment (BWA): 4-8 GB RAM
- Variant calling (bcftools): 2-4 GB RAM
- GWAS analysis: 2-8 GB RAM

## Alternative: Docker Container

If you prefer containerization:

```bash
docker pull biocontainers/bwa:v0.7.17-3-deb_cv1
docker pull biocontainers/samtools:v1.9-4-deb_cv1
docker pull biocontainers/bcftools:v1.9-1-deb_cv1
```

## Troubleshooting

### "Command not found" errors
Ensure tools are in your PATH:
```bash
export PATH=$PATH:/usr/local/bin
```

### Permission denied
Use sudo for installation or install to user directory with conda

### Out of memory
Reduce thread count in scripts (use -t 4 instead of -t 8)

## Next Steps

After installation, verify with:
```bash
bash scripts/check_tools.sh
```

Then proceed with genome-scale workflow:
```bash
bash scripts/download_genome_scale_data.sh  # Already running
bash scripts/align_samples.sh
bash scripts/call_variants.sh
```






