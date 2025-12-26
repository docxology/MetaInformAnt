# Genome Setup Guide

**User Guide** - Step-by-step instructions for setting up genomes and kallisto indexes.

> **For technical API reference**: See [genome_preparation.md](genome_preparation.md)  
> **For command reference**: See [commands.md](commands.md)  
> **For API functions**: See [API Reference](../API.md#genome-preparation-functions)

This guide provides step-by-step instructions for setting up genomes and kallisto indexes for all species in the METAINFORMANT RNA-seq workflow.

## Overview

The genome setup process consists of three main steps:

1. **Download Genomes**: Download NCBI genome packages for each species
2. **Prepare Transcriptomes**: Extract RNA FASTA files from downloaded genomes
3. **Build Indexes**: Build kallisto indexes from transcriptome FASTA files

## Prerequisites

- Python 3.11+
- METAINFORMANT installed or accessible
- `kallisto` installed (for index building): `conda install -c bioconda kallisto` or `sudo apt-get install -y kallisto`
- NCBI datasets CLI (optional, but recommended): `conda install -c bioconda ncbi-datasets-cli`

**Note**: For Python packages, use `uv pip install` (primary method). System tools like kallisto can use conda or system package managers.
- Sufficient disk space (genomes: 100MB-1GB each, indexes: 20-50MB each)

## Quick Start

### Automated Complete Setup

The fastest way to set up all species:

```bash
# Run the master orchestrator
python3 scripts/rna/orchestrate_genome_setup.py
```

This will:
- Verify current status
- Download missing genomes
- Prepare transcriptomes
- Build kallisto indexes
- Verify final status

### Step-by-Step Setup

For more control, run each step individually:

```bash
# 1. Check current status
python3 scripts/rna/verify_genomes_and_indexes.py

# 2. Download missing genomes
python3 scripts/rna/download_missing_genomes.py

# 3. Prepare transcriptomes
python3 scripts/rna/prepare_transcriptomes.py

# 4. Build kallisto indexes
python3 scripts/rna/build_kallisto_indexes.py

# 5. Verify final status
python3 scripts/rna/verify_genomes_and_indexes.py
```

## Detailed Instructions

### Step 1: Verification

Before starting, check what's already done:

```bash
python3 scripts/rna/verify_genomes_and_indexes.py
```

**Output**: Status report showing:
- Which genomes are downloaded
- Which transcriptomes are prepared
- Which indexes are built
- Recommended next steps

**Example Output**:
```
================================================================================
VERIFICATION SUMMARY
================================================================================
Total species: 24

Genome Downloads:
  ✓ Downloaded: 9/24 (37.5%)
  ✗ Missing: 15/24

RNA FASTA Files (of downloaded genomes):
  ✓ Found: 3/9 (33.3%)
  ✗ Missing: 6/9

Kallisto Indexes (of prepared transcriptomes):
  ✓ Built: 4/3 (133.3%)
  ✗ Missing: 0/3

================================================================================

Next Steps:
  → Download 15 missing genomes:
     python3 scripts/rna/download_missing_genomes.py
  → Prepare 6 transcriptomes:
     python3 scripts/rna/prepare_transcriptomes.py
```

### Step 2: Download Genomes

Download genomes for species that are missing them:

```bash
# Download all missing genomes
python3 scripts/rna/download_missing_genomes.py

# Download for specific species
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus

# Dry run (see what would be downloaded)
python3 scripts/rna/download_missing_genomes.py --dry-run
```

**What it does**:
- Scans all config files in `config/amalgkit/`
- Identifies species with missing genome downloads
- Downloads using best-effort methods (CLI → API → FTP)
- Extracts genome packages
- Logs progress and results

**Output Files**:
- `output/amalgkit/{species}/genome/download_record.json` - Download metadata
- `output/genome_download_results.json` - Summary report

**Time**: ~5-30 minutes per genome depending on size and network

### Step 3: Prepare Transcriptomes

Extract RNA FASTA files from downloaded genomes:

```bash
# Prepare all transcriptomes
python3 scripts/rna/prepare_transcriptomes.py

# Prepare for specific species
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus

# Dry run
python3 scripts/rna/prepare_transcriptomes.py --dry-run
```

**What it does**:
- Finds RNA FASTA files in extracted genome packages
- Decompresses if needed (handles `.gz` files)
- Copies to `work_dir/fasta/{Species_Name}_rna.fasta`
- Skips already-prepared transcriptomes

**Output Files**:
- `output/amalgkit/{species}/work/fasta/{Species_Name}_rna.fasta` - Prepared transcriptome
- `output/transcriptome_preparation_results.json` - Summary report

**Time**: ~1-5 minutes per species

### Step 4: Build Kallisto Indexes

Build kallisto indexes from prepared transcriptomes:

```bash
# Build all indexes
python3 scripts/rna/build_kallisto_indexes.py

# Build for specific species
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus

# Custom k-mer size (for short reads <100bp)
python3 scripts/rna/build_kallisto_indexes.py --kmer-size 23

# Dry run
python3 scripts/rna/build_kallisto_indexes.py --dry-run
```

**What it does**:
- Checks if kallisto is available
- Finds prepared transcriptome FASTA files
- Builds kallisto indexes with specified k-mer size
- Skips already-built indexes

**Output Files**:
- `output/amalgkit/{species}/work/index/{Species_Name}_transcripts.idx` - Kallisto index
- `output/kallisto_index_build_results.json` - Summary report

**Time**: ~1-10 minutes per species depending on transcriptome size

**k-mer Size Guidelines**:
- `31` (default): Standard reads ≥100bp
- `23`: Short reads <100bp

### Step 5: Final Verification

Verify that everything is set up correctly:

```bash
python3 scripts/rna/verify_genomes_and_indexes.py
```

**Expected Output** (when complete):
```
================================================================================
VERIFICATION SUMMARY
================================================================================
Total species: 24

Genome Downloads:
  ✓ Downloaded: 24/24 (100.0%)
  ✗ Missing: 0/24

RNA FASTA Files (of downloaded genomes):
  ✓ Found: 24/24 (100.0%)
  ✗ Missing: 0/24

Kallisto Indexes (of prepared transcriptomes):
  ✓ Built: 24/24 (100.0%)
  ✗ Missing: 0/24

================================================================================

Next Steps:
  ✓ All species are fully prepared!
```

## Single Species Setup

For setting up a single species:

```bash
# Complete setup for one species
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus

# Or step-by-step:
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus
```

## Troubleshooting

### Issue: Genome Download Fails

**Symptoms**: Download script reports errors

**Solutions**:
1. Check network connectivity
2. Verify NCBI datasets CLI is installed: `datasets --version`
3. Check disk space: `df -h`
4. Try downloading single species to isolate the issue
5. Check download logs: `output/amalgkit/{species}/genome/download_record.json`

### Issue: RNA FASTA Not Found

**Symptoms**: Prepare script reports "RNA FASTA not found"

**Solutions**:
1. Verify genome download completed: Check `download_record.json`
2. Ensure `include: ["rna"]` in genome config
3. Manually search for RNA FASTA:
   ```bash
   find output/amalgkit/{species}/genome -name "rna.fna*"
   ```
4. Re-download genome if RNA FASTA is missing

### Issue: Kallisto Index Build Fails

**Symptoms**: Build script reports errors

**Solutions**:
1. Verify kallisto is installed: `which kallisto`
2. Check transcriptome FASTA exists:
   ```bash
   ls -lh output/amalgkit/{species}/work/fasta/*.fasta
   ```
3. Verify FASTA file is valid:
   ```bash
   head output/amalgkit/{species}/work/fasta/*.fasta
   ```
4. Check disk space
5. Try building with different k-mer size

### Issue: Script Can't Find Config Files

**Symptoms**: "Config directory not found" error

**Solutions**:
1. Run from repository root
2. Specify config directory explicitly:
   ```bash
   python3 scripts/rna/verify_genomes_and_indexes.py --config-dir config/amalgkit
   ```
3. Check that config files exist:
   ```bash
   ls config/amalgkit/amalgkit_*.yaml
   ```

## Best Practices

### 1. Incremental Setup

Set up genomes incrementally as needed:

```bash
# Check status first
python3 scripts/rna/verify_genomes_and_indexes.py

# Download only what's missing
python3 scripts/rna/download_missing_genomes.py

# Prepare only what's ready
python3 scripts/rna/prepare_transcriptomes.py

# Build only what's needed
python3 scripts/rna/build_kallisto_indexes.py
```

### 2. Dry Run First

Always test with `--dry-run` before actual operations:

```bash
python3 scripts/rna/download_missing_genomes.py --dry-run
python3 scripts/rna/prepare_transcriptomes.py --dry-run
python3 scripts/rna/build_kallisto_indexes.py --dry-run
```

### 3. Verify Regularly

Check status regularly to track progress:

```bash
python3 scripts/rna/verify_genomes_and_indexes.py
```

### 4. Monitor Disk Space

Genome downloads can be large:

```bash
# Check space usage
du -sh output/amalgkit/*/genome
du -sh output/amalgkit/*/work/fasta
du -sh output/amalgkit/*/work/index
```

### 5. Keep Results Files

All scripts generate JSON result files. Keep these for:
- Tracking progress
- Debugging issues
- Reporting status

## Script Reference

| Script | Purpose | Key Options |
|--------|---------|-------------|
| `verify_genomes_and_indexes.py` | Check status | `--species`, `--output` |
| `download_missing_genomes.py` | Download genomes | `--species`, `--dry-run` |
| `prepare_transcriptomes.py` | Prepare FASTA | `--species`, `--dry-run` |
| `build_kallisto_indexes.py` | Build indexes | `--species`, `--kmer-size`, `--dry-run` |
| `orchestrate_genome_setup.py` | Complete pipeline | `--species`, `--skip-*`, `--dry-run` |

## File Locations

### Genome Files
- Download location: `output/amalgkit/{species}/genome/`
- RNA FASTA (in package): `output/amalgkit/{species}/genome/ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna`

### Prepared Transcriptomes
- Location: `output/amalgkit/{species}/work/fasta/{Species_Name}_rna.fasta`

### Kallisto Indexes
- Location: `output/amalgkit/{species}/work/index/{Species_Name}_transcripts.idx`

### Results Files
- Verification: `output/genome_index_status.json`
- Downloads: `output/genome_download_results.json`
- Preparation: `output/transcriptome_preparation_results.json`
- Index building: `output/kallisto_index_build_results.json`

## See Also

### Technical Reference
- **[genome_preparation.md](genome_preparation.md)** - Technical API reference
- **[commands.md](commands.md)** - Command reference for all scripts

### Documentation
- **[API Reference](../API.md#genome-preparation-functions)** - Complete function documentation
- **[Function Index](FUNCTIONS.md)** - Quick function lookup
- **[Quantification Step](steps/quant.md)** - Using indexes for quantification
- **[Pipeline Overview](amalgkit.md)** - Complete workflow documentation
- **[Main Index](../README.md)** - RNA domain master index

