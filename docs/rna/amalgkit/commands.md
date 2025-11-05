# Genome Setup Commands

This document lists all commands to run the genome setup scripts sequentially.

## Quick Reference

### Option 1: Run All Steps Sequentially (Recommended)

```bash
# Run the complete pipeline script
bash scripts/rna/run_genome_setup.sh
```

### Option 2: Use Master Orchestrator

```bash
# Run everything in one command
python3 scripts/rna/orchestrate_genome_setup.py
```

### Option 3: Run Steps Individually

Run each step manually for more control:

```bash
# Step 1: Check current status
python3 scripts/rna/verify_genomes_and_indexes.py

# Step 2: Download missing genomes
python3 scripts/rna/download_missing_genomes.py

# Step 3: Prepare transcriptomes
python3 scripts/rna/prepare_transcriptomes.py

# Step 4: Build kallisto indexes
python3 scripts/rna/build_kallisto_indexes.py

# Step 5: Verify final status
python3 scripts/rna/verify_genomes_and_indexes.py
```

## Detailed Command Reference

### Step 1: Verification

**Check status of all species:**

```bash
python3 scripts/rna/verify_genomes_and_indexes.py
```

**Check specific species:**

```bash
python3 scripts/rna/verify_genomes_and_indexes.py --species camponotus_floridanus
```

**Custom output location:**

```bash
python3 scripts/rna/verify_genomes_and_indexes.py --output output/my_status.json
```

### Step 2: Download Genomes

**Download all missing genomes:**

```bash
python3 scripts/rna/download_missing_genomes.py
```

**Download specific species:**

```bash
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus
```

**Dry run (see what would be downloaded):**

```bash
python3 scripts/rna/download_missing_genomes.py --dry-run
```

### Step 3: Prepare Transcriptomes

**Prepare all transcriptomes:**

```bash
python3 scripts/rna/prepare_transcriptomes.py
```

**Prepare specific species:**

```bash
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
```

**Dry run:**

```bash
python3 scripts/rna/prepare_transcriptomes.py --dry-run
```

### Step 4: Build Kallisto Indexes

**Build all indexes:**

```bash
python3 scripts/rna/build_kallisto_indexes.py
```

**Build specific species:**

```bash
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus
```

**Custom k-mer size (for short reads):**

```bash
python3 scripts/rna/build_kallisto_indexes.py --kmer-size 23
```

**Dry run:**

```bash
python3 scripts/rna/build_kallisto_indexes.py --dry-run
```

### Step 5: Master Orchestrator

**Run complete pipeline:**

```bash
python3 scripts/rna/orchestrate_genome_setup.py
```

**Run for specific species:**

```bash
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus
```

**Skip specific steps:**

```bash
# Skip download step (if genomes already downloaded)
python3 scripts/rna/orchestrate_genome_setup.py --skip-download

# Skip preparation step
python3 scripts/rna/orchestrate_genome_setup.py --skip-prepare

# Skip index building
python3 scripts/rna/orchestrate_genome_setup.py --skip-build

# Skip initial verification
python3 scripts/rna/orchestrate_genome_setup.py --skip-verify-initial

# Skip final verification
python3 scripts/rna/orchestrate_genome_setup.py --skip-verify-final
```

**Dry run:**

```bash
python3 scripts/rna/orchestrate_genome_setup.py --dry-run
```

## Example Workflow

### Complete Setup for All Species

```bash
# Option A: Automated pipeline script
bash scripts/rna/run_genome_setup.sh

# Option B: Master orchestrator
python3 scripts/rna/orchestrate_genome_setup.py

# Option C: Manual step-by-step
python3 scripts/rna/verify_genomes_and_indexes.py
python3 scripts/rna/download_missing_genomes.py
python3 scripts/rna/prepare_transcriptomes.py
python3 scripts/rna/build_kallisto_indexes.py
python3 scripts/rna/verify_genomes_and_indexes.py
```

### Setup for Single Species

```bash
# Complete setup for one species
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus

# Or step-by-step:
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus
```

### Incremental Setup

If you only need to do specific steps:

```bash
# Check what's missing
python3 scripts/rna/verify_genomes_and_indexes.py

# Download only missing genomes
python3 scripts/rna/download_missing_genomes.py

# Prepare only newly downloaded genomes
python3 scripts/rna/prepare_transcriptomes.py

# Build only newly prepared transcriptomes
python3 scripts/rna/build_kallisto_indexes.py
```

## Output Files

All scripts write JSON result files to `output/`:

- `output/genome_index_status.json` - Verification results
- `output/genome_download_results.json` - Download results
- `output/transcriptome_preparation_results.json` - Preparation results
- `output/kallisto_index_build_results.json` - Index build results

## Troubleshooting

### If a script fails midway

You can resume by running only the failed steps:

```bash
# Check current status
python3 scripts/rna/verify_genomes_and_indexes.py

# Continue from where it left off
python3 scripts/rna/download_missing_genomes.py  # Only downloads missing
python3 scripts/rna/prepare_transcriptomes.py     # Only prepares missing
python3 scripts/rna/build_kallisto_indexes.py     # Only builds missing
```

### Check progress during long downloads

```bash
# Check status while downloads are running (in another terminal)
python3 scripts/rna/verify_genomes_and_indexes.py

# Check download progress files
cat output/amalgkit/*/genome/download.progress.txt
```

## See Also

- [Genome Setup Guide](genome_setup_guide.md) - Detailed setup instructions
- [Genome Preparation](genome_preparation.md) - Technical documentation and API reference
- [Quantification Step](steps/quant.md) - Using kallisto indexes for quantification
- [Workflow Guide](amalgkit.md) - Complete workflow documentation
- [Amalgkit README](README.md) - Overview of amalgkit integration
- [RNA Domain README](../README.md) - RNA domain overview

