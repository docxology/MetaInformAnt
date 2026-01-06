# Genome Preparation and Kallisto Index Building

**Technical API Reference** - Complete function documentation for genome preparation.

> **For step-by-step user guide**: See [genome_setup_guide.md](genome_setup_guide.md)  
> **For command reference**: See [commands.md](commands.md)  
> **For API functions**: See [API Reference](../API.md#genome-preparation-functions)

## Overview

This document describes the genome preparation pipeline for RNA-seq quantification workflows. The process automatically downloads genomes, extracts transcriptome FASTA files, and builds kallisto indexes for reproducible transcript abundance estimation.

**This is the technical reference** - for user-friendly step-by-step instructions, see [genome_setup_guide.md](genome_setup_guide.md).

## Pipeline Flow

```
1. Download Genome Package
   ↓
2. Extract Transcriptome FASTA
   ↓
3. Prepare FASTA for Kallisto
   ↓
4. Build Kallisto Index (optional)
   ↓
5. Ready for Quantification
```

## Configuration

### Genome Configuration

Each species configuration includes a `genome` section:

```yaml
genome:
  accession: GCF_003227725.1  # NCBI assembly accession
  assembly_name: Cflo_v7.5
  dest_dir: output/amalgkit/camponotus_floridanus/genome
  include:
    - genome           # Genomic sequences
    - gff3            # Gene annotations
    - rna             # RNA sequences (required for kallisto)
    - cds             # CDS sequences
    - protein         # Protein sequences
    - seq-report      # Sequence report
  ftp_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/725/...
```

### Quantification Configuration

Enable automatic index building:

```yaml
steps:
  quant:
    # CRITICAL: out_dir should be work_dir (not separate quant_dir)
    # This allows quant to find getfastq output in {out_dir}/getfastq/
    out_dir: output/amalgkit/camponotus_floridanus/work
    build_index: yes  # Automatically build kallisto index
    threads: 12
```

## Automatic Workflow Integration

The genome preparation is automatically integrated into the workflow execution:

1. **Genome Download**: Downloads genome package from NCBI using best-effort methods (CLI → API → FTP)
2. **Transcriptome Extraction**: Finds RNA FASTA in extracted genome directory
3. **FASTA Preparation**: Copies transcriptome to `work_dir/fasta/{Species_Name}_rna.fasta`
4. **Index Building**: Builds kallisto index if `build_index: yes` is set

### Workflow Execution

```python
from metainformant.rna.workflow import load_workflow_config, execute_workflow

config = load_workflow_config("config/amalgkit/camponotus_floridanus.yaml")
return_codes = execute_workflow(config)
```

The workflow automatically:
- Downloads genome if not present
- Prepares transcriptome FASTA
- Builds kallisto index (if enabled)
- Sets `fasta_dir` and `index_dir` in quant parameters

## File Structure

### Genome Directory

```
output/amalgkit/{species}/genome/
├── download_record.json          # Download metadata
├── ncbi_dataset_api.zip          # Downloaded package (optional)
└── ncbi_dataset_api_extracted/   # Extracted files
    └── ncbi_dataset/
        └── data/
            └── {accession}/
                ├── rna.fna       # RNA transcriptome FASTA
                ├── genomic.fna   # Genome sequence
                ├── genomic.gff  # Gene annotations
                └── ...
```

### Work Directory

```
output/amalgkit/{species}/work/
├── fasta/
│   └── {Species_Name}_rna.fasta  # Prepared transcriptome
├── index/
│   └── {Species_Name}_transcripts.idx  # Kallisto index
└── ...
```

## Manual Genome Preparation

You can also prepare genomes manually using the `genome_prep` module:

```python
from metainformant.rna.genome_prep import (
    prepare_transcriptome_for_kallisto,
    build_kallisto_index,
    prepare_genome_for_quantification,
)
from pathlib import Path

genome_dir = Path("output/amalgkit/camponotus_floridanus/genome")
work_dir = Path("output/amalgkit/camponotus_floridanus/work")
species_name = "Camponotus_floridanus"
accession = "GCF_003227725.1"

# Option 1: Complete preparation
result = prepare_genome_for_quantification(
    genome_dir,
    species_name,
    work_dir,
    accession=accession,
    build_index=True,
    kmer_size=31,
)

if result["success"]:
    print(f"FASTA: {result['fasta_path']}")
    print(f"Index: {result['index_path']}")
else:
    print(f"Error: {result['error']}")

# Option 2: Step-by-step
fasta_path = prepare_transcriptome_for_kallisto(
    genome_dir,
    species_name,
    work_dir,
    accession=accession,
)

if fasta_path:
    index_path = work_dir / "index" / f"{species_name}_transcripts.idx"
    build_kallisto_index(fasta_path, index_path, kmer_size=31)
```

## Orchestration Scripts

METAINFORMANT provides thin orchestrator scripts to manage genome setup at scale. These scripts can be run individually or via the master orchestrator.

### Verification Script

**Script**: `scripts/rna/verify_genomes_and_indexes.py`

Check the status of all species:

```bash
# Check all species
python3 scripts/rna/verify_genomes_and_indexes.py

# Check specific species
python3 scripts/rna/verify_genomes_and_indexes.py --species camponotus_floridanus

# Custom output location
python3 scripts/rna/verify_genomes_and_indexes.py --output output/my_status.json
```

**Output**: Detailed status report with percentages and next-step recommendations

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

### Download Missing Genomes

**Script**: `scripts/rna/download_missing_genomes.py`

Download genomes for species that are missing them:

```bash
# Download all missing genomes
python3 scripts/rna/download_missing_genomes.py

# Download for specific species
python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus

# Dry run (see what would be downloaded)
python3 scripts/rna/download_missing_genomes.py --dry-run
```

**Features**:
- Skips already-downloaded genomes
- Uses best-effort download (CLI → API → FTP)
- Logs progress and results
- Writes JSON report with download status

**Output**: `output/genome_download_results.json`

### Prepare Transcriptomes

**Script**: `scripts/rna/prepare_transcriptomes.py`

Extract and prepare RNA FASTA files for species with downloaded genomes:

```bash
# Prepare all transcriptomes
python3 scripts/rna/prepare_transcriptomes.py

# Prepare for specific species
python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus

# Dry run
python3 scripts/rna/prepare_transcriptomes.py --dry-run
```

**Features**:
- Finds RNA FASTA in downloaded genome packages
- Decompresses if needed
- Copies to `work_dir/fasta/{Species_Name}_rna.fasta`
- Skips already-prepared transcriptomes

**Output**: `output/transcriptome_preparation_results.json`

### Build Kallisto Indexes

**Script**: `scripts/rna/build_kallisto_indexes.py`

Build kallisto indexes for species with prepared transcriptomes:

```bash
# Build all indexes
python3 scripts/rna/build_kallisto_indexes.py

# Build for specific species
python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus

# Custom k-mer size (for short reads)
python3 scripts/rna/build_kallisto_indexes.py --kmer-size 23

# Dry run
python3 scripts/rna/build_kallisto_indexes.py --dry-run
```

**Features**:
- Requires kallisto on PATH
- Uses k-mer size 31 by default (use 23 for short reads)
- Skips already-built indexes
- Verifies index creation

**Output**: `output/kallisto_index_build_results.json`

### Master Orchestrator

**Script**: `scripts/rna/orchestrate_genome_setup.py`

Runs the complete pipeline sequentially:

```bash
# Complete setup for all species
python3 scripts/rna/orchestrate_genome_setup.py

# Setup for specific species
python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus

# Skip specific steps
python3 scripts/rna/orchestrate_genome_setup.py --skip-download --skip-prepare

# Dry run
python3 scripts/rna/orchestrate_genome_setup.py --dry-run
```

**Pipeline Steps**:
1. Initial verification (optional: `--skip-verify-initial`)
2. Download missing genomes (optional: `--skip-download`)
3. Prepare transcriptomes (optional: `--skip-prepare`)
4. Build kallisto indexes (optional: `--skip-build`)
5. Final verification (optional: `--skip-verify-final`)

**Output Files**:
- `output/genome_index_status_initial.json`
- `output/genome_download_results.json`
- `output/transcriptome_preparation_results.json`
- `output/kallisto_index_build_results.json`
- `output/genome_index_status_final.json`

## Verification

### Verify All Species

Use the verification script to check genome downloads and indexes:

```bash
python3 scripts/rna/verify_genomes_and_indexes.py
```

This will:
- Scan all config files in `config/amalgkit/`
- Check genome download status
- Verify RNA FASTA availability
- Check kallisto index existence
- Generate a summary report with percentages
- Suggest next steps

### Output

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

Results are written to `output/genome_index_status.json` with detailed information for each species.

## Troubleshooting

### Issue: RNA FASTA Not Found

**Symptoms**: `rna_fasta_found: false` in verification output

**Solutions**:
1. Check if genome download completed successfully:
   ```bash
   cat output/amalgkit/{species}/genome/download_record.json
   ```

2. Verify RNA was included in download:
   ```yaml
   genome:
     include:
       - rna  # Must be present
   ```

3. Manually search for RNA FASTA:
   ```bash
   find output/amalgkit/{species}/genome -name "rna.fna*" -o -name "*rna*.fna*"
   ```

### Issue: Kallisto Index Not Built

**Symptoms**: `kallisto_index_found: false` in verification output

**Solutions**:
1. Ensure `build_index: yes` in quant config
2. Check if kallisto is installed:
   ```bash
   which kallisto
   kallisto version
   ```

3. Manually build index:
   ```python
   from metainformant.rna.genome_prep import build_kallisto_index
   from pathlib import Path
   
   fasta_path = Path("output/amalgkit/{species}/work/fasta/{Species}_rna.fasta")
   index_path = Path("output/amalgkit/{species}/work/index/{Species}_transcripts.idx")
   build_kallisto_index(fasta_path, index_path)
   ```

### Issue: Wrong Species Name in Index

**Symptoms**: Index exists but has different naming

**Cause**: Species name in config doesn't match index filename pattern

**Solution**: Index files must match: `{Species_Name}_transcripts.idx` where spaces are replaced with underscores

## Best Practices

### 1. Consistent Naming

- Use underscores in species names: `Camponotus_floridanus`
- Index files follow pattern: `{Species_Name}_transcripts.idx`
- FASTA files follow pattern: `{Species_Name}_rna.fasta`

### 2. Index Building

- Build indexes once per genome version
- Use k-mer size 31 for standard reads (≥100bp)
- Use k-mer size 23 for short reads (<100bp)

### 3. Verification

- Run verification script regularly to track preparation status
- Check verification output before running large quantification jobs
- Document any manual interventions

### 4. Disk Space

- Genome packages can be large (100MB - 1GB+)
- Indexes are typically 20-50MB
- FASTA files are typically 10-100MB
- Consider cleaning up extracted genome packages after index building

## API Reference

### `prepare_genome_for_quantification()`

Complete genome preparation pipeline.

```python
prepare_genome_for_quantification(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
    build_index: bool = True,
    kmer_size: int = 31,
) -> dict[str, Any]
```

**Returns**:
- `success`: bool - Whether preparation succeeded
- `fasta_path`: str | None - Path to prepared FASTA
- `index_path`: str | None - Path to built index
- `error`: str | None - Error message if failed

### `prepare_transcriptome_for_kallisto()`

Extract and prepare transcriptome FASTA.

```python
prepare_transcriptome_for_kallisto(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
) -> Path | None
```

**Returns**: Path to prepared FASTA file or None if failed

### `build_kallisto_index()`

Build kallisto index from FASTA.

```python
build_kallisto_index(
    fasta_path: Path,
    index_path: Path,
    *,
    kmer_size: int = 31,
    check_existing: bool = True,
) -> bool
```

**Returns**: True if index was built successfully or already exists

## Quick Start Guide

### Complete Setup for All Species

```bash
# Step 1: Check current status
python3 scripts/rna/verify_genomes_and_indexes.py

# Step 2: Run complete setup (or use master orchestrator)
python3 scripts/rna/orchestrate_genome_setup.py

# Step 3: Verify final status
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

If you need to download/prepare/build only missing items:

```bash
# Download missing genomes only
python3 scripts/rna/download_missing_genomes.py

# Prepare transcriptomes for newly downloaded genomes
python3 scripts/rna/prepare_transcriptomes.py

# Build indexes for newly prepared transcriptomes
python3 scripts/rna/build_kallisto_indexes.py
```

## See Also

### User Guides
- **[genome_setup_guide.md](genome_setup_guide.md)** - Step-by-step user guide
- **[commands.md](commands.md)** - Command reference for all scripts

### Documentation
- **[API Reference](../API.md#genome-preparation-functions)** - Complete function documentation
- **[Function Index](FUNCTIONS.md)** - Quick function lookup
- **[Quantification Step](steps/quant.md)** - Using kallisto indexes
- **[Pipeline Overview](amalgkit.md)** - Complete workflow documentation
- **[Main Index](../README.md)** - RNA domain master index

