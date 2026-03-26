# amalgkit quant: Transcript Abundance Estimation

## Purpose

Quantifies transcript abundance from FASTQ files using kallisto pseudoalignment. This step performs the **core quantification** of gene expression levels from raw sequencing reads.

## Overview

The `quant` step:
- Runs kallisto quant on each FASTQ file
- Performs pseudoalignment to reference transcriptome
- Generates abundance estimates (TPM, counts)
- Optionally builds kallisto indices from FASTA files
- Supports batch processing for HPC environments
- Automatically cleans up FASTQs after successful quantification (configurable)

## Function Signature

### Python API

```python
from metainformant.rna import amalgkit

# High-level wrapper
def quant(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]

# Step runner (low-level)
from metainformant.rna.engine.workflow_steps import run_quant

def run_quant(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**See Also**: [API Reference](../../API.md#quant) | [Function Index](../FUNCTIONS.md)

## Usage

### Basic Usage

```bash
amalgkit quant \
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/pivot_qualified.tsv \
  --index_dir output/amalgkit/work/index \
  --threads 8 \
  --clean_fastq yes
```

### Python API

```python
from metainformant.rna import amalgkit

result = amalgkit.quant(
    out_dir="output/amalgkit/work",
    metadata="output/amalgkit/work/metadata/pivot_qualified.tsv",
    index_dir="output/amalgkit/work/index",
    threads=8,
    clean_fastq=True
)
```

### Configuration File

```yaml
steps:
  quant:
    # CRITICAL: out_dir must be work_dir (not separate quant_dir)
    # This allows quant to find getfastq output in {out_dir}/getfastq/
    out_dir: output/amalgkit/amellifera/work
    metadata: output/amalgkit/amellifera/work/metadata/pivot_qualified.tsv
    index_dir: output/amalgkit/amellifera/work/index
    threads: 16
    clean_fastq: yes
    build_index: yes
    fasta_dir: output/amalgkit/amellifera/work/fasta
```

## Parameters

### CRITICAL: out_dir Configuration

**IMPORTANT**: For `amalgkit quant`, the `out_dir` parameter should be set to the `work_dir` (not a separate quant directory). This is because:

1. **FASTQ Location**: `amalgkit quant` looks for FASTQ files in `{out_dir}/getfastq/{sample_id}/`. If `getfastq` used `out_dir: output/amalgkit/{species}/fastq`, then quant's `out_dir` must be set to `work_dir` so it can find the getfastq output via a symlink or the same directory structure.

2. **Output Location**: Quantification results are written to `{out_dir}/quant/{sample_id}/abundance.tsv`. If `out_dir` is `work_dir`, results will be in `work_dir/quant/`.

3. **Path Resolution**: Setting `out_dir` to `work_dir` ensures that all workflow steps can find files in expected locations relative to the working directory.

**Correct Configuration**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work  # ← Use work_dir, not separate quant_dir
    threads: 16
```

**Incorrect Configuration** (will fail to find FASTQ files):
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/quant  # ← WRONG: quant can't find getfastq output
    threads: 16
```

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory for intermediate and output files. **Should be set to `work_dir`** (not a separate quant directory) so quant can find getfastq output in `{out_dir}/getfastq/`. |
| `--metadata` | PATH | `inferred` | Path to metadata.tsv. Default: `out_dir/metadata/metadata.tsv` |
| `--threads` | INT | `1` | Number of threads for kallisto. |
| `--redo` | yes/no | `no` | Re-quantify even if output exists. |
| `--batch` | INT | `None` | Process only one SRA record (1-based index). For HPC array jobs. |
| `--index_dir` | PATH | `None` | Path to directory containing kallisto index files. Required if not `out_dir/index/` |
| `--clean_fastq` | yes/no | `yes` | **Remove FASTQ files** after successful quantification. Saves massive disk space! |
| `--fasta_dir` | PATH | `inferred` | Path to reference transcriptome FASTA files. Default: `out_dir/fasta` |
| `--build_index` | yes/no | `no` | Build kallisto index from FASTA if index doesn't exist. |

## Input Requirements

### Prerequisites

- **FASTQ Files**: From `amalgkit getfastq` (in `{out_dir}/getfastq/{sample_id}/`)
  - **Path Resolution**: `amalgkit quant` looks for FASTQ files in `{out_dir}/getfastq/{sample_id}/` relative to its `out_dir` parameter
  - **Critical**: If `getfastq` used `out_dir: output/amalgkit/{species}/fastq`, then quant's `out_dir` must be set to `work_dir` (e.g., `output/amalgkit/{species}/work`) so it can find the getfastq output
  - The workflow may create symlinks from `{work_dir}/getfastq/` to `{fastq_dir}/getfastq/` to ensure quant can find the files
- **Kallisto Index**: Pre-built index for the species OR FASTA file with `--build_index yes`
  - Index location: `{out_dir}/index/{Scientific_Name}_transcripts.idx`
  - FASTA location: `{out_dir}/fasta/{Scientific_Name}_rna.fasta`
- **kallisto**: Installed and available on PATH

### Kallisto Index Requirements

**Option 1: Pre-built Index**
```bash
# Index files must be named:
# {Scientific_Name}_transcripts.idx
# Example: Apis_mellifera_transcripts.idx

mkdir -p output/amalgkit/work/index
# Place index file in this directory
```

**Option 2: Build from FASTA**
```bash
# FASTA files must be named:
# {Scientific_Name}*.fasta or {Scientific_Name}*.fa
# Example: Apis_mellifera_transcripts.fasta

mkdir -p output/amalgkit/work/fasta
# Place FASTA file in this directory

# Enable index building
--build_index yes --fasta_dir output/amalgkit/work/fasta
```

**Option 3: Automated Genome Setup (Recommended)**

METAINFORMANT provides automated scripts to download genomes, prepare transcriptomes, and build kallisto indexes:

```bash
# Complete setup for all species
python3 scripts/rna/orchestrate_genome_setup.py

# Or step-by-step:
python3 scripts/rna/verify_genomes_and_indexes.py        # Check status
python3 scripts/rna/download_missing_genomes.py          # Download genomes
python3 scripts/rna/prepare_transcriptomes.py            # Prepare FASTA
python3 scripts/rna/build_kallisto_indexes.py           # Build indexes
```

These scripts automatically:
- Download genome packages from NCBI (best-effort: CLI → API → FTP)
- Extract RNA FASTA files from downloaded genomes
- Prepare transcriptome FASTA files in the correct location
- Build kallisto indexes with proper naming

**See Also:**
- **[genome_preparation.md](../genome_preparation.md)**: Technical documentation for automatic index building
- **[genome_setup_guide.md](../genome_setup_guide.md)**: Complete step-by-step guide for genome setup
- **[commands.md](../commands.md)**: Command reference for all genome setup scripts

The workflow execution in `metainformant.rna.engine.workflow` automatically integrates genome preparation when `build_index: yes` is set in the quant configuration. See [genome_preparation.md](../genome_preparation.md) for details on automatic workflow integration.

### System Dependencies

| Tool | Purpose | Installation |
|------|---------|--------------|
| **kallisto** | Pseudoalignment & quantification | `conda install -c bioconda kallisto` or `sudo apt-get install -y kallisto` |

### Reference Transcriptome Sources

- **NCBI RefSeq**: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
- **Ensembl**: https://ftp.ensembl.org/pub/
- **NCBI Datasets**: `datasets download genome taxon <TAXON_ID> --include rna`

## Output Files

### Directory Structure

```
out_dir/quant/
└── SRR12345678/
    ├── abundance.tsv               # ⭐ Transcript-level abundance estimates
    ├── abundance.h5                 # HDF5 format abundance
    ├── run_info.json                # Kallisto run information
    └── kallisto_output.log          # Kallisto stdout/stderr
```

### Primary Output: `abundance.tsv`

Tab-delimited file with columns:

| Column | Description |
|--------|-------------|
| `target_id` | Transcript ID (from reference FASTA) |
| `length` | Transcript length (bp) |
| `eff_length` | Effective length (accounts for fragment distribution) |
| `est_counts` | Estimated read counts |
| `tpm` | Transcripts Per Million (normalized expression) |

**Example**:
```
target_id               length  eff_length  est_counts  tpm
XM_006566778.3          1821    1672.00     1234.00     15.23
XM_006566779.2          2145    1996.00     567.50      5.87
XM_006566780.1          987     838.00      89.00       2.19
```

### Run Information: `run_info.json`

JSON file with quantification metadata:
```json
{
    "n_targets": 20672,
    "n_bootstraps": 0,
    "n_processed": 5876543,
    "n_pseudoaligned": 3789012,
    "n_unique": 3245678,
    "p_pseudoaligned": 64.5,
    "p_unique": 55.2,
    "kallisto_version": "0.46.2",
    "index_version": 10,
    "start_time": "Thu Oct 24 21:17:56 2024",
    "call": "kallisto quant ..."
}
```

**Key Metrics**:
- `n_processed`: Total reads processed
- `n_pseudoaligned`: Reads successfully aligned
- `p_pseudoaligned`: **Alignment rate** (%)
- `p_unique`: Unique alignment rate (%)

## Workflow Integration

### Position in Pipeline

```mermaid
flowchart LR
    A[getfastq] --> B[integrate]
    B --> C[quant]
    C --> D[merge]
```

**quant** runs **after getfastq/integrate**, **before merge**.

### Downstream Dependencies

| Step | Dependency | Description |
|------|------------|-------------|
| `merge` | `abundance.tsv` files | Combines all samples into expression matrix |
| `curate` | Merged matrix | QC and batch correction |

## Performance Considerations

### Runtime

**Per Sample** (typical 5M-20M reads):
- **Single-end, 8 threads**: 2-5 minutes
- **Paired-end, 8 threads**: 5-10 minutes
- **Large sample (50M+ reads)**: 15-30 minutes

**For 100 Samples**:
- **Serial**: 8-16 hours
- **Parallel** (10 concurrent): 1-2 hours
- **HPC array**: 10-30 minutes (all samples simultaneously)

### Resource Usage

**Per Sample**:
- **CPU**: Scales well up to 8-16 threads
- **Memory**: 4-8GB (depends on index size)
- **Disk I/O**: Heavy read from FASTQs

**Index Loading**:
- Index is loaded into memory once per run
- ~2-4GB RAM for typical transcriptome

### Parallelization Strategies

**Strategy 1: Multi-threading per sample**
```bash
# Good for single-sample processing
--threads 16
```

**Strategy 2: Multiple samples in parallel**
```bash
# Better for many samples: run 4 samples × 4 threads each
for i in {1..4}; do
    amalgkit quant --batch $i --threads 4 &
done
wait
```

**Strategy 3: HPC Array Jobs**
```bash
# Best: Process all samples simultaneously on cluster
sbatch --array=1-100 quant_job.sh
# Each job: --batch ${SLURM_ARRAY_TASK_ID} --threads 8
```

## Common Use Cases

### 1. Quantify All Samples

```bash
amalgkit quant \
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/pivot_qualified.tsv \
  --index_dir output/amalgkit/work/index \
  --threads 8 \
  --clean_fastq yes
```

**Result**: Quantifies all samples, removes FASTQs to save space

### 2. Build Index and Quantify

```bash
# First time with new species
amalgkit quant \
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/pivot_qualified.tsv \
  --fasta_dir output/amalgkit/work/fasta \
  --build_index yes \
  --threads 12
```

**Result**: Builds kallisto index, then quantifies all samples

### 3. Re-quantify Failed Samples

```bash
# Get list of samples without abundance files
# Then re-run
amalgkit quant \
  --out_dir output/amalgkit/work \
  --metadata failed_samples.tsv \
  --redo yes \
  --threads 8
```

### 4. HPC Array Job

```bash
# SLURM job script
#SBATCH --array=1-100
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00

amalgkit quant \
  --out_dir output/amalgkit/work \
  --batch ${SLURM_ARRAY_TASK_ID} \
  --threads 8 \
  --clean_fastq yes
```

**Result**: 100 samples quantified in parallel, ~1 hour total

### 5. Keep FASTQs for Later Analysis

```bash
# Don't delete FASTQs
amalgkit quant \
  --out_dir output/amalgkit/work \
  --threads 8 \
  --clean_fastq no
```

**Warning**: Requires substantial disk space!

## Kallisto Index Building

### Automatic Index Building

```bash
# amalgkit will build index if:
# 1. --build_index yes
# 2. Index doesn't exist at index_dir/{Species_Name}_transcripts.idx
# 3. FASTA exists at fasta_dir/{Species_Name}*.fasta

amalgkit quant \
  --out_dir output/work \
  --fasta_dir output/work/fasta \
  --build_index yes
```

### Manual Index Building

```bash
# Build index manually for more control
kallisto index \
  -i output/work/index/Apis_mellifera_transcripts.idx \
  -k 31 \
  output/work/fasta/Apis_mellifera_rna.fasta

# Then quantify
amalgkit quant \
  --out_dir output/work \
  --index_dir output/work/index
```

**Recommended k-mer size**:
- `-k 31`: Standard (reads ≥100bp)
- `-k 23`: Short reads (<100bp)

### Index Naming Convention

**CRITICAL**: Index files must match species name in metadata:

| Species in metadata | Expected index filename |
|---------------------|-------------------------|
| `Apis mellifera` | `Apis_mellifera_transcripts.idx` |
| `Pogonomyrmex barbatus` | `Pogonomyrmex_barbatus_transcripts.idx` |
| `Drosophila melanogaster` | `Drosophila_melanogaster_transcripts.idx` |

**Pattern**: Replace spaces with underscores, append `_transcripts.idx`

## Troubleshooting

See **[06_quant_troubleshooting.md](06_quant_troubleshooting.md)** for common issues:
- "kallisto: command not found"
- Index file not found
- Low alignment rate
- Out of memory
- FASTQ files missing

## Advanced Usage & Best Practices

See **[06_quant_advanced.md](06_quant_advanced.md)** for:
- Automatic FASTQ cleanup after quantification
- Manual and batch cleanup workflows
- Best practices (threading, verification, monitoring)
- Real-world examples (single species, HPC, workflow integration)

## See Also

### Related Steps
- **[04_getfastq.md](04_getfastq.md)** - Previous: Download FASTQ files
- **[07_merge.md](07_merge.md)** - Next: Merge quantification results

### Documentation
- **[API Reference](../../API.md#quant)** - Complete function documentation
- **[Function Index](../FUNCTIONS.md)** - Quick function lookup
- **[Pipeline Overview](../amalgkit.md)** - Complete pipeline documentation
- **[Steps Index](README.md)** - All step documentation

---

**Last Updated**: March 8, 2026
**AMALGKIT Version**: 0.16.0
**kallisto Version**: ≥0.50.0
**Status**: ✅ Production-ready


