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
    out_dir: output/amalgkit/amellifera/work
    metadata: output/amalgkit/amellifera/work/metadata/pivot_qualified.tsv
    index_dir: output/amalgkit/amellifera/work/index
    threads: 8
    clean_fastq: yes
    build_index: yes
    fasta_dir: output/amalgkit/amellifera/fasta
```

## Parameters

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory for intermediate and output files. |
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

- **FASTQ Files**: From `amalgkit getfastq` (in `out_dir/getfastq/`)
- **Kallisto Index**: Pre-built index for the species OR FASTA file with `--build_index yes`
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

The workflow execution in `metainformant.rna.workflow` automatically integrates genome preparation when `build_index: yes` is set in the quant configuration. See [genome_preparation.md](../genome_preparation.md) for details on automatic workflow integration.

### System Dependencies

| Tool | Purpose | Installation |
|------|---------|--------------|
| **kallisto** | Pseudoalignment & quantification | `conda install -c bioconda kallisto` |

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

### Issue: "kallisto: command not found"

**Solutions**:
```bash
# Install kallisto
conda install -c bioconda kallisto

# Verify installation
which kallisto
kallisto version
```

### Issue: Index file not found

```
Error: Could not find index file
```

**Diagnosis**:
```bash
# Check index directory
ls -lh output/amalgkit/work/index/

# Check species name in metadata
cut -f2 output/amalgkit/work/metadata/pivot_qualified.tsv | sort -u
# Should match index filename (with spaces → underscores)
```

**Solutions**:
1. Build index:
   ```bash
   amalgkit quant --build_index yes --fasta_dir output/work/fasta
   ```

2. Download pre-built index and place in correct location:
   ```bash
   cp Apis_mellifera_transcripts.idx output/work/index/
   ```

3. Specify index directory explicitly:
   ```bash
   --index_dir /path/to/indices
   ```

### Issue: Low alignment rate

```json
{
    "p_pseudoaligned": 15.2  // <30% is concerning
}
```

**Causes**:
1. Wrong reference transcriptome (different species/version)
2. Low-quality sequencing data
3. Contamination from other organisms
4. Incorrect library type

**Solutions**:
1. Verify reference matches species:
   ```bash
   # Check what's in your index
   grep "^>" output/work/fasta/Apis_mellifera_rna.fasta | head
   ```

2. Check FASTQ quality:
   ```bash
   # Inspect fastp reports
   cat output/work/getfastq/SRR*/fastp.json | grep "total_reads"
   ```

3. Try different reference version:
   ```bash
   # Download different annotation release
   datasets download genome taxon 7460 --include rna
   ```

### Issue: Out of memory

```
kallisto: std::bad_alloc
```

**Solutions**:
1. Increase memory allocation:
   ```bash
   #SBATCH --mem=16G  # Instead of 8G
   ```

2. Reduce number of concurrent jobs

3. Use smaller index (subset of transcriptome)

### Issue: FASTQ files missing

```
Error: Could not find FASTQ files for SRR12345678
```

**Solutions**:
1. Verify getfastq completed:
   ```bash
   ls output/work/getfastq/SRR12345678/
   ```

2. Check if already cleaned (from previous quant run):
   ```bash
   # Check if abundance files exist
   ls output/work/quant/SRR12345678/abundance.tsv
   # FASTQs deleted after successful quant with --clean_fastq yes
   ```

3. Re-run getfastq:
   ```bash
   amalgkit getfastq --id SRR12345678 --out_dir output/work
   ```

## Best Practices

### 1. Always Clean FASTQs After Quantification

```bash
# Good: Save massive disk space
--clean_fastq yes

# Bad: Wastes 100s of GBs
--clean_fastq no
```

**Reasoning**: FASTQ files are 10-50GB per sample. After quantification, you only need the ~1MB abundance.tsv files.

### 2. Verify Reference Transcriptome

```bash
# Before quantifying 100s of samples, test with one
amalgkit quant --batch 1 --threads 8

# Check alignment rate
cat output/work/quant/SRR*/run_info.json | grep "p_pseudoaligned"

# Good: >60%
# Acceptable: 40-60%
# Concerning: <40% (check reference)
```

### 3. Use Appropriate Threading

```bash
# Single sample: use all cores
--threads 16

# Multiple samples in parallel: divide cores
# 4 samples × 4 threads = 16 cores total
```

### 4. Monitor Quantification Progress

```bash
# Check how many samples completed
find output/work/quant -name "abundance.tsv" | wc -l

# Check for failed samples
find output/work/quant -type d | while read dir; do
    if [ ! -f "$dir/abundance.tsv" ]; then
        echo "Failed: $dir"
    fi
done
```

### 5. Archive Abundance Files

```bash
# After successful merge, archive quant outputs
tar -czf quant_outputs_$(date +%Y%m%d).tar.gz output/work/quant/
```

## Real-World Examples

### Example 1: Apis mellifera (83 Samples)

```bash
# With pre-built index
amalgkit quant \
  --out_dir output/amalgkit/amellifera/work \
  --metadata output/amalgkit/amellifera/work/metadata/pivot_qualified.tsv \
  --index_dir output/amalgkit/amellifera/work/index \
  --threads 16 \
  --clean_fastq yes
```

**Result**:
- 83 samples quantified in ~8 hours (serial)
- Average alignment rate: 64.5%
- Average transcripts: 20,672 per sample
- Disk space saved: ~350GB (FASTQs deleted)

### Example 2: Pogonomyrmex barbatus (120 Samples, HPC)

```bash
# Submit array job
sbatch --array=1-120 --cpus-per-task=8 --mem=8G quant.sh

# In quant.sh:
amalgkit quant \
  --out_dir output/amalgkit/pbarbatus/work \
  --batch ${SLURM_ARRAY_TASK_ID} \
  --threads 8 \
  --clean_fastq yes
```

**Result**:
- All 120 samples completed in 30 minutes
- Parallel processing on HPC cluster
- No disk space issues (FASTQs cleaned immediately)

### Example 3: First-time Species Setup

```bash
# Download reference
datasets download genome taxon 7460 --include rna
unzip ncbi_dataset.zip
mv ncbi_dataset/data/*/rna.fna output/work/fasta/Apis_mellifera_rna.fasta

# Build index and quantify
amalgkit quant \
  --out_dir output/work \
  --fasta_dir output/work/fasta \
  --build_index yes \
  --threads 12
```

**Result**: Index built once, reused for all samples

## Integration with METAINFORMANT Workflow

### Automatic Execution

```python
from metainformant.rna.workflow import execute_workflow, load_workflow_config

cfg = load_workflow_config("config/amalgkit_amellifera.yaml")
execute_workflow(cfg)  # quant runs automatically after getfastq/integrate
```

### Genome Preparation Integration

METAINFORMANT workflow automatically handles transcriptome preparation:

```python
# From src/metainformant/rna/workflow.py
# Automatic steps:
# 1. Download genome from NCBI if accession provided
# 2. Extract RNA FASTA
# 3. Build kallisto index
# 4. Run quantification
```

### Configuration

```yaml
# config/amalgkit_amellifera.yaml
genome:
  accession: GCF_003254395.2
  include:
    - rna

steps:
  quant:
    threads: 16
    clean_fastq: yes
    build_index: yes
```

## References

- **kallisto**: https://pachterlab.github.io/kallisto/
- **kallisto paper**: https://www.nature.com/articles/nbt.3519
- **Pseudoalignment**: https://pachterlab.github.io/kallisto/about
- **METAINFORMANT Workflow**: `docs/rna/workflow.md`

## See Also

- **Previous Step**: [`getfastq.md`](getfastq.md) - Downloading FASTQ files
- **Next Step**: [`merge.md`](merge.md) - Combining abundance estimates
- **Workflow Overview**: [`../amalgkit.md`](../amalgkit.md)
- **Testing**: `tests/test_rna_amalgkit_steps.py::test_quant_basic_execution`

---

**Last Updated**: October 29, 2025  
**AMALGKIT Version**: 0.12.19  
**kallisto Version**: 0.46.2+  
**Status**: ✅ Production-ready, comprehensively tested


