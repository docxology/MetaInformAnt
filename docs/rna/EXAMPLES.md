# RNA Workflow Examples

This document collects narrative walkthroughs and quick references for running METAINFORMANT RNA workflows in real-world settings.

## Example Overview

METAINFORMANT RNA workflows support complete end-to-end transcriptomic analysis from raw SRA data to publication-ready results. This document provides real-world examples demonstrating:

- Complete analysis workflows
- Configuration patterns
- Troubleshooting approaches
- Best practices

## Complete Analysis Example: Pogonomyrmex barbatus

**Date**: October 28, 2025  
**Species**: *Pogonomyrmex barbatus* (red harvester ant)  
**Tissue**: Brain  
**Directory**: `output/amalgkit/pbarbatus/`

### Directory Structure

```
output/amalgkit/pbarbatus/
├── work/
│   ├── metadata/
│   │   └── metadata.tsv              # 83 samples discovered from NCBI SRA
│   ├── index/
│   │   └── Pogonomyrmex_barbatus.idx # Kallisto index (20.3MB, ready to use)
│   ├── fasta/
│   │   └── Pogonomyrmex_barbatus.fasta # Transcriptome reference (51.3MB)
│   └── quant/
│       ├── SRR14740487/              # Quantified sample 1
│       │   └── abundance.tsv         # Gene expression values (826KB)
│       └── SRR14740488/              # Quantified sample 2
│           └── abundance.tsv         # Gene expression values (828KB)
├── genome/                           # P. barbatus genome reference
├── fastq/                            # Downloaded FASTQ files (cleaned after quant)
├── logs/                             # Pipeline logs
└── merged/                           # Merged expression matrix (after all samples)
```

### Current Status

| Component | Status | Details |
|-----------|--------|---------|
| **Metadata** | ✅ Complete | 83 samples from NCBI SRA |
| **Genome Reference** | ✅ Downloaded | GCF_000187915.1 |
| **Transcriptome** | ✅ Prepared | 20,672 transcripts |
| **Kallisto Index** | ✅ Built | Ready for quantification |
| **Samples Downloaded** | 2/83 (2.4%) | SRR14740487, SRR14740488 |
| **Samples Quantified** | 2/83 (2.4%) | ~90% expression rate |
| **Expression Matrix** | ⏹️ Pending | Need all 83 samples |

### Completed Samples

#### SRR14740487 (Brain Sample 1)
```
Location:   output/amalgkit/pbarbatus/work/quant/SRR14740487/abundance.tsv
Transcripts: 20,672 total, 18,592 expressed (89.9%)
Reads:       26,078,116 processed, 16,838,576 mapped (64.6%)
Quality:     ✅ Excellent mapping and expression rates

Top 3 genes:
  1. XM_011631231.1 - TPM: 24,372.3
  2. XM_011648336.2 - TPM:  8,560.4
  3. XM_025219454.1 - TPM:  4,694.1
```

#### SRR14740488 (Brain Sample 2)
```
Location:   output/amalgkit/pbarbatus/work/quant/SRR14740488/abundance.tsv
Transcripts: 20,672 total, 18,719 expressed (90.6%)
Reads:       27,085,159 processed, 17,828,830 mapped (65.8%)
Quality:     ✅ Excellent mapping and expression rates

Top 3 genes:
  1. XM_011631231.1 - TPM: 29,342.9
  2. XM_011632595.2 - TPM: 28,000.0
  3. XM_011648336.2 - TPM: 10,236.0
```

### Processing Remaining Samples

**Recommended**: Use `run_workflow.py` for complete end-to-end execution:

```bash
# Full end-to-end workflow (all steps)
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Check status at any time
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status
```

The workflow automatically:
- Downloads remaining samples
- Quantifies each sample immediately
- Deletes FASTQs after quantification
- Merges all samples into expression matrix

### After All Samples Quantified

#### Merge into Expression Matrix

```bash
cd output/amalgkit/pbarbatus
amalgkit merge --out_dir work
```

This will create:
```
output/amalgkit/pbarbatus/work/merge/Pogonomyrmex_barbatus/
├── Pogonomyrmex_barbatus_tc.tsv    # Estimated counts matrix
├── Pogonomyrmex_barbatus_tpm.tsv   # TPM-normalized matrix
└── Pogonomyrmex_barbatus_eff_len.tsv # Effective lengths
```

**Format**: One row per transcript, one column per sample, TPM-normalized values

### Accessing Expression Data

#### Load Individual Samples

```python
import pandas as pd
from pathlib import Path

base = Path("output/amalgkit/pbarbatus/work/quant")

# Load one sample
df = pd.read_csv(base / "SRR14740487" / "abundance.tsv", sep='\t')

# Columns: target_id, length, eff_length, est_counts, tpm
# Use 'tpm' for expression comparisons
```

#### Load Merged Matrix (After All Samples)

```python
import pandas as pd

# Load merged expression matrix
expr = pd.read_csv(
    "output/amalgkit/pbarbatus/work/merge/Pogonomyrmex_barbatus/Pogonomyrmex_barbatus_tpm.tsv",
    sep="\t",
    index_col=0  # transcript_id as index
)

# Shape: (20672 transcripts, 83 samples)
print(expr.shape)

# Filter expressed genes (TPM > 1 in at least 10% of samples)
min_samples = int(0.1 * expr.shape[1])
expressed = expr[(expr > 1).sum(axis=1) >= min_samples]
```

### Storage Management

#### Current Usage
```
Total directory:     599MB
├── Metadata:        ~1MB
├── Genome:          ~520MB
├── Index:           20.3MB
├── Transcriptome:   51.3MB
├── Quantification:  2.2MB (2 samples × ~1MB)
└── Logs:            ~5MB
```

#### After All 83 Samples
```
Expected total:      ~690MB
├── Metadata:        ~1MB
├── Genome:          ~520MB
├── Index:           20.3MB
├── Transcriptome:   51.3MB
├── Quantification:  ~90MB (83 samples × ~1MB)
├── Merged matrix:   ~5MB
└── Logs:            ~10MB

FASTQ files:         0MB (deleted after quantification)
```

**Key Point**: FASTQs are deleted immediately after each successful quantification, keeping storage minimal.

### Quality Metrics

#### Mapping Rates
- **Target**: >60% mapping rate
- **Achieved**: 64.6-65.8% (excellent)

#### Expression Rates
- **Target**: >80% of transcripts detected
- **Achieved**: 89.9-90.6% (excellent)

#### Library Quality
- **Fragment length**: 234.9-236.7bp (consistent)
- **Read depth**: 26-27M reads per sample (good)

### Key Files Reference

#### Metadata
```
output/amalgkit/pbarbatus/work/metadata/metadata.tsv
- 83 samples from NCBI SRA
- Columns: run_accession, scientific_name, tissue, etc.
```

#### Kallisto Index
```
output/amalgkit/pbarbatus/work/index/Pogonomyrmex_barbatus.idx
- Pre-built, ready for quantification
- 20,672 transcripts indexed
- K-mer length: 31
```

#### Expression Values
```
output/amalgkit/pbarbatus/work/quant/{sample_id}/abundance.tsv
- Per-sample gene expression
- Columns: target_id, length, eff_length, est_counts, tpm
- Use TPM for normalized expression
```

#### Merged Matrix (after completion)
```
output/amalgkit/pbarbatus/work/merge/Pogonomyrmex_barbatus/Pogonomyrmex_barbatus_tpm.tsv
- All samples in one matrix
- 20,672 rows (transcripts) × 83 columns (samples)
- TPM-normalized values
```

### Workflow Summary

```
1. ✅ Metadata Discovery       [COMPLETE] 83 samples found
2. ✅ Genome Download          [COMPLETE] Reference genome
3. ✅ Transcriptome Prep       [COMPLETE] 20,672 transcripts
4. ✅ Index Building           [COMPLETE] Kallisto index ready
5. 🔄 Sample Download          [2/83]     Process in batches
6. 🔄 Sample Quantification    [2/83]     Using existing index
7. ⏹️ Merge Expression Matrix  [PENDING]  After all samples
8. ⏹️ Downstream Analysis      [PENDING]  DESeq2, WGCNA, etc.
```

### Consolidated Folder Benefits

1. **Single Source of Truth**: All P. barbatus data in one place
2. **Organized Structure**: Clear hierarchy for all pipeline stages
3. **Reusable Index**: Built once, use for all 83 samples
4. **Efficient Storage**: ~690MB total (vs ~360GB if keeping FASTQs)
5. **Ready to Scale**: Validated workflow for remaining 81 samples

**Status**: ✅ **READY FOR BATCH PROCESSING**  
**Next Step**: Download and quantify remaining 81 samples using existing index  
**Estimated Time**: 10-15 hours for all samples (can run overnight)

## Usage Patterns

### Pattern 1: Single Species End-to-End

```bash
# Complete workflow for one species
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
```

**Use Case**: Focused analysis on a single species with all samples.

### Pattern 2: Multi-Species Comparative

```bash
# Process each species separately
for species in pbarbatus cfloridanus amellifera; do
    python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_${species}.yaml
done

# Then run cross-species analysis
# (See amalgkit steps: cstmm, csca)
```

**Use Case**: Comparative transcriptomics across multiple species.

### Pattern 3: Tissue-Specific Analysis

```yaml
# In config file
filters:
  require_tissue: true
  tissue_types: ["brain", "antenna"]

steps:
  select:
    sample_group: brain,antenna
```

**Use Case**: Focused analysis on specific tissues.

### Pattern 4: Quality-Filtered Analysis

```yaml
# In config file
filters:
  min_spots: 10000000  # Minimum 10M reads
  library_layout: PAIRED
  platform: Illumina

steps:
  select:
    min_nspots: 10000000
    max_sample: 100
```

**Use Case**: High-quality samples only for publication-ready analyses.

## Maintenance

- Update the guides whenever workflows add new steps or configuration parameters
- Cross-link to `docs/rna/amalgkit/steps/README.md` and `docs/rna/workflow.md` when functionality changes
- Ensure code snippets stay executable with `python3 scripts/rna/run_workflow.py`

## Related Documentation

- **[workflow.md](workflow.md)**: Workflow planning and execution
- **[amalgkit/steps/README.md](amalgkit/steps/README.md)**: Individual step documentation
- **[CONFIGURATION.md](CONFIGURATION.md)**: Configuration management
- **[GETTING_STARTED.md](GETTING_STARTED.md)**: Setup and installation

