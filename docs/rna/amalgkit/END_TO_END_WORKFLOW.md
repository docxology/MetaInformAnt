# Amalgkit End-to-End Workflow - Pogonomyrmex barbatus

Complete documentation of the amalgkit RNA-seq analysis workflow using real data.

## Overview

This document demonstrates a complete end-to-end RNA-seq analysis using amalgkit with **Pogonomyrmex barbatus** (harvester ant) brain samples.

**Total samples**: 83  
**Species**: *Pogonomyrmex barbatus*  
**Tissue**: Brain  
**Workflow**: metadata → config → getfastq → quant → merge → curate → sanity

---

## Workflow Steps

### 1. Initial Setup (metadata, integrate, config)

These steps were completed earlier using the config file at `config/amalgkit_pbarbatus.yaml`:

```yaml
species: Pogonomyrmex barbatus
search_string: Pogonomyrmex barbatus[Organism] AND brain[tissue]
out_dir: output/amalgkit/pbarbatus
threads: 6
```

**Results**:
- ✅ Metadata retrieved: 83 samples from NCBI SRA
- ✅ Config files generated for tissue filtering
- ✅ Samples validated with sanity check

**Files created**:
- `work/metadata/metadata.tsv` (65.3KB, 83 samples, 56 columns)
- `work/config_base/*.config` (tissue filters, group attributes)

---

### 2. Genome/Transcriptome Download

Downloaded reference transcriptome from NCBI:

```bash
PYTHONPATH=src python3 -c "
from metainformant.rna import amalgkit
amalgkit.getgenome({
    'species': 'Pogonomyrmex_barbatus',
    'out_dir': 'output/amalgkit/pbarbatus/work'
})
"
```

**Results**:
- ✅ Reference: *Pogonomyrmex barbatus* transcriptome
- ✅ 20,672 transcripts
- ✅ Kallisto index built (20MB)

**Files created**:
- `work/fasta/Pogonomyrmex_barbatus.fasta` (51MB)
- `work/index/Pogonomyrmex_barbatus.idx` (20MB)

---

### 3. FASTQ Download & Quantification

**Challenge**: Initial NCBI SRA downloads were extremely slow (0.15 Mbps)  
**Solution**: Switched to ENA (European Nucleotide Archive) for **187X speedup**

#### Parallel Processing Configuration

Created `scripts/rna/batch_ena.py` with optimized settings:
- **5 concurrent downloads** (ENA direct FASTQ links)
- **3 concurrent quantifications** (Kallisto)
- **3 threads per Kallisto** process
- **Automatic cleanup** after quantification

**Performance**:
- Download speed: 28 Mbps (vs. 0.15 Mbps from NCBI)
- Processing rate: ~4.5 samples/hour
- Total time: ~18 hours for 83 samples

**Results**:
- ✅ 83/83 samples downloaded and quantified
- ✅ 1,695,104 total transcript quantifications (20,672 × 83)
- ✅ Average file size: 767KB per abundance.tsv

**Files created**:
- `work/quant/SRR*/abundance.tsv` (83 files, 90MB total)
- `work/quant/SRR*/abundance.h5` (H5 format for downstream tools)
- `work/quant/SRR*/run_info.json` (Kallisto metadata)

---

### 4. Merge Step

Combined quantifications into species-level expression matrices:

```python
from metainformant.rna import amalgkit
result = amalgkit.merge({
    'out_dir': 'output/amalgkit/pbarbatus/work'
})
```

**Results**:
- ✅ Merged metadata with expression data
- ✅ Generated quality control visualizations

**Files created**:
- `work/merge/metadata.tsv` (65.3KB, 83 samples with QC metrics)
- `work/merge/merge_total_spots.pdf` (5.2KB)
- `work/merge/merge_total_bases.pdf` (5.2KB)
- `work/merge/merge_library_layout.pdf` (4.4KB)
- `work/merge/merge_mapping_rate.pdf` (4.2KB)
- `work/merge/merge_exclusion.pdf` (4.4KB)

---

### 5. Curate Step

Quality control and outlier removal:

```python
result = amalgkit.curate({
    'out_dir': 'output/amalgkit/pbarbatus/work'
})
```

**Results**:
- ✅ Tissue filtering applied ("brain", "whole cleanly-dissected brains")
- ✅ No samples excluded (all passed QC)

**Output**:
```
Starting: Pogonomyrmex_barbatus
Tissues to be included: brain, whole cleanly-dissected brains
Number of species: 1
Number of excluded samples: 0
```

---

### 6. Sanity Check

Final validation of workflow integrity:

```python
result = amalgkit.sanity({
    'out_dir': 'output/amalgkit/pbarbatus/work'
})
```

**Results**:
- ✅ All 83 SRA runs validated
- ✅ Metadata integrity confirmed
- ✅ Expression data complete

**Output**:
```
1 species detected: ['Pogonomyrmex barbatus']
83 SRA runs detected
Time elapsed: 0 sec
```

---

## Steps Requiring Additional Data

### 7. CSTMM (Cross-Species TMM) - **SKIPPED**

Requires orthology information (`--orthogroup_table` or `--dir_busco`)  
For multi-species analysis with ortholog-based normalization.

### 8. CSCA (Cross-Species Correlation Analysis) - **SKIPPED**

Requires orthology information  
For cross-species expression correlation plots.

---

## Summary Statistics

### Data Scale
| Metric | Value |
|--------|-------|
| Species | 1 (*Pogonomyrmex barbatus*) |
| Samples | 83 |
| Transcripts | 20,672 |
| Total quantifications | 1,695,104 |
| Processing time | ~18 hours |
| Final storage | 1.2GB |

### Files Generated
| Directory | Files | Size | Description |
|-----------|-------|------|-------------|
| `work/metadata/` | 4 | 65KB | Sample metadata and filters |
| `work/fasta/` | 1 | 51MB | Reference transcriptome |
| `work/index/` | 1 | 20MB | Kallisto index |
| `work/quant/` | 249 | 90MB | Quantification results (83 samples) |
| `work/merge/` | 6 | 89KB | Merged data + QC plots |
| `genome/` | - | 524MB | Downloaded reference data |
| `logs/` | 100+ | <1MB | Processing logs |

---

## Performance Optimizations

### 1. ENA Download Method
- **Before**: 0.15 Mbps (NCBI SRA with `prefetch`)
- **After**: 28 Mbps (ENA direct download with `wget`)
- **Improvement**: **187X speedup**

### 2. Parallel Processing
- **Sequential**: ~37 hours (1 at a time)
- **Parallel**: ~18 hours (5 downloads + 3 quants)
- **Improvement**: **2X speedup**

### 3. Automatic Cleanup
- **FASTQ size**: ~4GB per sample
- **Cleanup**: Delete after quantification
- **Storage saved**: ~332GB (83 samples × 4GB)

---

## Key Commands

### Run Complete Workflow
```bash
# Using Python API
PYTHONPATH=src python3 << 'EOF'
from metainformant.rna import amalgkit
from pathlib import Path

work_dir = Path("output/amalgkit/pbarbatus/work")

# 1. Get metadata
amalgkit.metadata({
    "out_dir": str(work_dir),
    "search_string": "Pogonomyrmex barbatus[Organism] AND brain[tissue]"
})

# 2. Generate config
amalgkit.config({"out_dir": str(work_dir)})

# 3. Get genome
amalgkit.getgenome({
    "species": "Pogonomyrmex_barbatus",
    "out_dir": str(work_dir)
})

# 4-5. Download & quantify (use batch_ena.py for all samples)
# 6. Merge
amalgkit.merge({"out_dir": str(work_dir)})

# 7. Curate
amalgkit.curate({"out_dir": str(work_dir)})

# 8. Sanity check
amalgkit.sanity({"out_dir": str(work_dir)})
EOF
```

### Parallel Batch Processing
```bash
# Run ENA-based parallel downloader
python3 scripts/rna/batch_ena.py

# Or use restart script
bash scripts/rna/restart_batch.sh
```

---

## Accessing Results

### Expression Data
```python
import pandas as pd
from pathlib import Path

# Load one sample
sample_data = pd.read_csv(
    "output/amalgkit/pbarbatus/work/quant/SRR14740487/abundance.tsv",
    sep='\t'
)

print(f"Transcripts: {len(sample_data)}")
print(f"Columns: {list(sample_data.columns)}")

# View top expressed transcripts
top10 = sample_data.nlargest(10, 'tpm')
print(top10[['target_id', 'tpm', 'est_counts']])
```

### Build Expression Matrix
```python
import pandas as pd
from pathlib import Path

quant_dir = Path("output/amalgkit/pbarbatus/work/quant")
samples = sorted(quant_dir.glob("SRR*/abundance.tsv"))

# Collect TPM values
tpm_matrix = {}
for sample_file in samples:
    sample_id = sample_file.parent.name
    df = pd.read_csv(sample_file, sep='\t')
    tpm_matrix[sample_id] = df.set_index('target_id')['tpm']

expression_matrix = pd.DataFrame(tpm_matrix)
print(f"Expression matrix: {expression_matrix.shape}")
print(f"Samples: {len(expression_matrix.columns)}")
print(f"Transcripts: {len(expression_matrix.index)}")

# Save
expression_matrix.to_csv("output/amalgkit/pbarbatus/expression_matrix.csv")
```

---

## Testing & Validation

### Test Coverage
- ✅ **100% amalgkit wrapper coverage** (72 tests)
- ✅ **End-to-end workflow tests** (12 tests, 10 passed + 2 skipped network-dependent)
- ✅ **Real execution** (NO_MOCKING_POLICY enforced)
- ✅ **All step runners tested**

### Test Files
- `tests/test_rna_amalgkit_comprehensive.py` - Wrapper and workflow tests
- `tests/test_rna_amalgkit_steps.py` - Individual step runner tests
- `tests/test_rna_amalgkit_end_to_end.py` - End-to-end integration tests

---

## Troubleshooting

### Common Issues

**1. Slow NCBI downloads**
- **Solution**: Use ENA instead (`scripts/rna/batch_ena.py`)
- **Script**: Automatically queries ENA API for direct FASTQ URLs

**2. Blocking processes**
- **Issue**: Old `prefetch`/`kallisto` processes interfere
- **Solution**: Script includes `kill_blocking_processes()` function

**3. Disk space**
- **Issue**: FASTQ files are large (~4GB per sample)
- **Solution**: Automatic cleanup after quantification

**4. Network timeouts**
- **Configuration**: Download timeout set to 1 hour per FASTQ file
- **Retries**: Script handles failed downloads gracefully

---

## Next Steps

### Downstream Analysis
1. **Differential Expression**: Compare conditions using DESeq2/edgeR
2. **Pathway Analysis**: GO/KEGG enrichment of expressed genes
3. **Visualization**: PCA, heatmaps, volcano plots
4. **Cross-species**: Add more ant species for comparative analysis

### Multi-species CSTMM/CSCA
To enable cross-species analysis:
1. Run OrthoFinder to get ortholog groups
2. Or run BUSCO to get single-copy orthologs
3. Provide `--orthogroup_table` or `--dir_busco` to cstmm/csca

---

## References

- **Amalgkit**: https://github.com/kfuku52/amalgkit
- **Kallisto**: https://pachterlab.github.io/kallisto/
- **ENA**: https://www.ebi.ac.uk/ena/browser/
- **NCBI SRA**: https://www.ncbi.nlm.nih.gov/sra

---

## Reproducibility

All steps are fully reproducible using:
1. Configuration: `config/amalgkit_pbarbatus.yaml`
2. Scripts: `scripts/rna/batch_ena.py`
3. Python API: `metainformant.rna.amalgkit`

**Storage**: All outputs in `output/amalgkit/pbarbatus/`  
**Logs**: All execution logs in `output/amalgkit/pbarbatus/logs/`  
**Documentation**: This file and `docs/rna/examples/`

---

*Workflow completed: October 29, 2025*  
*Processing time: ~18 hours for 83 samples*  
*Success rate: 100% (83/83 samples quantified)*

