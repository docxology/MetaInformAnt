# Amalgkit Path Resolution Guide

**Complete reference for file path resolution in amalgkit RNA-seq workflows**

## Overview

This guide explains how amalgkit resolves file paths for each workflow step, based on the actual working implementation. Understanding these path relationships is critical for configuring workflows correctly.

## Key Principles

1. **getfastq creates subdirectories**: The `getfastq` step automatically creates a `getfastq/` subdirectory within its specified `out_dir`
2. **quant uses work_dir**: The `quant` step should use `out_dir = work_dir` (not a separate quant_dir) so it can find getfastq output
3. **Path auto-adjustment**: The workflow automatically adjusts paths when possible, but explicit configuration is recommended
4. **Relative paths**: All paths in config files are resolved relative to the repository root

## Directory Structure

```
output/amalgkit/{species}/
├── work/                          ← Working directory (used as quant.out_dir)
│   ├── fasta/                     ← Prepared transcriptome
│   │   └── {Species_Name}_rna.fasta
│   ├── index/                      ← Kallisto index
│   │   └── {Species_Name}_transcripts.idx
│   ├── metadata/                  ← Sample metadata
│   │   ├── metadata.tsv
│   │   ├── metadata_selected.tsv
│   │   └── metadata_integrated.tsv
│   ├── quant/                      ← Quantification results (when quant.out_dir = work_dir)
│   │   └── {sample_id}/
│   │       └── abundance.tsv
│   ├── getfastq/                   ← Symlink to fastq/getfastq (if needed for quant)
│   │   └── {sample_id}/
│   │       └── {sample_id}_1.fastq.gz
│   └── merge/                      ← Merged expression matrices (if merge.out_dir = work_dir)
│       └── {Scientific_Name}/
│           └── {Scientific_Name}_tpm.tsv
├── fastq/                          ← FASTQ files (getfastq out_dir)
│   └── getfastq/                   ← Automatically created by amalgkit getfastq
│       └── {sample_id}/
│           ├── {sample_id}_1.fastq.gz
│           └── {sample_id}_2.fastq.gz
└── merged/                         ← Merged results (if merge uses separate out_dir)
    └── merge/
        └── {Scientific_Name}/
            └── {Scientific_Name}_tpm.tsv
```

## Step-by-Step Path Resolution

### 1. getfastq Step

**Configuration**:
```yaml
steps:
  getfastq:
    out_dir: output/amalgkit/{species}/fastq
```

**Actual Output Location**:
- FASTQ files are created in: `{out_dir}/getfastq/{sample_id}/`
- Example: `output/amalgkit/{species}/fastq/getfastq/SRR12345678/SRR12345678_1.fastq.gz`

**Key Points**:
- The `getfastq/` subdirectory is **automatically created** by amalgkit
- FASTQ files are NOT placed directly in `{out_dir}/`, but in `{out_dir}/getfastq/`
- This subdirectory structure is required by amalgkit's internal logic

**Why This Matters**:
- The `integrate` step must point to `{out_dir}/getfastq/` (the subdirectory), not just `{out_dir}/`
- The `quant` step looks for FASTQ files in `{quant_out_dir}/getfastq/{sample_id}/`

### 2. integrate Step

**Configuration**:
```yaml
steps:
  integrate:
    fastq_dir: output/amalgkit/{species}/fastq/getfastq
```

**Path Resolution**:
- **Required**: `fastq_dir` must point to the `getfastq/` subdirectory created by the `getfastq` step
- **Auto-adjustment**: The workflow automatically adjusts `fastq_dir` if the `getfastq` subdirectory exists
- **Inference**: If `fastq_dir` is not specified, the workflow infers it from the `getfastq` step's `out_dir`

**Example**:
- If `getfastq.out_dir: output/amalgkit/{species}/fastq`
- Then `integrate.fastq_dir` should be: `output/amalgkit/{species}/fastq/getfastq`

**What integrate does**:
- Scans `fastq_dir` for FASTQ files
- Adds file paths and statistics to metadata
- Creates `metadata_integrated.tsv` with local FASTQ information

### 3. quant Step

**Configuration**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work  # ← CRITICAL: Use work_dir, not separate quant_dir
```

**Path Resolution**:
- **CRITICAL**: `out_dir` should be set to `work_dir` (not a separate quant_dir)
- **FASTQ Location**: `amalgkit quant` looks for FASTQ files in `{out_dir}/getfastq/{sample_id}/`
- **Output Location**: Quantification results are written to `{out_dir}/quant/{sample_id}/abundance.tsv`
- **Index Location**: Kallisto index should be in `{out_dir}/index/{Scientific_Name}_transcripts.idx`
- **FASTA Location**: Transcriptome FASTA should be in `{out_dir}/fasta/{Scientific_Name}_rna.fasta`

**Why work_dir is Required**:
1. `amalgkit quant` expects to find FASTQ files in `{out_dir}/getfastq/{sample_id}/`
2. If `getfastq` used `out_dir: output/amalgkit/{species}/fastq`, then quant's `out_dir` must be `work_dir` so it can find the getfastq output
3. The workflow may create a symlink from `{work_dir}/getfastq/` to `{fastq_dir}/getfastq/` to ensure quant can find the files

**Incorrect Configuration** (will fail):
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/quant  # ← WRONG: quant can't find getfastq output
```

**Correct Configuration**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work  # ← CORRECT: quant can find getfastq output
```

**Output Files**:
- `{out_dir}/quant/{sample_id}/abundance.tsv` - Transcript abundance estimates
- `{out_dir}/quant/{sample_id}/abundance.h5` - Binary abundance file
- `{out_dir}/quant/{sample_id}/run_info.json` - Run metadata

### 4. merge Step

**Configuration**:
```yaml
steps:
  merge:
    out_dir: output/amalgkit/{species}/work  # ← Use same work_dir as quant, OR
    # out_dir: output/amalgkit/{species}/merged  # ← Use separate merged directory (requires copying files)
```

**Path Resolution**:
- **Input Location**: `amalgkit merge` looks for abundance files in `{out_dir}/quant/{sample_id}/{sample_id}_abundance.tsv`
- **Output Location**: Merged matrices are written to `{out_dir}/merge/{Scientific_Name}/{Scientific_Name}_tpm.tsv`

**Two Configuration Options**:

**Option 1: Use same work_dir as quant** (Recommended):
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work
  merge:
    out_dir: output/amalgkit/{species}/work  # Same as quant
```
- Merge finds quant output directly in `{work_dir}/quant/`
- No file copying needed

**Option 2: Use separate merged directory**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work
  merge:
    out_dir: output/amalgkit/{species}/merged  # Different from quant
```
- Requires copying abundance files to `{merged_dir}/quant/{sample_id}/{sample_id}_abundance.tsv`
- The workflow may handle this automatically, but explicit configuration is clearer

**Output Files**:
- `{out_dir}/merge/{Scientific_Name}/{Scientific_Name}_tc.tsv` - Estimated counts matrix
- `{out_dir}/merge/{Scientific_Name}/{Scientific_Name}_tpm.tsv` - TPM (normalized) matrix
- `{out_dir}/merge/{Scientific_Name}/{Scientific_Name}_eff_len.tsv` - Effective lengths matrix
- `{out_dir}/merge/{Scientific_Name}/{Scientific_Name}_sample_info.tsv` - Sample metadata

**R Dependency**:
- The merge step requires R and `ggplot2` for generating visualization plots
- If R/ggplot2 is not available, merge will still create expression matrices but will skip plot generation
- See [R_INSTALLATION.md](R_INSTALLATION.md) for installation instructions

## Common Path Issues and Solutions

### Issue 1: integrate can't find FASTQ files

**Error**:
```
ValueError: No detected fastq files (...) in: output/amalgkit/{species}/fastq
```

**Cause**: `fastq_dir` points to the base fastq directory, not the `getfastq` subdirectory

**Solution**:
```yaml
steps:
  integrate:
    fastq_dir: output/amalgkit/{species}/fastq/getfastq  # ← Add /getfastq
```

### Issue 2: quant can't find FASTQ files

**Error**:
```
getfastq output not found in: {out_dir}/getfastq/{sample_id}
```

**Cause**: `quant.out_dir` is set to a separate quant_dir instead of work_dir

**Solution**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work  # ← Use work_dir, not quant_dir
```

### Issue 3: merge can't find abundance files

**Error**:
```
quant outfile not found: {out_dir}/quant/{sample_id}/{sample_id}_abundance.tsv
```

**Cause**: `merge.out_dir` doesn't match `quant.out_dir`, or files aren't in expected location

**Solution**:
```yaml
steps:
  quant:
    out_dir: output/amalgkit/{species}/work
  merge:
    out_dir: output/amalgkit/{species}/work  # ← Use same out_dir as quant
```

## Configuration Best Practices

### Recommended Configuration Pattern

```yaml
work_dir: output/amalgkit/{species}/work

steps:
  getfastq:
    out_dir: output/amalgkit/{species}/fastq
    # Creates: output/amalgkit/{species}/fastq/getfastq/{sample_id}/
  
  integrate:
    fastq_dir: output/amalgkit/{species}/fastq/getfastq
    # Points to getfastq subdirectory
  
  quant:
    out_dir: output/amalgkit/{species}/work  # ← Use work_dir
    # Looks for FASTQ in: {work_dir}/getfastq/{sample_id}/
    # Writes output to: {work_dir}/quant/{sample_id}/abundance.tsv
  
  merge:
    out_dir: output/amalgkit/{species}/work  # ← Use same as quant
    # Looks for abundance in: {work_dir}/quant/{sample_id}/{sample_id}_abundance.tsv
    # Writes output to: {work_dir}/merge/{Scientific_Name}/{Scientific_Name}_tpm.tsv
```

## Path Resolution in Code

The workflow automatically handles some path adjustments:

1. **integrate fastq_dir adjustment**: `src/metainformant/rna/workflow.py::plan_workflow()` automatically adjusts `fastq_dir` to include the `getfastq` subdirectory if it exists

2. **vdb-config repository path**: `src/metainformant/rna/workflow.py::execute_workflow()` sets `vdb-config` repository root to the getfastq directory so `prefetch` downloads to the correct location

3. **Symlink creation**: The workflow may create symlinks from `{work_dir}/getfastq/` to `{fastq_dir}/getfastq/` to ensure quant can find FASTQ files

## See Also

- **[getfastq.md](steps/04_getfastq.md)**: Complete getfastq step documentation
- **[integrate.md](steps/05_integrate.md)**: Complete integrate step documentation
- **[quant.md](steps/06_quant.md)**: Complete quant step documentation
- **[merge.md](steps/07_merge.md)**: Complete merge step documentation
- **[FILE_PATH_STORAGE.md](../FILE_PATH_STORAGE.md)**: Complete file path storage reference
- **[R_INSTALLATION.md](R_INSTALLATION.md)**: R and ggplot2 installation guide

---

**Last Updated**: January 2026  
**Status**: ✅ Production-ready, comprehensively documented


