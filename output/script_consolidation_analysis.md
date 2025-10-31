# Scripts Consolidation Analysis

**Date**: October 31, 2025  
**Status**: Analysis Complete + Action Items

---

## Executive Summary

**Current State:**
- 24 scripts across `scripts/`, `scripts/package/`, `scripts/rna/`, and `scripts/rna/amalgkit/`
- Many scripts contain reusable functions that should be in `src/metainformant/rna/`
- `src/metainformant/rna/` already has good structure but scripts have duplicate logic

**Recommendations:**
1. ✅ **Keep current structure** - Scripts are thin wrappers, most logic already in `src/`
2. 🔄 **Minor consolidation** - A few functions can move to `src/`
3. ⚠️ **Action needed** - Quantify 24 unquantified samples, cleanup quantified FASTQs

---

## Analysis: Scripts vs. Source Modules

### ✅ Good Consolidation (Already Done)

**`src/metainformant/rna/steps/`** already contains:
- `sequential_process.py` - Sequential download/quant/delete logic
- `batched_process.py` - Batched processing with skip logic
- `parallel_download.py` - Parallel download with workers
- Individual step modules (`quant.py`, `getfastq.py`, etc.)

These modules **already implement** core logic that scripts use!

### Scripts That Are Properly Thin Wrappers

#### ✅ `run_multi_species_amalgkit.py` and `run_multi_species_sequential.py`
- **Status**: GOOD - Thin wrappers around `src/metainformant/rna/workflow.py`
- **Uses**: `load_workflow_config()`, `execute_workflow()`, `run_amalgkit()`
- **Purpose**: CLI orchestration, not business logic
- **Recommendation**: **Keep as-is**

#### ✅ `monitor_workflow.py`
- **Status**: GOOD - Display logic, not business logic
- **Purpose**: Real-time monitoring dashboard (UI concerns)
- **Recommendation**: **Keep as-is**

#### ✅ Test Scripts (`test_*.py`)
- **Status**: GOOD - Test code doesn't belong in `src/`
- **Recommendation**: **Keep as-is**

### Scripts with Duplicated Logic (Minor Consolidation Opportunities)

#### 🔄 `batch_ena.py` - Some consolidation possible
**Current functions:**
- `get_ena_fastq_urls(sample_id)` - Query ENA API
- `download_fastq_file(url, dest_path, sample_id)` - Download with wget
- `quantify_sample(sample_id, base_dir)` - Run kallisto
- `cleanup_sample(sample_id, base_dir)` - Delete FASTQs
- `process_sample_pipeline(sample_id, base_dir)` - Full pipeline

**Consolidation opportunity:**
Move to `src/metainformant/rna/ena.py`:
```python
# New module: src/metainformant/rna/ena.py
def get_fastq_urls_from_ena(sample_id: str) -> list[str] | None:
    """Query ENA API for direct FASTQ URLs."""
    # Move logic from batch_ena.py
    
def download_ena_fastq(sample_id: str, output_dir: Path) -> bool:
    """Download FASTQ files from ENA."""
    # Move download logic
```

Then `batch_ena.py` becomes:
```python
from metainformant.rna.ena import get_fastq_urls_from_ena, download_ena_fastq
# Much thinner script
```

**Benefit**: Reusable ENA download logic for other scripts/workflows

#### 🔄 `quant_downloaded_samples.py` - Already mostly using src/

**Current functions:**
- `find_unquantified_samples(fastq_dir, quant_dir)` - Find samples needing quant
- `delete_sample_fastqs(sample_id, fastq_dir)` - Delete FASTQs
- `quantify_samples(config_path, species_name)` - Quantify loop

**Status**: **Already uses `src/`** modules:
- Uses `load_workflow_config()` from `src/metainformant/rna/workflow.py`
- Uses `run_amalgkit()` from `src/metainformant/rna/amalgkit.py`
- Uses `read_delimited()`, `write_delimited()` from `src/metainformant/core/io.py`

**Consolidation opportunity:**
Small utility functions could move to `src/metainformant/rna/utils.py`:
```python
# New or existing module: src/metainformant/rna/utils.py
def find_unquantified_samples(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples with FASTQs that haven't been quantified."""
    # Move logic
    
def delete_sample_fastqs(sample_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a sample."""
    # Move logic
```

**Benefit**: Reusable utilities across multiple scripts

### Bash Scripts - Keep as Command-Line Tools

#### ✅ `cleanup_quantified_sra.sh`, `list_unquantified.sh`, `verify_skip_logic.sh`
- **Status**: GOOD - Shell scripts for ad-hoc operations
- **Purpose**: Quick command-line tools, not library code
- **Recommendation**: **Keep as-is**

#### ✅ `force_fasterq.sh`, `force_fasterq_parallel.sh`, `process_one_srr.sh`
- **Status**: GOOD - Shell wrappers for sra-tools commands
- **Purpose**: Bash optimizations, parallel job control
- **Recommendation**: **Keep as-is** (bash is appropriate here)

---

## Recommended Consolidation Plan

### Phase 1: Create New Modules (Low-Hanging Fruit)

#### 1. Create `src/metainformant/rna/ena.py`
```python
"""ENA (European Nucleotide Archive) download utilities."""

from pathlib import Path
import urllib.request
import subprocess

def get_fastq_urls(sample_id: str, timeout: int = 30) -> list[str] | None:
    """Query ENA API for direct FASTQ URLs.
    
    Args:
        sample_id: SRR/ERR/DRR accession
        timeout: HTTP timeout in seconds
        
    Returns:
        List of HTTP URLs for FASTQ files, or None if not found
    """
    # Move logic from batch_ena.py:get_ena_fastq_urls()
    
def download_fastq(url: str, dest_path: Path, timeout: int = 3600) -> tuple[bool, str]:
    """Download a single FASTQ file from ENA.
    
    Args:
        url: Direct HTTP URL to FASTQ file
        dest_path: Destination path
        timeout: Download timeout in seconds
        
    Returns:
        (success: bool, size_mb_or_error: str)
    """
    # Move logic from batch_ena.py:download_fastq_file()
    
def download_sample(sample_id: str, output_dir: Path) -> tuple[bool, dict]:
    """Download all FASTQ files for a sample from ENA.
    
    Args:
        sample_id: Sample accession (SRR/ERR/DRR)
        output_dir: Output directory (will create sample_id/ subdir)
        
    Returns:
        (success: bool, log_entry: dict)
    """
    # Move logic from batch_ena.py:download_sample_ena()
```

#### 2. Create `src/metainformant/rna/utils.py` (or add to existing)
```python
"""RNA-seq workflow utilities."""

from pathlib import Path

def find_unquantified_samples(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples with FASTQs that haven't been quantified.
    
    Args:
        fastq_dir: Directory containing FASTQ files (e.g., fastq/getfastq/)
        quant_dir: Directory containing quantification results
        
    Returns:
        List of sample IDs needing quantification
    """
    # Move from quant_downloaded_samples.py

def delete_sample_fastqs(sample_id: str, fastq_dir: Path) -> bool:
    """Delete FASTQ files for a sample after successful quantification.
    
    Args:
        sample_id: Sample accession
        fastq_dir: Base FASTQ directory
        
    Returns:
        True if deletion successful
    """
    # Move from quant_downloaded_samples.py

def find_quantified_samples_with_fastqs(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples that have been quantified but still have FASTQs.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        quant_dir: Directory containing quantification results
        
    Returns:
        List of sample IDs with quantification and FASTQs (ready for cleanup)
    """
    # New function for cleanup operations
```

### Phase 2: Update Scripts to Use New Modules

#### Update `batch_ena.py`:
```python
from metainformant.rna.ena import get_fastq_urls, download_sample
# Remove duplicated functions, use imported ones
```

#### Update `quant_downloaded_samples.py`:
```python
from metainformant.rna.utils import find_unquantified_samples, delete_sample_fastqs
# Remove duplicated functions, use imported ones
```

### Phase 3: Benefits

✅ **Code Reusability**: ENA download logic available to all scripts and workflows  
✅ **Testing**: Functions in `src/` can be unit tested  
✅ **Maintainability**: Update once, applies everywhere  
✅ **Discoverability**: Users can import from `metainformant.rna.ena`  
✅ **Scripts Stay Simple**: Scripts focus on CLI/orchestration

---

## Current Sample Status

### Unquantified Samples (Need Quantification)

#### P. barbatus
- **Count**: 2 samples
- **Size**: 0.57 GB
- **Samples**: DRR029869, DRR048431
- **Action**: Quantify with kallisto, then delete FASTQs

#### C. floridanus
- **Count**: 22 samples
- **Size**: 56.27 GB
- **Samples**: SRR10960077, SRR2060722, SRR22031369, SRR22031377, SRR22031378, SRR22862547, SRR22862548, SRR22862549, SRR25496790, SRR25496798, SRR25496806, SRR25496814, SRR25496822, SRR25496830, SRR32143976, SRR32143984, SRR32143985, SRR32143986, SRR32143992, SRR5931468, SRR5931475, SRR5931476
- **Action**: Quantify with kallisto, then delete FASTQs

### Quantified Samples (Check for Cleanup)

#### P. barbatus
- **Quantified**: 83 samples
- **FASTQs present**: Yes (found some in `fastq/getfastq/`)
- **Action**: Delete FASTQs for quantified samples to reclaim space

---

## Immediate Action Items

### 1. Quantify Unquantified Samples (24 total)

**P. barbatus (2 samples):**
```bash
python3 scripts/rna/quant_downloaded_samples.py
# Or use the existing quant_downloaded_samples.py script
# which already uses src/ modules appropriately
```

**C. floridanus (22 samples):**
```bash
# Run quantification for cfloridanus
python3 scripts/rna/quant_downloaded_samples.py
```

### 2. Cleanup Quantified Sample FASTQs

**P. barbatus (83 quantified):**
```bash
# Delete FASTQs for all 83 quantified samples
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

**Expected reclaimed space:** ~20-30 GB for pbarbatus

### 3. Implement Consolidation (Optional)

**Priority 1** - High value, low effort:
- Create `src/metainformant/rna/ena.py` with ENA download functions
- Update `batch_ena.py` to use it

**Priority 2** - Medium value, low effort:
- Add utility functions to `src/metainformant/rna/utils.py`
- Update `quant_downloaded_samples.py` to use them

**Priority 3** - Nice to have:
- Add monitoring utilities to `src/metainformant/rna/monitoring.py`
- Move more bash logic to Python where appropriate

---

## Conclusion

### Current Architecture is Good ✅

The repository already has excellent separation:
- `src/metainformant/rna/` contains core business logic
- `scripts/` contains thin CLI wrappers and orchestration
- Most scripts **already use** modules from `src/`

### Minor Improvements Recommended 🔄

1. **Create `src/metainformant/rna/ena.py`** - Reusable ENA download logic
2. **Enhance `src/metainformant/rna/utils.py`** - Common sample utilities
3. **Update 2-3 scripts** to use new modules

### Immediate Priorities ⚠️

1. **Quantify 24 unquantified samples** (2 pbarbatus + 22 cfloridanus)
2. **Delete FASTQs for 83 quantified pbarbatus samples** (reclaim ~20-30 GB)
3. **Optional**: Implement consolidation changes

---

**Summary**: Architecture is already well-designed. Minor consolidation of ENA/utility functions would be beneficial but not critical. Immediate focus should be on processing the 24 unquantified samples and cleaning up quantified FASTQs.

