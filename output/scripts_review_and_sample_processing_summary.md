# Scripts Review & Sample Processing Summary

**Date**: October 31, 2025  
**Status**: ✅ Analysis Complete | 🔄 Quantification In Progress | ✅ Cleanup Complete

---

## Executive Summary

### 1. Scripts Consolidation Analysis ✅

**Verdict**: **Current architecture is excellent - minimal consolidation needed**

#### What We Found:
- ✅ **24 scripts properly organized** in logical subfolders
- ✅ **Good separation**: `src/metainformant/rna/` contains business logic
- ✅ **Scripts are thin wrappers**: Most already use `src/` modules
- ✅ **Existing modules cover most needs**: `sequential_process.py`, `batched_process.py`, `parallel_download.py`

#### Minor Improvements Recommended:
1. **Create `src/metainformant/rna/ena.py`** - Reusable ENA download functions from `batch_ena.py`
2. **Enhance `src/metainformant/rna/utils.py`** - Sample utility functions from `quant_downloaded_samples.py`

#### Scripts Keeping As-Is:
- All bash scripts (appropriate for shell operations)
- Test scripts (test code doesn't belong in `src/`)
- Orchestration scripts (CLI wrappers, not business logic)
- Monitoring scripts (display/UI logic)

**Recommendation**: Architecture is well-designed. Optional consolidation of ENA/utility functions would be beneficial but not critical.

---

### 2. Sample Processing Actions 🔄

#### Cleanup Completed ✅

**P. barbatus SRA cleanup:**
- ✅ Deleted 18 SRA files for quantified samples
- ✅ Reclaimed **19.83 GB** disk space
- ✅ Skipped 2 unquantified samples (DRR029869, DRR048431)
- 📊 New storage: **35 GB** (down from ~55 GB)

**C. floridanus SRA cleanup:**
- ⏭️  Skipped 16 unquantified samples (awaiting quantification)
- 📊 Current storage: **25 GB**

**Total reclaimed**: **~19.83 GB**

#### Quantification In Progress 🔄

**Running**: `python3 scripts/rna/quant_downloaded_samples.py`

**Samples needing quantification:**

**P. barbatus (2 samples, 0.57 GB):**
1. DRR029869
2. DRR048431

**C. floridanus (16 samples remaining after script check):**
- SRR10960077, SRR2060722, SRR22031369, SRR22031378
- SRR22862549, SRR25496790, SRR25496798, SRR25496806
- SRR25496814, SRR25496822, SRR25496830, SRR32143976
- SRR32143986, SRR32143992, SRR5931468, SRR5931475

**Note**: The quantification script will:
1. Quantify each sample with kallisto
2. Automatically delete FASTQs after successful quantification
3. Log results to `output/quant_unquantified_*.log`

---

## Detailed Analysis

### Scripts Architecture Review

#### ✅ Excellent Design Patterns Found:

**1. Workflow modules properly separated:**
```
src/metainformant/rna/
├── workflow.py              # Workflow orchestration
├── amalgkit.py              # Amalgkit CLI wrapper
├── steps/                   # Individual step modules
│   ├── sequential_process.py
│   ├── batched_process.py
│   └── parallel_download.py
```

**2. Scripts use src/ appropriately:**
- `run_multi_species_*.py` → Uses `workflow.py`, `amalgkit.py`
- `quant_downloaded_samples.py` → Uses `workflow.py`, `io.py`
- Test scripts → Use all relevant modules

**3. Bash scripts handle shell-specific tasks:**
- `force_fasterq.sh` - Parallel job control with bash
- `cleanup_quantified_sra.sh` - File operations with shell tools
- `verify_skip_logic.sh` - Shell-based verification

#### 🔄 Optional Consolidation Opportunities:

**Create `src/metainformant/rna/ena.py`:**
```python
"""ENA (European Nucleotide Archive) download utilities."""

def get_fastq_urls(sample_id: str) -> list[str] | None:
    """Query ENA API for direct FASTQ URLs."""
    # From batch_ena.py:get_ena_fastq_urls()
    
def download_sample(sample_id: str, output_dir: Path) -> tuple[bool, dict]:
    """Download all FASTQ files for a sample from ENA."""
    # From batch_ena.py:download_sample_ena()
```

**Benefits:**
- Reusable ENA download logic across all scripts
- Testable with unit tests
- Discoverable via `from metainformant.rna.ena import ...`

**Enhance `src/metainformant/rna/utils.py`:**
```python
"""RNA-seq workflow utilities."""

def find_unquantified_samples(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find samples with FASTQs that haven't been quantified."""
    # From quant_downloaded_samples.py
    
def find_quantified_with_fastqs(fastq_dir: Path, quant_dir: Path) -> list[str]:
    """Find quantified samples that still have FASTQs (ready for cleanup)."""
    # New utility function
```

**Benefits:**
- Common utilities available to all scripts
- Consistent behavior across workflows
- Easier maintenance

---

## Current Repository Status

### Storage After Cleanup

```
output/amalgkit/
├── pbarbatus/      35 GB  (was ~55 GB, saved 19.83 GB)
├── cfloridanus/    25 GB  (has unquantified samples)
├── mpharaonis/    825 MB
└── sinvicta/      908 MB
```

### Samples Status

| Species | Quantified | Unquantified | With FASTQs | Status |
|---------|------------|--------------|-------------|---------|
| pbarbatus | 83 | 2 | 0 (cleaned) | ✅ Cleanup done, 2 quantifying |
| cfloridanus | ? | 16 | 16 | 🔄 Quantification in progress |
| mpharaonis | ? | 0 | 0 | ✅ No action needed |
| sinvicta | ? | 0 | 0 | ✅ No action needed |

---

## Next Steps

### Immediate (In Progress) 🔄

1. **Monitor quantification**: `tail -f output/quant_unquantified_*.log`
2. **After quantification completes**: FASTQs will be automatically deleted
3. **Verify completion**: Run `bash scripts/rna/list_unquantified.sh` to confirm

### Optional Improvements 📋

1. **Create `src/metainformant/rna/ena.py`**
   - Extract ENA functions from `batch_ena.py`
   - Add unit tests
   - Update `batch_ena.py` to use new module

2. **Enhance `src/metainformant/rna/utils.py`**
   - Add `find_unquantified_samples()`
   - Add `find_quantified_with_fastqs()`
   - Update `quant_downloaded_samples.py` to use utilities

3. **Documentation**
   - Add examples to `docs/rna/ena_downloads.md`
   - Document utility functions in module docstrings

---

## Files Created

1. **`output/script_consolidation_analysis.md`** - Detailed consolidation analysis
2. **`output/quant_unquantified_*.log`** - Quantification process log (in progress)
3. This summary document

---

## Recommendations

### For Immediate Use ✅
- Current architecture is excellent
- Scripts are well-organized and maintainable
- Use existing scripts as-is for all workflows

### For Future Enhancement 🔄
- Consider creating `ena.py` module for reusability
- Consider adding utilities to `utils.py`
- Both are optional improvements, not critical needs

### Maintenance 🔧
- Keep scripts thin (orchestration only)
- Keep business logic in `src/metainformant/`
- Continue pattern of importing from `src/` in scripts

---

**Summary**: Repository has excellent script organization. Successfully reclaimed 19.83 GB by cleaning up quantified samples. Quantification of remaining 18 unquantified samples in progress. Minor consolidation opportunities identified but current architecture is well-designed and maintainable.

