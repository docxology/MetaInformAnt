# FASTQ Extraction Improvements - Complete Summary

## Date: 2026-01-04

## Problem Summary

The RNA workflow was failing to extract FASTQ files from SRA files, resulting in:
- ✅ SRA files downloaded successfully
- ❌ 0 FASTQ files extracted (all samples showed "Sum of fastp input reads: 0 bp")
- ❌ Workflow continued to downstream steps despite extraction failure
- ❌ Cascading failures in integrate, quant, merge steps

## Root Causes Identified

### 1. Disk Space Issues in `/tmp`
- `/tmp` is a 16GB tmpfs (RAM-based filesystem) that was **100% full**
- `fasterq-dump` uses `/tmp` for temporary files during extraction
- When `/tmp` is full, fasterq-dump fails with "disk-limit exceeded!" error
- This caused extraction to fail silently, leaving SRA files but no FASTQ files

### 2. LITE SRA Files
- Some samples have LITE SRA files (metadata-only, no sequence data)
- LITE files contain "lite" in their URLs (e.g., `SRR123456.lite.1`)
- Even AWS has some LITE files (2 out of 15 samples in test dataset)
- LITE files cannot be converted to FASTQ (no sequence data)

### 3. Missing Early Exit Logic
- Workflow continued to integrate/quant/merge steps even when getfastq failed
- No validation check to detect missing FASTQ files
- Cascading failures in all downstream steps

## Solutions Implemented

### 1. Automatic TMPDIR and VDB_CONFIG Configuration ✅

**Files Modified**:
- `src/metainformant/rna/amalgkit.py` (lines 299-310, 383-393)
- `src/metainformant/rna/workflow.py` (lines 422-444)

**Implementation**:
- Automatically sets `TMPDIR`, `TEMP`, and `TMP` environment variables to repository's `.tmp/fasterq-dump` directory before running getfastq step
- Sets `VDB_CONFIG` environment variable to repository's `.tmp/vdb` directory
- Configures vdb-config repository path at workflow initialization (non-interactive)
- Works for both monitored and non-monitored subprocess calls
- Creates directories if they don't exist

**Benefits**:
- Prevents "disk-limit exceeded" errors when `/tmp` is full
- Uses external drive space (5.5TB) instead of limited tmpfs
- Works automatically without user configuration
- No manual intervention required

### 2. LITE File Filtering ✅

**Files Modified**:
- `src/metainformant/rna/metadata_filter.py` (lines 18-24, 86-95)
- `src/metainformant/rna/workflow.py` (line 664-668)

**Implementation**:
- Added `exclude_lite_files` parameter to `filter_selected_metadata()` function (default: `True`)
- Automatically detects LITE files by checking `AWS_Link`, `NCBI_Link`, and `GCP_Link` columns for "lite" in URLs
- Filters out LITE files before getfastq step
- Logs warning when LITE files are excluded

**Benefits**:
- Prevents downloading metadata-only SRA files
- Reduces wasted time and disk space
- Clear logging of excluded samples
- Works automatically after select step

### 3. Early Exit Logic ✅

**Files Modified**:
- `src/metainformant/rna/workflow.py` (lines 534-550)

**Implementation**:
- Validation runs automatically after getfastq step completes
- Detects when 0 samples have FASTQ files extracted
- Stops workflow immediately with clear error message
- Provides remediation steps

**Benefits**:
- Prevents cascading failures in downstream steps
- Saves time by stopping early
- Clear error messages with actionable remediation steps

### 4. Enhanced Pre-Step Validation ✅

**Files Modified**:
- `src/metainformant/rna/workflow.py` (lines 488-639)

**Implementation**:
- Pre-step validation checks before integrate, quant, and merge steps
- Checks for FASTQ files before integrate step
- Checks for FASTQ files and quantification tools before quant step
- Checks for quantification output and R dependencies before merge step
- Provides clear error messages with remediation steps

**Benefits**:
- Catches missing prerequisites early
- Prevents confusing errors in downstream steps
- Clear guidance on how to fix issues

## Testing Status

### ✅ Completed
- Early exit logic triggers correctly when validation detects 0 FASTQ files
- TMPDIR and VDB_CONFIG environment variables are set correctly
- vdb-config repository path can be configured non-interactively
- LITE file filtering works correctly (tested on metadata with 2 LITE files)
- Pre-step validation checks work correctly

### ⏳ Pending
- **Full workflow test**: Re-run getfastq step to verify FASTQ extraction works with TMPDIR/VDB_CONFIG fixes
- **LITE file filtering test**: Verify that LITE files are excluded from getfastq step
- **End-to-end test**: Run complete workflow from metadata → getfastq → quant to verify all fixes work together

## Next Steps

1. **Re-run getfastq step** to verify FASTQ extraction works:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus_5sample.yaml --steps getfastq
   ```

2. **Verify LITE files are excluded**:
   - Check metadata_selected.tsv after select step
   - Verify that samples with LITE files are not in the list
   - Confirm getfastq only processes full SRA files

3. **Monitor extraction**:
   - Check that TMPDIR is set correctly during execution
   - Verify vdb-config repository path is configured
   - Monitor disk space usage in `.tmp/fasterq-dump` and `.tmp/vdb`

## Files Changed

1. `src/metainformant/rna/amalgkit.py` - TMPDIR and VDB_CONFIG environment setup
2. `src/metainformant/rna/workflow.py` - vdb-config initialization, LITE filtering, early exit logic
3. `src/metainformant/rna/metadata_filter.py` - LITE file filtering logic
4. `docs/rna/amalgkit/steps/04_getfastq.md` - Updated troubleshooting documentation

## Expected Behavior After Fixes

1. **Workflow Initialization**:
   - vdb-config repository path is automatically configured to `.tmp/vdb`
   - Log message: "Configured vdb-config repository path: /path/to/repo/.tmp/vdb"

2. **Metadata Filtering** (after select step):
   - LITE files are automatically excluded
   - Log message: "Excluding X LITE SRA files (metadata-only, no sequence data)"

3. **Getfastq Step**:
   - TMPDIR, TEMP, TMP set to `.tmp/fasterq-dump`
   - VDB_CONFIG set to `.tmp/vdb`
   - FASTQ files are successfully extracted
   - No "disk-limit exceeded" errors

4. **Validation** (after getfastq):
   - FASTQ files are detected and counted
   - If 0 FASTQ files: workflow stops with clear error message
   - If FASTQ files exist: workflow continues to next step

5. **Pre-Step Validation**:
   - Before integrate: checks for FASTQ files
   - Before quant: checks for FASTQ files and tools
   - Before merge: checks for quantification output and R

## Verification Commands

```bash
# Check vdb-config repository path
vdb-config -i  # Should show repository path set to .tmp/vdb

# Check TMPDIR during execution
echo $TMPDIR  # Should show .tmp/fasterq-dump

# Check for LITE files in metadata
python3 -c "
import pandas as pd
df = pd.read_csv('output/amalgkit/pbarbatus_test5/work/metadata/metadata_selected.tsv', sep='\t')
lite_count = df['AWS_Link'].astype(str).str.contains('lite', case=False, na=False).sum()
print(f'LITE files in selected metadata: {lite_count}')  # Should be 0
"

# Check FASTQ files after getfastq
find output/amalgkit/pbarbatus_test5/fastq/getfastq -name "*.fastq*" | wc -l  # Should be > 0
```

## Summary

All critical fixes have been implemented:
- ✅ TMPDIR and VDB_CONFIG environment variables set automatically
- ✅ vdb-config repository path configured at workflow initialization
- ✅ LITE file filtering prevents downloading metadata-only files
- ✅ Early exit logic stops workflow when extraction fails
- ✅ Pre-step validation catches missing prerequisites early

**Status**: Ready for testing. Re-run getfastq step to verify FASTQ extraction works.

