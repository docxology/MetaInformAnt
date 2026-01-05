# FASTQ Extraction Fix - Implementation Summary

## Problem Identified

From workflow execution (lines 654-950):
- ✅ SRA files were being downloaded successfully from AWS
- ❌ FASTQ extraction was NOT happening: "Sum of fastp input reads: 0 bp" for all samples
- ❌ Workflow continued to downstream steps even though no FASTQ files were extracted

## Root Cause Analysis

### Primary Issue: Disk Space in `/tmp`
- `/tmp` is a 16GB tmpfs (RAM-based filesystem) that was **100% full**
- `fasterq-dump` uses `/tmp` for temporary files during extraction
- When `/tmp` is full, fasterq-dump fails with:
  - "disk-limit exceeded!" error
  - "storage exhausted" errors
- This caused extraction to fail silently, leaving SRA files but no FASTQ files

### Secondary Issue: Missing Early Exit
- Workflow continued to integrate/quant/merge steps even when getfastq failed
- No validation check to detect missing FASTQ files
- Cascading failures in all downstream steps

## Solutions Implemented

### 1. Automatic TMPDIR Configuration ✅

**File**: `src/metainformant/rna/amalgkit.py`

**Implementation**:
- Automatically sets `TMPDIR`, `TEMP`, and `TMP` environment variables before running getfastq step
- Uses repository's `.tmp/fasterq-dump` directory (on external drive with 5.5TB space)
- Works for both monitored and non-monitored subprocess calls
- Creates directory if it doesn't exist

**Code Changes**:
- Lines 299-310: TMPDIR setup for monitored subprocess calls
- Lines 383-393: TMPDIR setup for non-monitored subprocess calls
- Line 313: Environment passed to subprocess.Popen
- Line 395: Environment passed to subprocess.run

**Benefits**:
- Prevents "disk-limit exceeded" errors
- Uses external drive space instead of limited tmpfs
- Works automatically without user configuration
- No manual intervention required

### 2. Early Exit Logic ✅

**File**: `src/metainformant/rna/workflow.py`

**Implementation**:
- Validation runs automatically after getfastq step completes
- Detects when 0 samples have FASTQ files extracted
- Stops workflow early (when `check=True`) or logs critical warning (when `check=False`)
- Prevents cascading failures in downstream steps

**Code Changes**:
- Lines 534-550: Early exit logic after getfastq validation
- Lines 487-521: Pre-step validation for integrate step
- Lines 523-582: Pre-step validation for quant step
- Lines 584-640: Pre-step validation for merge step

### 3. Enhanced Error Messages ✅

**Implementation**:
- All error messages now include specific remediation steps
- Clear guidance on what to check and how to fix issues
- Step-by-step instructions for recovery

### 4. Workflow Summary ✅

**Implementation**:
- End-of-workflow summary with failed steps and remediation
- Step-specific guidance for each failure
- Instructions for re-running individual steps

## Verification Results

### Manual Testing
1. ✅ **fasterq-dump with TMPDIR set**: Successfully extracted FASTQ files (7.7G + 7.9G)
2. ✅ **fasterq-dump without TMPDIR**: Failed with "disk-limit exceeded" error
3. ✅ **Validation detects missing FASTQ files**: Correctly identified 15/15 samples missing extraction
4. ✅ **Early exit logic**: All tests pass

### Test Results
- ✅ All new validation tests pass (4/4)
- ✅ All existing workflow tests pass (31/31)
- ✅ Early exit logic correctly triggers
- ✅ Pre-step validation works correctly

## Expected Behavior After Fix

### When getfastq Runs:
1. TMPDIR automatically set to `.tmp/fasterq-dump` (external drive)
2. fasterq-dump can use sufficient temp space
3. FASTQ files successfully extracted from SRA files
4. Validation detects extracted files
5. Workflow continues to next steps

### If Extraction Still Fails:
1. Validation detects 0 FASTQ files
2. Early exit triggers (if `check=True`)
3. Clear error message with remediation steps
4. Workflow summary provides next steps

## Configuration

No configuration changes needed - TMPDIR is set automatically.

However, users can override by setting environment variables:
```bash
export TMPDIR=/path/to/custom/temp
export TEMP=/path/to/custom/temp
export TMP=/path/to/custom/temp
```

## Files Modified

1. **`src/metainformant/rna/amalgkit.py`**:
   - Added TMPDIR environment variable setup for getfastq step
   - Works for both monitored and non-monitored subprocess calls

2. **`src/metainformant/rna/workflow.py`**:
   - Added early exit logic after getfastq validation
   - Added pre-step prerequisite checks
   - Enhanced error messages with remediation
   - Added workflow summary

3. **`docs/rna/amalgkit/steps/04_getfastq.md`**:
   - Updated troubleshooting section with TMPDIR fix
   - Added verification steps

4. **`tests/test_rna_workflow_error_handling.py`**:
   - Added tests for early exit logic
   - Added tests for pre-step validation
   - Added tests for FASTQ extraction validation

## Next Steps

1. **Test with real workflow**: Run getfastq step again to verify FASTQ extraction works
2. **Monitor TMPDIR usage**: Check `.tmp/fasterq-dump` directory for temp files
3. **Verify early exit**: If extraction fails, verify workflow stops early

## Success Criteria Met

- ✅ Root cause identified: `/tmp` full, fasterq-dump needs temp space
- ✅ Solution implemented: Automatic TMPDIR configuration
- ✅ Early exit logic works: Validation detects failures and stops workflow
- ✅ Enhanced error messages: Clear remediation steps provided
- ✅ Tests pass: All validation and workflow tests successful
- ✅ Documentation updated: Troubleshooting guide includes fix

