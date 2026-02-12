# Workflow Verification Checklist

## Install Optional Dependencies

To remove warnings about missing optional dependencies:

```bash
# Install ncbi-datasets-pylib (biopython is already in main dependencies)
uv pip install ncbi-datasets-pylib --python /tmp/metainformant_venv/bin/python3
```

## Convert Existing SRA Files (Parallel Processing)

The conversion script now supports parallel processing with detailed progress logging:

```bash
/tmp/metainformant_venv/bin/python3 scripts/rna/convert_existing_sra.py
```

**Features**:
- **Parallel processing**: Up to 10 samples converted simultaneously
- **Progress tracking**: Shows elapsed time, estimated remaining time, file counts
- **Detailed logging**: Each conversion shows SRA size, output size, and timing
- **Thread-safe**: Safe concurrent execution with proper locking
- **Timeout protection**: Automatic timeouts prevent stuck processes (30min-6 hours based on file size)

## Quick Status Check

Run these commands to verify everything is working:

```bash
# 1. Check workflow status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status

# 2. Check active processes
ps aux | grep -E "amalgkit|fasterq-dump|run_workflow" | grep -v grep

# 3. Count files
echo "FASTQ files: $(find output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq -name '*.fastq*' 2>/dev/null | wc -l)"
echo "SRA files: $(find output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq -name '*.sra' 2>/dev/null | wc -l)"
echo "Quantified: $(find output/amalgkit/pogonomyrmex_barbatus/quant -name 'abundance.tsv' 2>/dev/null | wc -l)"

# 4. Verify wrapper exists
test -f output/amalgkit/pogonomyrmex_barbatus/fastq/temp/fasterq-dump && echo "✅ Wrapper exists" || echo "❌ Wrapper missing"

# 5. Check recent logs
tail -50 output/amalgkit/pogonomyrmex_barbatus/logs/workflow_final_*.log 2>/dev/null | grep -E "Downloading|Quantifying|finished|✅|⚠️|❌" | tail -15

# 6. Check for errors
tail -100 output/amalgkit/pogonomyrmex_barbatus/logs/workflow_final_*.log 2>/dev/null | grep -i "error\|failed\|exception" | tail -5
```

## Code Verification (Already Complete)

✅ **Automatic SRA-to-FASTQ Conversion**
- Location: `src/metainformant/rna/engine/sra_extraction.py::run()` (lines 274-295)
- Status: Implemented and working

✅ **Wrapper Script for amalgkit**
- Location: `src/metainformant/rna/engine/pipeline.py::_download_worker()` (lines 259-270)
- Status: Implemented and working
- Wrapper file: `output/amalgkit/pogonomyrmex_barbatus/fastq/temp/fasterq-dump` exists

✅ **Direct Binary Detection**
- Location: `src/metainformant/rna/engine/sra_extraction.py::convert_sra_to_fastq()` (lines 578-614)
- Status: Implemented and working

✅ **Documentation**
- User docs: `docs/rna/amalgkit/steps/getfastq.md` - Updated
- API docs: `docs/rna/API.md` - Updated
- Technical docs: `docs/rna/amalgkit/steps/EXTRACTION_FIXES.md` - Complete

✅ **Configuration**
- `config/amalgkit/amalgkit_pbarbatus.yaml` - AWS-only enabled

✅ **No Linter Errors**
- All code passes linting

## Expected Behavior

1. **Downloads**: 10 parallel workers downloading samples
2. **Extraction**: Automatic conversion of SRA to FASTQ
3. **Quantification**: Proceeds once FASTQ files are available
4. **Cleanup**: FASTQ files deleted after quantification

## What to Look For

- ✅ Active `amalgkit getfastq` processes
- ✅ SRA files being downloaded
- ✅ FASTQ files being created
- ✅ Quantification running
- ✅ No "disk-limit exceeded" errors
- ✅ No duplicate flag errors
- ✅ Wrapper script exists and is executable

## If Issues Found

1. Check logs: `tail -f output/amalgkit/pogonomyrmex_barbatus/logs/workflow_final_*.log`
2. Verify wrapper: `ls -l output/amalgkit/pogonomyrmex_barbatus/fastq/temp/fasterq-dump`
3. Check processes: `ps aux | grep amalgkit`
4. Review documentation: `docs/rna/amalgkit/steps/EXTRACTION_FIXES.md`

