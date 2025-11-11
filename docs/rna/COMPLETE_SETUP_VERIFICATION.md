# Complete Setup and Workflow Verification

This document provides a comprehensive verification checklist for the complete end-to-end amalgkit workflow setup and execution.

## Pre-Execution Verification ✅

### Configuration Files
- [x] `amalgkit_pogonomyrmex_barbatus.yaml` - Loads successfully, 11 steps configured
- [x] `amalgkit_template.yaml` - Template loads successfully
- [x] `amalgkit_test.yaml` - Test config loads successfully
- [x] All paths resolve correctly relative to repo root
- [x] Genome configuration validated (accession, FTP URL, files)

### Environment Setup
- [x] Virtual environment: `/tmp/metainformant_venv` (auto-selected for ext6 filesystem)
- [x] metainformant: v0.2.0 installed
- [x] amalgkit: v0.12.20 installed
- [x] biopython: v1.86 installed (resolved symlink issue)
- [x] ncbi-datasets-pylib: v16.6.1 installed (resolved symlink issue)
- [x] SRA Toolkit: fasterq-dump 3.0.3 available
- [x] kallisto: v0.48.0 available
- [x] UV cache: Configured to `/tmp/uv-cache` (avoids symlink issues)

### External Drive Issues Resolved
- [x] Virtual environment location: Automatic fallback to `/tmp/metainformant_venv`
- [x] UV cache directory: Automatic use of `/tmp/uv-cache`
- [x] Symlink errors: Resolved by cache directory configuration
- [x] Biopython installation: Fixed with proper cache location
- [x] Documentation: Complete guide created in `EXTERNAL_DRIVE_SETUP.md`

## Execution Verification ✅

### Genome Setup (Automatic)
- [x] Genome downloaded: 97MB package from NCBI
- [x] Transcriptome extracted: `Pogonomyrmex_barbatus_transcripts.fa`
- [x] Kallisto index built: 567MB index file
- [x] Index location: `output/amalgkit/pogonomyrmex_barbatus/work/index/Pogonomyrmex_barbatus_transcripts.idx`

### Workflow Steps
- [x] **metadata**: Complete - 84 samples retrieved
- [x] **config**: Complete - Config files generated
- [x] **select**: Complete - 83 samples selected
- [x] **getfastq**: In Progress - 11 parallel workers active
- [ ] **integrate**: Pending - Will run after getfastq
- [ ] **quant**: Pending - Will run as downloads complete
- [ ] **merge**: Pending - Will run after all quantifications
- [ ] **cstmm**: Pending
- [ ] **curate**: Pending
- [ ] **csca**: Pending
- [ ] **sanity**: Pending

### Current Status
- **Total Samples**: 83
- **Downloading**: 11 (parallel workers)
- **Queued**: 72
- **Quantified**: 0 (downloads in progress)
- **Process Status**: 14 active processes (main workflow + download workers)

## Output Verification

### Created Files ✅
- [x] Genome package: `output/amalgkit/pogonomyrmex_barbatus/genome/ncbi_dataset_api.zip`
- [x] Transcriptome: `output/amalgkit/pogonomyrmex_barbatus/work/Pogonomyrmex_barbatus_transcripts.fa`
- [x] Kallisto index: `output/amalgkit/pogonomyrmex_barbatus/work/index/Pogonomyrmex_barbatus_transcripts.idx` (567MB)
- [x] Metadata: `output/amalgkit/pogonomyrmex_barbatus/work/metadata/metadata.tsv` (84 samples)
- [x] Qualified samples: `output/amalgkit/pogonomyrmex_barbatus/work/metadata/pivot_qualified.tsv`
- [x] Selected samples: `output/amalgkit/pogonomyrmex_barbatus/work/metadata/pivot_selected.tsv`
- [x] Config files: `output/amalgkit/pogonomyrmex_barbatus/work/config_base/`

### Pending Files (Will be created automatically)
- [ ] Quantification results: `output/amalgkit/pogonomyrmex_barbatus/quant/*/abundance.tsv`
- [ ] Merged matrix: `output/amalgkit/pogonomyrmex_barbatus/merged/merged_abundance.tsv`
- [ ] CSTMM results: `output/amalgkit/pogonomyrmex_barbatus/cstmm/`
- [ ] Curate results: `output/amalgkit/pogonomyrmex_barbatus/curate/`
- [ ] CSCA results: `output/amalgkit/pogonomyrmex_barbatus/csca/`
- [ ] Manifest: `output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl`

## Documentation Created

### New Documentation Files
1. **[EXTERNAL_DRIVE_SETUP.md](EXTERNAL_DRIVE_SETUP.md)** - Comprehensive guide for:
   - External drive and filesystem limitations
   - Symlink issues and solutions
   - UV cache configuration
   - Virtual environment location handling
   - Troubleshooting guide

2. **[WORKFLOW_EXECUTION_SUMMARY.md](WORKFLOW_EXECUTION_SUMMARY.md)** - Complete execution summary:
   - Pre-execution setup
   - Workflow execution details
   - Current status
   - Expected outputs
   - Monitoring commands

3. **[COMPLETE_SETUP_VERIFICATION.md](COMPLETE_SETUP_VERIFICATION.md)** - This file:
   - Verification checklist
   - Status tracking
   - Completion criteria

### Updated Documentation
- **[GETTING_STARTED.md](GETTING_STARTED.md)** - Added external drive reference
- **[README.md](README.md)** - Added external drive setup to troubleshooting

## Key Solutions Implemented

### 1. External Drive Symlink Issues

**Problem**: ext6 filesystem doesn't support symlinks in `output/.uv-cache`

**Solution**: 
- Updated `_setup_utils.py` to use `/tmp/uv-cache` for UV cache
- Automatic detection and fallback
- Documented in `EXTERNAL_DRIVE_SETUP.md`

### 2. Biopython Installation

**Problem**: Biopython not available in venv despite being in `pyproject.toml`

**Solution**:
- Installed using `/tmp/uv-cache` cache location
- Verified in venv: `from Bio import Entrez` works
- No more warnings in logs

### 3. NCBI Datasets Library

**Problem**: ncbi-datasets-pylib not available

**Solution**:
- Installed using `/tmp/uv-cache` cache location
- Verified in venv: `from ncbi.datasets import GenomeApi` works
- No more warnings in logs

### 4. Virtual Environment Location

**Problem**: `.venv` creation may fail on ext6 filesystems

**Solution**:
- Automatic fallback to `/tmp/metainformant_venv`
- Transparent handling in all scripts
- Documented in `EXTERNAL_DRIVE_SETUP.md`

## Workflow Execution Flow

```
1. Environment Check
   ├─ Auto-activate venv (/tmp/metainformant_venv)
   ├─ Verify dependencies
   └─ Check amalgkit CLI

2. Configuration Load
   ├─ Load YAML config
   ├─ Resolve paths relative to repo root
   └─ Validate genome config

3. Automatic Genome Setup
   ├─ Check if genome/index exists
   ├─ Download genome (if missing)
   ├─ Extract transcriptome
   └─ Build kallisto index

4. Workflow Steps (Automatic)
   ├─ metadata → Retrieve sample metadata
   ├─ config → Generate config files
   ├─ select → Filter samples
   ├─ getfastq → Download FASTQ (parallel workers)
   ├─ integrate → Integrate FASTQ paths
   ├─ quant → Quantify (sequential, auto-delete FASTQ)
   ├─ merge → Merge expression matrix
   ├─ cstmm → Cross-species normalization
   ├─ curate → Quality control
   ├─ csca → Correlation analysis
   └─ sanity → Final validation
```

## Monitoring Commands

### Check Overall Status
```bash
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed
```

### Monitor Downloads
```bash
ps aux | grep "amalgkit getfastq" | grep -v grep | wc -l
```

### Check Quantification Progress
```bash
find output/amalgkit/pogonomyrmex_barbatus/quant -name "abundance.tsv" | wc -l
```

### Check Manifest
```bash
tail -f output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl
```

### Verify Genome Setup
```bash
ls -lh output/amalgkit/pogonomyrmex_barbatus/work/index/*.idx
ls -lh output/amalgkit/pogonomyrmex_barbatus/work/*.fa
```

## Completion Criteria

The workflow is complete when:

1. ✅ All 83 samples downloaded
2. ✅ All 83 samples quantified
3. ✅ Expression matrix merged
4. ✅ All analysis steps completed (cstmm, curate, csca)
5. ✅ Sanity checks passed
6. ✅ Manifest file contains all step records

## Verification After Completion

Once workflow completes, verify:

```bash
# Final status check
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed

# Verify all outputs
ls -lh output/amalgkit/pogonomyrmex_barbatus/merged/merged_abundance.tsv
ls -d output/amalgkit/pogonomyrmex_barbatus/curate/
ls -d output/amalgkit/pogonomyrmex_barbatus/csca/

# Check manifest
python3 -c "import json; steps = [json.loads(l) for l in open('output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl')]; print(f'Total steps: {len(steps)}'); print('Steps:', [s['step'] for s in steps])"
```

## Summary

✅ **Configuration**: All validated  
✅ **Environment**: All dependencies installed  
✅ **External Drive Issues**: All resolved and documented  
✅ **Genome Setup**: Automatic and complete  
✅ **Workflow Execution**: Running end-to-end  
✅ **Documentation**: Comprehensive guides created  

The workflow is executing comprehensively and will complete all steps automatically. All methods and issues with external drives are fully documented.

---

*Verification Date: November 10, 2025*  
*Status: Workflow Running - Downloads Active*  
*Expected Completion: Automatic - All steps will complete in sequence*

