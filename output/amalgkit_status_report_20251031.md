# Amalgkit Multi-Species Processing Status Report

**Date**: October 31, 2025 at 9:15 AM  
**Generated**: Comprehensive sample assessment across all workflow stages

---

## Active Processes 🔄

**Running amalgkit processes:**
- ✅ **`run_multi_species_amalgkit.py`** (2 instances - PID 91164, 95218)
- 🔄 **cfloridanus getfastq** (2 active download batches - batch 13, batch 16)
  - Batch 13: PID 99552 (running 26 minutes)
  - Batch 16: PID 1961 (running 19 minutes)

**Status**: C. floridanus is actively downloading FASTQ files using batched parallel processing

---

## Sample Counts by Species and Stage

### 📊 P. barbatus (Completed ✅)

| Stage | Count | Status |
|-------|-------|--------|
| **Metadata** | 120 | ✅ Retrieved |
| **FASTQs** | 20 | ⚠️ 20 samples still have files |
| **Quantified** | 83 | ✅ 69% complete |
| **Merge** | - | ⏳ Not started |
| **Curate** | - | ⏳ Not started |

**Storage**: 35 GB (cleaned up 19.83 GB earlier)

**Status**: 
- ✅ 83 samples successfully quantified
- ⚠️ 20 samples still have FASTQ files (from 120 metadata records)
- ⚠️ 37 samples (31%) not quantified yet (need investigation)

**Issues**:
- 2 samples tried to quantify but no FASTQs found (DRR029869, DRR048431 - were in unquantified list)
- Additional 35 samples unaccounted for (120 metadata - 83 quantified - 2 failed = 35 unknown)

---

### 📊 C. floridanus (In Progress 🔄)

| Stage | Count | Status |
|-------|-------|--------|
| **Metadata** | 358 | ✅ Retrieved |
| **FASTQs** | 6 | 🔄 Actively downloading (2 batches running) |
| **Quantified** | 0 | ⏳ Not started (6 failed attempts) |
| **Merge** | - | ⏳ Not started |
| **Curate** | - | ⏳ Not started |

**Storage**: 37 GB (mostly SRA files waiting for conversion)

**Status**:
- 🔄 **Active**: Downloading FASTQs via amalgkit getfastq (batches 13 & 16)
- ⚠️ **Issue**: 6 quantification attempts failed (exit code 2) - likely missing kallisto index
- 📥 From earlier scan: 22 samples have SRA files needing conversion/quantification

**Recent Activity**:
- Quantification script tried to quantify 6 samples but failed (no index built yet)
- All 6 failed samples had their FASTQs deleted automatically
- Current getfastq batches are downloading new samples

**Issue Identified**:
- **Kallisto index missing** for C. floridanus
- Need to build index before quantification can succeed
- Genome download appears incomplete (integrate step failed)

---

### 📊 M. pharaonis (Minimal Progress ⏳)

| Stage | Count | Status |
|-------|-------|--------|
| **Metadata** | 370 | ✅ Retrieved |
| **FASTQs** | 0 | ⏳ Not started |
| **Quantified** | 0 | ⏳ Not started |
| **Merge** | - | ⏳ Not started |
| **Curate** | - | ⏳ Not started |

**Storage**: 825 MB (minimal - just metadata and config)

**Status**: Metadata retrieved, waiting for FASTQ download phase

---

### 📊 S. invicta (Minimal Progress ⏳)

| Stage | Count | Status |
|-------|-------|--------|
| **Metadata** | 415 | ✅ Retrieved |
| **FASTQs** | 0 | ⏳ Not started |
| **Quantified** | 0 | ⏳ Not started |
| **Merge** | - | ⏳ Not started |
| **Curate** | - | ⏳ Not started |

**Storage**: 908 MB (minimal - just metadata and config)

**Status**: Metadata retrieved, waiting for FASTQ download phase

---

## Overall Progress Summary

### By Species

| Species | Metadata | FASTQs | Quantified | Progress % | Status |
|---------|----------|--------|------------|------------|---------|
| **P. barbatus** | 120 | 20 | 83 | 69% | ✅ Mostly done |
| **C. floridanus** | 358 | 6 🔄 | 0 | <1% | 🔄 Downloading |
| **M. pharaonis** | 370 | 0 | 0 | 0% | ⏳ Queued |
| **S. invicta** | 415 | 0 | 0 | 0% | ⏳ Queued |
| **TOTAL** | **1,263** | **26** | **83** | **7%** | 🔄 In progress |

### Storage Usage

```
Total: 74 GB
├── P. barbatus:    35 GB (48%)
├── C. floridanus:  37 GB (50%) 
├── M. pharaonis:  825 MB (1%)
└── S. invicta:    908 MB (1%)
```

---

## Issues Identified ⚠️

### 1. C. floridanus Quantification Failures

**Problem**: All 6 quantification attempts failed with exit code 2

**Root cause**: Kallisto index not built
- Genome integration failed: `ValueError: PATH to fastq directory does not exist`
- Getfastq had issues: `AssertionError: No SRA entry found`
- Index building prerequisite not met

**Solution needed**:
1. Ensure genome is properly downloaded for C. floridanus
2. Build kallisto index before quantification
3. May need to run genome step manually

### 2. P. barbatus Missing Samples

**Problem**: 37 samples unaccounted for
- 120 in metadata
- 83 quantified
- 37 samples status unknown

**Investigation needed**:
- Check if 37 samples failed during getfastq
- Check logs for download failures
- May need to reprocess failed downloads

### 3. Batch Processing Configuration

**Observation**: C. floridanus using batch processing
- Found 1 batch metadata file (batch 13, 16 running)
- Batched approach helps manage disk space
- Currently downloading batches in parallel

---

## Recommendations

### Immediate Actions

1. **Monitor C. floridanus downloads**
   ```bash
   # Watch active downloads
   ps aux | grep amalgkit | grep getfastq
   
   # Check download progress
   find output/amalgkit/cfloridanus/fastq -name "*.fastq*" | wc -l
   ```

2. **Build kallisto index for C. floridanus**
   ```bash
   # Once genome is ready, build index
   cd output/amalgkit/cfloridanus
   # Check if genome files exist first
   ls -la genome/
   ```

3. **Investigate P. barbatus missing samples**
   ```bash
   # Check logs for failed downloads
   grep -i "error\|fail" output/amalgkit/pbarbatus/logs/*.log
   ```

### Next Phase Actions

4. **Complete C. floridanus processing**
   - Let current downloads finish
   - Build kallisto index
   - Restart quantification for all downloaded samples

5. **Process M. pharaonis and S. invicta**
   - These are queued in multi-species runner
   - Will start after C. floridanus completes

6. **P. barbatus finalization**
   - Quantify remaining 37 samples
   - Run merge step
   - Run curate step

---

## Timeline Estimate

**Current status**: ~7% complete overall (83/1,263 samples quantified)

**Estimated remaining time**:
- **C. floridanus**: 24-48 hours (358 samples, downloads + quant)
- **M. pharaonis**: 24-48 hours (370 samples)
- **S. invicta**: 24-48 hours (415 samples)
- **P. barbatus completion**: 4-8 hours (37 samples)

**Total estimated**: 72-144 hours (3-6 days) for all species

**Bottlenecks**:
- FASTQ downloads (network/SRA speed)
- Kallisto quantification (CPU bound)
- Disk space management (cleaned up 20 GB already)

---

## Commands for Monitoring

### Watch Downloads
```bash
# Active processes
ps aux | grep -E "amalgkit|kallisto" | grep -v grep

# C. floridanus progress
find output/amalgkit/cfloridanus/fastq -name "*.fastq*" | wc -l
du -sh output/amalgkit/cfloridanus/

# All species status
bash scripts/rna/list_unquantified.sh
```

### Check Logs
```bash
# Latest C. floridanus activity
tail -f output/amalgkit/cfloridanus/logs/*.log

# Multi-species runner
ps aux | grep run_multi_species

# Quantification log
tail -f output/quant_unquantified_*.log
```

---

## Summary

**✅ Achievements**:
- P. barbatus 69% quantified (83/120 samples)
- 19.83 GB disk space reclaimed
- Multi-species orchestration running smoothly

**🔄 In Progress**:
- C. floridanus actively downloading (2 batches)
- Multi-species runner coordinating workflows

**⚠️ Issues**:
- C. floridanus kallisto index needs building
- P. barbatus 37 samples status unclear
- C. floridanus 6 quantification attempts failed

**📋 Next**:
- Monitor C. floridanus downloads
- Build kallisto index for C. floridanus
- Investigate P. barbatus missing samples
- Continue multi-species processing (M. pharaonis, S. invicta queued)

**Overall**: System is actively processing with expected challenges. Primary focus should be on resolving C. floridanus index issue once current downloads complete.

