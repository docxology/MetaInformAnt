# AmalgKit Storage Management Action Plan

**Generated:** 2024-10-31
**Location:** `/Users/4d/Documents/GitHub/metainformant/output/amalgkit/`

## Executive Summary

Total storage in amalgkit output: **~91 GB**
- **P. barbatus**: 55 GB (54 GB in fastq)
- **C. floridanus**: 35 GB (34 GB in fastq)
- **M. pharaonis**: 825 MB
- **S. invicta**: 908 MB

### Immediate Action Available
**19.83 GB** can be reclaimed immediately by deleting quantified SRA files in P. barbatus.

### Total Potential Savings
**52.5 GB** can be reclaimed after completing all pending quantifications.

---

## Detailed Breakdown

### P. barbatus (55 GB)

#### ✅ Ready for Cleanup (18 samples)
These samples have completed quantification and SRA files can be safely deleted:

| Sample | Size | Quant Status |
|--------|------|--------------|
| SRR14740487 | 7.1G | ✓ Complete |
| SRR14740535 | 7.0G | ✓ Complete |
| SRR14740551 | 6.5G | ✓ Complete |
| SRR14740519 | 4.4G | ✓ Complete |
| SRR1817186 | 2.0G | ✓ Complete |
| SRR14740520 | 1.9G | ✓ Complete |
| SRR1817187 | 1.9G | ✓ Complete |
| SRR1817188 | 1.9G | ✓ Complete |
| SRR1817185 | 1.7G | ✓ Complete |
| SRR1817189 | 1.4G | ✓ Complete |
| SRR14740543 | 1.3G | ✓ Complete |
| SRR14740495 | 934M | ✓ Complete |
| SRR1817183 | 666M | ✓ Complete |
| SRR14740511 | 582M | ✓ Complete |
| SRR14740503 | 390M | ✓ Complete |
| SRR14740488 | 313M | ✓ Complete |
| SRR14740536 | 279M | ✓ Complete |
| SRR14740527 | 8.4M | ✓ Complete |

**Reclaimable:** 19.83 GB

#### ⚠️ Needs Quantification (2 samples)
These samples have SRA files but no quantification output:

| Sample | Size |
|--------|------|
| DRR029869 | 363M |
| DRR048431 | 226M |

**Additional reclaimable after quantification:** 0.57 GB

**Sample list saved to:** `pbarbatus_unquantified.txt`

---

### C. floridanus (35 GB)

#### ⚠️ All Samples Need Quantification (18 samples)

**None** of the C. floridanus samples have been quantified yet. All 18 samples with SRA files require processing:

| Sample | Size | Priority |
|--------|------|----------|
| SRR22031377 | 7.4G | High (largest) |
| SRR32143976 | 4.8G | High |
| SRR25496822 | 4.0G | High |
| SRR32143984 | 3.8G | Medium* |
| SRR25496830 | 3.4G | Medium |
| SRR5931476 | 1.8G | Medium |
| SRR10960077 | 1.6G | Medium |
| SRR5931475 | 1.1G | Low |
| SRR2060722 | 1.0G | Low |
| SRR5931468 | 1.0G | Low |
| SRR32143992 | 1.0G | Low |
| SRR32143986 | 795M | Low |
| SRR25496798 | 789M | Low |
| SRR25496806 | 553M | Low |
| SRR32143985 | 801M | Low |
| SRR25496814 | 326M | Low |
| SRR25496790 | 257M | Low |
| SRR22031369 | 12M | Low |

\* SRR32143984 has FASTQ files extracted but no quantification output

**Total reclaimable after quantification:** 32.12 GB

**Sample list saved to:** `cfloridanus_unquantified.txt`

---

## Action Scripts and Reports

Storage management utilities are available in the repository:

### 1. `output/amalgkit/storage_analysis.md`
Complete storage analysis report with detailed breakdown by species and sample.

### 2. `scripts/rna/cleanup_quantified_sra.sh` ✨
**Safe deletion of quantified SRA files**

```bash
# Preview what will be deleted (dry run)
bash scripts/rna/cleanup_quantified_sra.sh

# Execute cleanup
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

Features:
- Verifies quantification completion before deletion
- Provides detailed logging of all operations
- Shows space reclaimed after completion
- Skips unquantified samples with warnings

### 3. `scripts/rna/list_unquantified.sh`
**Generate lists of samples needing quantification**

```bash
# Run to create/update sample lists
bash scripts/rna/list_unquantified.sh
```

Outputs:
- `output/amalgkit/pbarbatus_unquantified.txt` (2 samples)
- `output/amalgkit/cfloridanus_unquantified.txt` (18 samples)
- Summary of sizes and counts

---

## Recommended Workflow

### Phase 1: Immediate Cleanup (P. barbatus)
```bash
cd /Users/4d/Documents/GitHub/metainformant

# Review what will be deleted
bash scripts/rna/cleanup_quantified_sra.sh

# Execute cleanup to reclaim 19.83 GB
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

**Expected result:** Free up ~20 GB immediately

### Phase 2: Complete P. barbatus
```bash
# Quantify remaining 2 samples
# Use existing pbarbatus workflow with:
cat output/amalgkit/pbarbatus_unquantified.txt
# DRR029869
# DRR048431

# After quantification completes, run cleanup again
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

**Expected result:** Additional 0.57 GB freed

### Phase 3: Process C. floridanus
```bash
# Review sample list
cat output/amalgkit/cfloridanus_unquantified.txt

# Run amalgkit quantification workflow for C. floridanus
# (Use existing multi-species or single-species scripts)

# After each batch completes, run cleanup
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

**Expected result:** Up to 32 GB freed after all quantifications

---

## Workflow Integration Recommendations

### Automatic Cleanup Strategy

Consider modifying the amalgkit workflow to automatically clean up after quantification:

1. **After FASTQ extraction:** Delete SRA file immediately
2. **After quantification:** Delete FASTQ files, keep only quant outputs
3. **Retention policy:** Keep quantification outputs permanently, discard intermediate files

### Storage Efficiency

Current workflow appears to retain all intermediate files:
- SRA files (original downloads)
- FASTQ files (extracted from SRA)
- Quantification outputs (final results)

**Recommendation:** Implement progressive cleanup to retain only final outputs.

### Monitoring

Create a cron job or scheduled task to periodically run:
```bash
bash scripts/rna/list_unquantified.sh
bash scripts/rna/cleanup_quantified_sra.sh --execute
```

This ensures continuous cleanup as quantifications complete.

---

## Space Recovery Timeline

| Phase | Action | Space Recovered | Cumulative |
|-------|--------|----------------|------------|
| 0 (Current) | - | 0 GB | 91 GB used |
| 1 (Immediate) | Delete P. barbatus quantified SRA | 19.83 GB | 71 GB used |
| 2 (Short-term) | Complete P. barbatus, cleanup | 0.57 GB | 70 GB used |
| 3 (Medium-term) | Quantify & cleanup C. floridanus | 32.12 GB | 38 GB used |

**Final state:** ~38 GB (42% of current usage)

---

## Notes

- All SRA files should be considered **ephemeral** once quantification is complete
- FASTQ files can also be deleted after successful quantification
- Only quantification outputs (`abundance.tsv`, `abundance.h5`, `run_info.json`) need long-term retention
- The analysis includes the merge step which creates final expression matrices

## Support Files and Scripts

### Scripts (in `scripts/rna/`)
- `cleanup_quantified_sra.sh` - Safe deletion utility
- `list_unquantified.sh` - Sample listing utility

### Reports (in `output/amalgkit/`)
- `storage_analysis.md` - Detailed size breakdown
- `pbarbatus_unquantified.txt` - P. barbatus samples needing quantification
- `cfloridanus_unquantified.txt` - C. floridanus samples needing quantification
- `action_plan.md` - This document

---

**Next Step:** Review this plan and execute Phase 1 cleanup to immediately reclaim ~20 GB.

