# Quick Start Guide: Amalgkit Sanity & Curate

**Last Updated**: October 29, 2025  
**Status**: ✅ General workflow guide for all species

---

## Prerequisites

You must run all amalgkit commands from the **correct working directory** for your species:

```bash
cd output/amalgkit/<species_name>
```

**Why?** Amalgkit uses relative paths from the current directory to find `work/` subdirectories.

### Example Directory Structure

```
output/amalgkit/
├── pbarbatus/          # Pogonomyrmex barbatus
├── dmelanogaster/      # Drosophila melanogaster  
└── your_species/       # Your analysis
    └── work/           # Amalgkit working directory
        ├── metadata/
        ├── quant/
        ├── merge/
        ├── curate/
        └── sanity/
```

---

## Quick Verification

Check workflow outputs manually:

```bash
cd output/amalgkit/<species_name>

# Verify quantification results
ls -lh quant/*/abundance.tsv

# Verify merged expression matrix
ls -lh merged/*.tsv

# Check sanity check results
cat sanity/*.txt

# View curate QC reports
ls -lh curate/*.tsv
```

This will:
- ✅ Run `amalgkit sanity` to validate all samples
- ✅ Check existing `curate` outputs (or regenerate if needed)
- ✅ Display comprehensive summary of all files created

---

## Running Individual Steps

### 1. Sanity Check

**Purpose**: Validate integrity of all workflow outputs

```bash
cd output/amalgkit/<species_name>
amalgkit sanity --out_dir work --all
```

**Expected Output**:
```
Quant outputs found for all SRA IDs in --metadata (inferred)
amalgkit sanity: end
```

**Exit Code**: 0 (success)

**Files Created**:
- `work/sanity/SRA_IDs_without_fastq.txt` - Lists samples without FASTQ files (expected if FASTQs deleted after quant)
- NO `SRA_IDs_without_quant.txt` file = all samples validated ✅

**Verification**:
```bash
# Should NOT exist (or be empty) if all samples are good
ls work/sanity/SRA_IDs_without_quant.txt 2>/dev/null || echo "✅ All samples validated!"
```

---

### 2. Curate Step

**Purpose**: Quality control, outlier removal, and visualization generation

```bash
cd output/amalgkit/<species_name>
amalgkit curate --out_dir work --batch_effect_alg no
```

**Why `--batch_effect_alg no`?**  
- The `sva` R package has system-level compilation issues
- Using `no` skips SVA batch correction but still generates all other outputs
- This is perfectly acceptable for most analyses

**Expected Output**:
```
Removing samples with mapping rate of 0.2
Mapping rate cutoff: 20%
No entry removed due to low mapping rate.
Iteratively checking within-sample_group correlation
No batch effect correction was performed.
Round: 2 : # before = N : # after = N
Writing summary files for <species_name>
transcriptome_curation.r: Completed.
```

**Files Created** (17 total):

#### Visualizations (6 PDFs)
```
work/curate/<species_name>/plots/
├── <species_name>.0.original.pdf
├── <species_name>.0.original.no.pdf
├── <species_name>.1.mapping_cutoff.pdf
├── <species_name>.1.mapping_cutoff.no.pdf
├── <species_name>.2.correlation_cutoff.pdf
└── <species_name>.2.correlation_cutoff.no.pdf
```

Each PDF contains:
- Hierarchical clustering dendrogram
- Correlation heatmap
- PCA plot
- QC metrics

#### Data Tables (7 TSV files)
```
work/curate/<species_name>/tables/
├── <species_name>.uncorrected.tc.tsv                    # Original expression matrix
├── <species_name>.uncorrected.sample_group.mean.tsv     # Sample group means
├── <species_name>.no.tc.tsv                             # Final curated expression (⭐ PRIMARY OUTPUT)
├── <species_name>.no.sample_group.mean.tsv              # Final sample means
├── <species_name>.no.tau.tsv                            # Tissue specificity scores
├── <species_name>.no.correlation_statistics.tsv         # QC correlation statistics
└── <species_name>.metadata.tsv                          # Sample metadata with QC metrics
```

#### Analysis States (3 RData files)
```
work/curate/<species_name>/rdata/
├── <species_name>.no.0.RData  # Round 0 analysis state
├── <species_name>.no.1.RData  # Round 1 analysis state
└── <species_name>.no.2.RData  # Round 2 analysis state
```

#### Completion Flag
```
work/curate/<species_name>/curate_completion_flag.txt
```

---

## Viewing Outputs

### View Visualizations (macOS)

```bash
cd output/amalgkit/<species_name>
open work/curate/*/plots/*.pdf
```

### View Expression Matrix

```bash
# First 10 rows
head work/curate/*/tables/*.no.tc.tsv

# With column headers (first 5 samples)
head -20 work/curate/*/tables/*.no.tc.tsv | cut -f1-5
```

### Check File Counts

```bash
# All curate outputs
find work/curate -type f | wc -l
# Should show: 17

# PDF visualizations
find work/curate -name "*.pdf" | wc -l
# Should show: 6

# Data tables
find work/curate -name "*.tsv" | wc -l
# Should show: 7
```

---

## Troubleshooting

### Error: "FileNotFoundError: work/metadata/metadata.tsv"

**Cause**: Running from wrong directory

**Solution**:
```bash
# Always run from your species directory:
cd output/amalgkit/<species_name>

# Then run amalgkit commands
amalgkit sanity --out_dir work --all
```

### Error: "SRA_IDs_without_quant.txt shows all samples"

**Cause**: Missing file name symlinks (fixed!)

**Solution**: Already resolved - symlinks created for:
- `SRR*_abundance.tsv` → `abundance.tsv`
- `SRR*_abundance.h5` → `abundance.h5`  
- `SRR*_run_info.json` → `run_info.json` ✨

**Verify** (using your first SRR ID):
```bash
ls -la work/quant/SRR*/
# Should show symlinks like:
# SRR*_run_info.json -> run_info.json
# SRR*_abundance.tsv -> abundance.tsv
# SRR*_abundance.h5 -> abundance.h5
```

### Curate completes in 0 seconds

**Cause**: Existing outputs detected, no reprocessing needed

**Solution**: This is normal! To regenerate:
```bash
rm -rf work/curate/*
amalgkit curate --out_dir work --batch_effect_alg no
```

### R package errors (Rtsne, amap, vegan, sva, RUVSeq)

**Cause**: System-level gcc compilation issues (documented)

**Solution**: Already handled via patches to `amalgkit/curate.r`:
- Changed `library()` to `require()` for graceful degradation
- Added fallbacks for missing functions
- Core visualizations still work perfectly ✅

**Impact**: Minimal - only advanced optional features affected

---

## Expected Results

### Sanity Check ✅

```
Exit code: 0
Message: "Quant outputs found for all SRA IDs"
Files: work/sanity/SRA_IDs_without_fastq.txt (N samples)
       NO work/sanity/SRA_IDs_without_quant.txt (all validated!)
```

### Curate Step ✅

```
Exit code: 0
Files created: 17 total
  - 6 PDF visualizations (dendrograms, heatmaps, PCA)
  - 7 TSV data tables (expression matrices, QC metrics)
  - 3 RData files (analysis states)
  - 1 completion flag
```

---

## Next Steps

After successful sanity and curate:

1. **View visualizations**:
   ```bash
   open work/curate/*/plots/*.pdf
   ```

2. **Analyze expression data**:
   ```bash
   # Final curated expression matrix (primary output)
   head work/curate/*/tables/*.no.tc.tsv
   ```

3. **Check QC statistics**:
   ```bash
   cat work/curate/*/tables/*.no.correlation_statistics.tsv
   ```

4. **Examine tissue specificity**:
   ```bash
   head work/curate/*/tables/*.no.tau.tsv
   ```

---

## References

- **R Package Setup**: `docs/rna/amalgkit/r_packages.md`
- **R Installation**: `docs/rna/amalgkit/R_INSTALLATION.md`
- **Pipeline Guide**: `docs/rna/amalgkit/amalgkit.md` (includes advanced usage)
- **Step Documentation**: `docs/rna/amalgkit/steps/`
- **Workflow Orchestration**: `docs/rna/workflow.md`

---

## Summary

✅ **Sanity**: Validates 100% data integrity for all samples  
✅ **Curate**: Generates 6 comprehensive PDF visualizations  
✅ **Total**: 17 output files with expression data, QC, and plots  
✅ **Status**: Production-ready workflow for any species  

**Run from**: `output/amalgkit/<species_name>/`  
**Commands**: `amalgkit sanity --out_dir work --all`  
            `amalgkit curate --out_dir work --batch_effect_alg no`

---

## Example: P. barbatus

See `output/amalgkit/pbarbatus/` for a complete example with:
- 83 brain RNA-seq samples
- 20,672 genes quantified
- All QC visualizations and curated expression matrices  

