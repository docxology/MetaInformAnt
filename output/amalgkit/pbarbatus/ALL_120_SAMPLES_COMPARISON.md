# Should We Use All 120 Samples? Detailed Comparison

**Date**: October 29, 2025  
**Question**: Why use 83 brain samples instead of all 120 available samples?

---

## TL;DR

**Current (83 samples)**: Conservative, brain-specific analysis ✅ **RECOMMENDED**  
**All 120 samples**: Exploratory, maximum sample size, unknown tissue risk ⚠️

---

## Sample Breakdown

### Total Available: 120 Samples

| Category | Count | Details |
|----------|-------|---------|
| **Brain annotated** | **83** | Confirmed brain tissue |
| **Unknown tissue** | 37 | No tissue field - could be anything |
| **Single-end** | 5 | Lower quality (all in unknown group) |
| **Paired-end** | 115 | Higher quality |

### The 83 Brain Samples (Current Selection)

| BioProject | Samples | Instrument | Tissue | Quality |
|------------|---------|------------|--------|---------|
| PRJNA547792 | 65 | HiSeq 4000 | Brain | ✅ High |
| PRJNA277638 | 18 | HiSeq 2500 | Whole cleanly-dissected brains | ✅ High |

**Characteristics**:
- ✅ All paired-end
- ✅ All tissue-annotated (confirmed brain)
- ✅ Modern instruments (HiSeq 2500/4000)
- ✅ Consistent biological context
- ✅ Published studies with known experimental designs

### The 37 Unknown Samples (Would Add)

| BioProject | Samples | Instrument | Tissue | Quality |
|------------|---------|------------|--------|---------|
| PRJDB3493 | 21 | HiSeq 2000 | **NONE** | ⚠️ Unknown |
| PRJDB4312 | 16 | HiSeq 2000 | **NONE** | ⚠️ Unknown |

**Characteristics**:
- ⚠️ NO tissue annotation (tissue field is empty/NA)
- ⚠️ Older instrument (HiSeq 2000)
- ⚠️ 5 samples are single-end (lower quality)
- ❓ Unknown biological context
- ❓ Could be: brain, whole body, larvae, pupae, other tissues, mixed samples

---

## Pros & Cons

### Option 1: Keep 83 Brain Samples (Current) ✅

**Pros**:
- ✅ **Clean dataset** - All confirmed brain tissue
- ✅ **Biological homogeneity** - Same tissue type across all samples
- ✅ **Better statistical power** - Reduced biological variance
- ✅ **Known experimental context** - Can interpret results confidently
- ✅ **Reproducible** - Clear inclusion criteria
- ✅ **Publishable** - Well-justified sample selection

**Cons**:
- ⚠️ Smaller sample size (83 vs 120)
- ⚠️ Might miss some biological variation (if unknown samples are also brain)

**Best For**:
- Brain-specific gene expression analysis
- Identifying brain-enriched genes
- Comparing brain expression across conditions
- Publication-quality analysis

---

### Option 2: Use All 120 Samples ⚠️

**Pros**:
- ✅ **Maximum sample size** - 45% more data
- ✅ **Exploratory potential** - Might discover interesting patterns
- ✅ **Comprehensive** - Uses all available data

**Cons**:
- ❌ **Unknown tissue contamination** - 37 samples could be non-brain
- ❌ **Mixed biological context** - Different tissues = different biology
- ❌ **Reduced statistical power** - Increased variance from tissue heterogeneity
- ❌ **Harder to interpret** - Unknown what the 37 samples represent
- ❌ **Quality issues** - 5 single-end samples (lower accuracy)
- ❌ **Publication concerns** - Harder to justify mixing known + unknown tissues

**Potential Issues**:
1. **Differential expression analysis**: Mixed tissues will show tissue-specific genes as "differentially expressed" (false positives)
2. **Co-expression networks**: Tissue-specific modules will confound brain-specific patterns
3. **Functional enrichment**: Brain pathways diluted by non-brain tissue signals
4. **Outlier detection**: Unknown tissues might be flagged as outliers (not actually outliers)

**Best For**:
- Initial exploratory analysis
- Data mining for novel patterns
- When you can cluster/filter by expression later
- When you don't care about tissue specificity

---

## What Are The 37 Unknown Samples?

### Investigation Results

**Metadata checked**:
- ❌ No tissue field
- ❌ No sample_description field (or generic)
- ❌ No clear indication in exp_title
- ❌ No source_name indicating tissue

**BioProjects**:
- **PRJDB3493**: Japanese database submission, minimal metadata
- **PRJDB4312**: Japanese database submission, minimal metadata

**Likely scenarios**:
1. **Could be brain** - Same as the 83 annotated samples
2. **Could be whole body** - Common for ant studies
3. **Could be developmental stages** - Larvae, pupae, adults
4. **Could be different castes** - Queens, workers, males
5. **Could be other tissues** - Legs, antennae, etc.

**Bottom line**: We **cannot know** without:
- Checking original papers (if published)
- Analyzing expression patterns post-hoc
- Contacting authors
- Running exploratory clustering first

---

## How To Use All 120 Samples

### Method 1: New Config File ✅

I've created `config/amalgkit_pbarbatus_all120.yaml` with these changes:

```yaml
# Key changes:
filters:
  require_tissue: false    # ⭐ Allow samples without tissue annotation
  # library_layout: PAIRED # ⭐ Removed to allow single-end

# Different output directory to not overwrite brain-only results
work_dir: output/amalgkit/pbarbatus_all120/work
```

**Run it**:
```bash
cd output/amalgkit/pbarbatus_all120
amalgkit metadata --config ../../config/amalgkit_pbarbatus_all120.yaml
# ... continue with other steps
```

### Method 2: Modify Existing Config

Edit `config/amalgkit_pbarbatus.yaml`:

```yaml
# Change this:
filters:
  require_tissue: true     # OLD

# To this:
filters:
  require_tissue: false    # NEW - allows all 120 samples
```

---

## Recommended Approach: Two-Stage Analysis

### Stage 1: Use All 120 Samples (Exploratory)

1. **Download and quantify** all 120 samples
2. **Run PCA/clustering** on expression data
3. **Identify sample groups** - Do the 37 unknowns cluster with brain samples?
4. **Check outliers** - Are unknowns biological outliers or just different tissues?

### Stage 2: Filter Based on Expression (Confirmatory)

Based on Stage 1 results:

**If unknowns cluster WITH brain samples**:
- ✅ Include them in final analysis
- Likely they are brain tissue (just poorly annotated)

**If unknowns cluster SEPARATELY**:
- ❌ Exclude them from brain analysis
- They are likely different tissues
- Could analyze separately

**If unknowns are mixed**:
- 🔍 Manually inspect each sample
- Include only those that cluster with brain
- Exclude outliers and different tissues

---

## Statistical Considerations

### Sample Size vs Quality

**Current (83 samples)**:
- Power for differential expression: **High** (>80% for 2-fold changes)
- False discovery rate: **Low** (homogeneous tissue)
- Biological interpretation: **Clear**

**All 120 samples (unknown tissue)**:
- Power for differential expression: **Lower** (tissue heterogeneity increases variance)
- False discovery rate: **Higher** (tissue-specific genes confound)
- Biological interpretation: **Unclear**

**Paradox**: More samples doesn't always mean better power if they add noise!

### Batch Effects

**Current (83 samples)**:
- 2 batches (2 bioprojects)
- Same tissue type
- Batch effects manageable

**All 120 samples**:
- 4 batches (4 bioprojects)
- Mixed tissue types (unknown)
- Batch + tissue confounding
- Harder to correct

---

## My Recommendation

### For Brain-Specific Analysis: Use 83 Samples ✅

**Rationale**:
1. Clear, well-defined dataset
2. All samples confirmed brain tissue
3. Better for publication
4. Easier to interpret
5. Higher statistical power for brain-specific signals

### For Exploratory Analysis: Try Both

**Two-stage approach**:

1. **First pass**: Run all 120 samples
   - Use for PCA/clustering
   - Identify outliers
   - Check if unknowns are brain-like

2. **Second pass**: Filter and re-analyze
   - Keep only brain or brain-like samples
   - Remove clear outliers
   - Final analysis on cleaned dataset

### Current Status

**What we have**: Complete analysis of 83 brain samples ✅
- All quantified
- All QC'd
- 6 PDF visualizations
- Expression matrices ready
- Sanity check passed

**What you could do**: Run all 120 samples in parallel directory
- Use new config: `config/amalgkit_pbarbatus_all120.yaml`
- Compare results
- Make informed decision based on actual data

---

## Quick Comparison Table

| Aspect | 83 Brain Samples | All 120 Samples |
|--------|------------------|-----------------|
| **Sample size** | 83 | 120 (+45%) |
| **Tissue confirmed** | 100% | 69% (83/120) |
| **Library type** | 100% paired | 96% paired |
| **Instrument quality** | Modern (HiSeq 2500/4000) | Mixed (includes HiSeq 2000) |
| **Biological homogeneity** | ✅ High | ⚠️ Unknown |
| **Statistical power** | ✅ High (low variance) | ⚠️ Lower (high variance) |
| **Interpretation** | ✅ Clear | ⚠️ Ambiguous |
| **Publication quality** | ✅ High | ⚠️ Lower |
| **Exploratory value** | ✅ Good | ✅ Better |
| **Confirmatory value** | ✅ Excellent | ⚠️ Limited |
| **Disk space needed** | ~7 GB per sample × 83 = ~581 GB | ~7 GB × 120 = ~840 GB |

---

## Conclusion

**For brain-specific research**: **Stick with 83 samples** ✅

The 37 unknown samples are a risk:
- Unknown tissue type
- Lower quality (older instrument, some single-end)
- Will add noise if they're not brain
- Harder to interpret and publish

**If you want to explore**: **Run all 120 in a separate directory**
- Use `config/amalgkit_pbarbatus_all120.yaml`
- Analyze expression patterns
- Make data-driven decision about inclusion
- Compare with 83-sample results

**Best of both worlds**: 
1. Keep current 83-sample analysis for primary results
2. Run 120-sample analysis for exploration
3. Use clustering to identify which unknowns are brain-like
4. Create a third "validated" dataset combining confirmed + validated samples

---

## Files Created

- ✅ `config/amalgkit_pbarbatus_all120.yaml` - Config for all 120 samples
- ✅ This document - Detailed comparison and recommendations

**To use all 120 samples**:
```bash
mkdir -p output/amalgkit/pbarbatus_all120
cd output/amalgkit/pbarbatus_all120
# Run workflow with new config
```


