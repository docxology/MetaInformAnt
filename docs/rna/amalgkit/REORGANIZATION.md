# Documentation Reorganization Summary

**Date**: October 29, 2025  
**Purpose**: Separate general amalgkit workflow documentation from species-specific outputs

---

## Rationale

The original documentation was in `output/amalgkit/pbarbatus/`, which is a species-specific output directory. Since METAINFORMANT will be used to analyze many species, general workflow documentation should be in a central location that all species can reference.

---

## Files Moved to `docs/rna/amalgkit/`

### 1. `r_packages.md` (formerly `R_PACKAGE_SETUP.md`)
**Purpose**: R package installation and troubleshooting  
**Why moved**: R setup is identical across all species  
**Content**: Installation procedures, compilation issues, workarounds

### 2. `quick_start.md` (generalized from `QUICK_START.md`)
**Purpose**: Step-by-step guide for sanity and curate steps  
**Why moved**: Workflow is identical for all species, only paths differ  
**Changes**: 
- Replaced `pbarbatus` with `<species_name>` placeholder
- Added directory structure examples
- Generalized sample counts and file paths

### 3. `verify_template.sh` (generalized from `VERIFY_ALL_STEPS.sh`)
**Purpose**: Template verification script for QC  
**Why moved**: Can be copied and adapted for any species  
**Changes**:
- Auto-detects species name from directory
- Generic paths using wildcards
- Better error messages for wrong directory

---

## Files Kept in Species Directories

### `output/amalgkit/<species_name>/`

Each species directory now contains:

1. **`README.md`** - Species-specific summary with links to general docs
2. **`QUICK_REFERENCE.md`** - Species-specific data (sample counts, gene counts, statistics)
3. **`verify_workflow.sh`** - Copy of template adapted for the species
4. **`work/`** - All amalgkit working files
5. **`analysis/`** - Pre-built expression matrices and visualizations

---

## Benefits of This Organization

### For Users
- ✅ Clear separation of general workflow vs. species-specific data
- ✅ No duplicate documentation across species directories
- ✅ Easy to find workflow guides in central documentation
- ✅ Species directories focus on results, not instructions

### For Maintenance
- ✅ Update workflow docs once, applies to all species
- ✅ Species directories are smaller and cleaner
- ✅ Easier to add new species without recreating docs
- ✅ Clear documentation hierarchy

### For Scalability
- ✅ Support many species without documentation bloat
- ✅ Consistent workflow across all analyses
- ✅ Easy to compare species-specific results
- ✅ Template-based approach for new species

---

## Directory Structure

```
docs/rna/amalgkit/
├── README.md                    # Overview and index
├── quick_start.md               # General workflow guide ⭐
├── r_packages.md                # R setup guide ⭐
├── verify_template.sh           # Verification script template ⭐
├── END_TO_END_WORKFLOW.md       # Complete workflow
├── comprehensive_guide.md       # Detailed guide
└── [other general docs...]

output/amalgkit/
├── pbarbatus/                   # P. barbatus example
│   ├── README.md                # Species summary + links to docs
│   ├── QUICK_REFERENCE.md       # Species-specific data
│   ├── verify_workflow.sh       # Copied from template
│   ├── analysis/                # Expression matrices
│   └── work/                    # Amalgkit outputs
│
├── dmelanogaster/               # Next species (future)
│   ├── README.md
│   ├── QUICK_REFERENCE.md
│   ├── verify_workflow.sh
│   └── work/
│
└── <your_species>/              # Your analysis
    ├── README.md                # Copy template, customize
    ├── QUICK_REFERENCE.md       # Add your data
    ├── verify_workflow.sh       # Copy from docs
    └── work/
```

---

## Usage for New Species

### Step 1: Setup Directory
```bash
mkdir -p output/amalgkit/my_species
cd output/amalgkit/my_species
```

### Step 2: Copy Templates
```bash
# Copy verification script
cp ~/docs/rna/amalgkit/verify_template.sh ./verify_workflow.sh
chmod +x verify_workflow.sh

# Copy README template (then customize)
cp ../pbarbatus/README.md ./README.md
```

### Step 3: Reference Workflow Docs
- Read `docs/rna/amalgkit/quick_start.md` for workflow steps
- Read `docs/rna/amalgkit/r_packages.md` for R setup
- Run `./verify_workflow.sh` to validate outputs

### Step 4: Create Quick Reference
After analysis completes, create `QUICK_REFERENCE.md` with:
- Sample count and sources
- Gene/transcript counts
- Expression statistics
- Key file locations
- Species-specific findings

---

## Migration Checklist

- [x] Move R package setup to `docs/rna/amalgkit/r_packages.md`
- [x] Generalize and move quick start to `docs/rna/amalgkit/quick_start.md`
- [x] Generalize and move verification script to `docs/rna/amalgkit/verify_template.sh`
- [x] Create concise species-specific README for pbarbatus
- [x] Keep species-specific QUICK_REFERENCE.md in pbarbatus
- [x] Update amalgkit README to reference new organization
- [x] Copy verification script to pbarbatus for continuity
- [x] Delete old files from pbarbatus output directory

---

## Example: P. barbatus

The P. barbatus directory now has:
- **Concise README** linking to general docs
- **QUICK_REFERENCE** with species-specific data summary
- **verify_workflow.sh** copied from template
- **All results** in work/ and analysis/

This serves as the reference example for new species analyses.

---

**Result**: Clean, scalable documentation structure ready for multi-species analyses! 🎉

