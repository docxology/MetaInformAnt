# amalgkit cstmm: Cross-Species TMM Normalization

## Purpose

Performs cross-species transcript abundance normalization using TMM (Trimmed Mean of M-values) on orthologous genes identified via orthogroups. This step enables **comparative expression analysis** across multiple species.

## Overview

The `cstmm` step:
- Loads merged expression matrices from multiple species
- Maps transcripts to orthogroups (via OrthoFinder orthogroups)
- Identifies single-copy orthologs (1:1:1 orthologs)
- Applies TMM normalization using single-copy genes as reference
- Generates normalized expression matrices for cross-species comparison

## Usage

### Basic Usage

```bash
amalgkit cstmm \
  --out_dir output/amalgkit/work \
  --metadata output/amalgkit/work/metadata/metadata.tsv \
  --orthogroup_table output/orthogroups/Orthogroups.tsv \
  --dir_busco output/busco
```

### Python API

```python
from metainformant.rna import amalgkit

result = amalgkit.cstmm(
    out_dir="output/amalgkit/work",
    metadata="output/amalgkit/work/metadata/metadata.tsv",
    orthogroup_table="output/orthogroups/Orthogroups.tsv",
    dir_busco="output/busco"
)
```

### Configuration File

```yaml
steps:
  cstmm:
    out_dir: output/amalgkit/work
    orthogroup_table: output/orthogroups/Orthogroups.tsv
    dir_busco: output/busco
    dir_count: output/amalgkit/work/merge
```

## Parameters

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | PATH | `./` | Directory for intermediate and output files. |
| `--metadata` | PATH | `inferred` | Path to metadata.tsv. Default: `out_dir/metadata/metadata.tsv` |
| `--orthogroup_table` | PATH | `None` | Path to orthogroup table (OrthoFinder `Orthogroups.tsv` or `N0.tsv`). Set to empty string (`""`) for single-species TMM normalization. |
| `--dir_busco` | PATH | `None` | Path to directory containing per-species BUSCO full tables. File names expected: `GENUS_SPECIES_MISC.tsv` |
| `--dir_count` | PATH | `inferred` | Path to merged count matrices. Default: `out_dir/merge` |

## Input Requirements

### Prerequisites

- **Orthogroup Table**: From OrthoFinder analysis
- **Merged Expression Matrices**: From `amalgkit merge` for each species
- **BUSCO Results** (Optional): For quality filtering

### Orthogroup Table Format

**From OrthoFinder**: `Orthogroups.tsv`

```
Orthogroup              Species1_transcript1    Species1_transcript2    Species2_transcript1
OG0000001               Species1_001            Species1_002            Species2_001
OG0000002               Species1_003                                   Species2_002
OG0000002               Species1_004                                   Species2_003
```

**Required Format**:
- First column: Orthogroup ID (e.g., `OG0000001`)
- Remaining columns: Transcript IDs per species
- Tab-delimited

### Merged Count Matrices

Expected structure:
```
dir_count/
├── Species1/
│   └── Species1_tc.tsv    # Count matrix
├── Species2/
│   └── Species2_tc.tsv    # Count matrix
└── Species3/
    └── Species3_tc.tsv    # Count matrix
```

### BUSCO Full Tables (Optional)

Directory structure:
```
dir_busco/
├── Species1_full_table.tsv
├── Species2_full_table.tsv
└── Species3_full_table.tsv
```

**File Naming**: `{Genus}_{Species}_full_table.tsv`  
(Spaces replaced with underscores)

## Output Files

### Directory Structure

```
out_dir/cstmm/
└── {Species_List}/
    ├── {Species_List}_cstmm.tc.tsv          # ⭐ TMM-normalized counts matrix
    ├── {Species_List}_cstmm.sample_info.tsv # Sample metadata
    └── cstmm_summary.txt                    # Normalization statistics
```

### Primary Output

**`{Species_List}_cstmm.tc.tsv`**: TMM-normalized count matrix

**Format**:
```
orthogroup_id    Species1_SRR1    Species1_SRR2    Species2_SRR1    Species2_SRR2
OG0000001        1234.5           2456.0           1123.8           2234.2
OG0000002        567.2           1123.5           498.7            987.3
```

**Key Features**:
- Rows: Orthogroups (not individual transcripts)
- Columns: All samples from all species
- Values: TMM-normalized expression levels
- **Cross-species comparable**: Expression values normalized across species boundaries

## Workflow Integration

### Position in Pipeline

```mermaid
flowchart LR
    AmergeSpecies1[merge Species 1] --> B[cstmm]
    CmergeSpecies2[merge Species 2] --> B
    DmergeSpecies3[merge Species 3] --> B
    B --> E[curate]
```

**cstmm** runs **after merge** for all species, **before curate**.

### Single-Species Mode

```bash
# TMM normalization without orthogroups
amalgkit cstmm \
  --orthogroup_table "" \
  --out_dir output/work
```

**Effect**: Performs standard TMM normalization (within-species only)

## Cross-Species Normalization Theory

### Why Cross-Species Normalization?

**Problem**: Expression values from different species are not directly comparable:
- Different transcriptome sizes
- Different sequencing depths
- Different library preparation methods
- Different gene copy numbers

**Solution**: TMM normalization using single-copy orthologs:
1. Identify 1:1:1 orthologs (single-copy genes in all species)
2. Use these as reference genes (assumed constant expression)
3. Normalize all samples relative to reference
4. Enable cross-species expression comparison

### TMM Algorithm

**Trimmed Mean of M-values** (Robinson & Oshlack, 2010):

1. **Calculate log fold changes** (M-values) between samples
2. **Trim extreme values** (remove top/bottom 30% M-values)
3. **Calculate mean** of remaining M-values
4. **Apply scaling factor** to all genes

**Reference Genes**:
- Single-copy orthologs (1:1:1 across all species)
- High expression (above threshold)
- Low variance across samples

## Performance Considerations

### Runtime

- **Small analysis** (2-3 species, <50 samples): 1-5 minutes
- **Medium analysis** (3-5 species, 50-200 samples): 5-15 minutes
- **Large analysis** (>5 species, >200 samples): 15-60 minutes

**Factors**:
- Number of orthogroups
- Number of samples per species
- Computing resources

### Memory Usage

- **Moderate**: 2-8GB for typical analyses
- Scales with (# orthogroups) × (# total samples)

## Common Use Cases

### 1. Multi-Species Comparative Analysis

```bash
# After merging all species
amalgkit cstmm \
  --out_dir output/amalgkit/multi_species \
  --orthogroup_table output/orthofinder/Orthogroups.tsv \
  --dir_count output/amalgkit/multi_species/merge
```

**Result**: TMM-normalized matrix with all samples from all species

### 2. Cross-Species Brain Expression Comparison

```bash
# Compare brain expression across species
amalgkit cstmm \
  --out_dir output/comparative_brain \
  --orthogroup_table output/orthogroups/Orthogroups.tsv \
  --metadata output/comparative_brain/metadata/brain_samples.tsv
```

**Result**: Normalized expression matrix enabling cross-species brain comparisons

### 3. Single-Species TMM (No Orthogroups)

```bash
# Standard TMM normalization
amalgkit cstmm \
  --out_dir output/amalgkit/work \
  --orthogroup_table "" \
  --dir_count output/amalgkit/work/merge
```

**Result**: Within-species TMM normalization (like edgeR default)

## Preparing Orthogroup Tables

### Using OrthoFinder

```bash
# Run OrthoFinder on protein sequences
orthofinder \
  -f protein_sequences/ \
  -t 16 \
  -o output/orthofinder

# Extract orthogroups table
cp output/orthofinder/Orthogroups/Orthogroups.tsv output/orthogroups/
```

### Mapping Transcripts to Orthogroups

**Requirement**: Transcript IDs in expression matrices must match protein IDs in orthogroup table.

**Mapping Strategy**:
1. Use transcript-to-gene mapping (GFF files)
2. Map genes to orthogroups
3. Aggregate transcript counts to gene level
4. Map genes to orthogroups

**Or**: Ensure transcript IDs match protein IDs used in OrthoFinder

## Troubleshooting

### Issue: Orthogroup table not found

```
Error: Could not find orthogroup table
```

**Solutions**:
1. Verify file exists:
   ```bash
   ls -lh output/orthogroups/Orthogroups.tsv
   ```

2. Check file format:
   ```bash
   head -5 output/orthogroups/Orthogroups.tsv
   # Should show: Orthogroup_ID, Species1_transcript1, Species1_transcript2, ...
   ```

3. Specify full path:
   ```bash
   --orthogroup_table /absolute/path/to/Orthogroups.tsv
   ```

### Issue: No orthologs matched

```
Warning: No orthologs found in expression matrices
```

**Causes**:
1. Transcript ID mismatch between matrices and orthogroup table
2. Different ID formats (e.g., "XM_001" vs "Species1_XM_001")
3. Orthogroups based on proteins, but matrices have transcripts

**Solutions**:
1. Check ID formats:
   ```bash
   # In orthogroup table
   head -1 output/orthogroups/Orthogroups.tsv | cut -f2
   
   # In expression matrix
   head -1 output/work/merge/Species1/Species1_tc.tsv | cut -f1
   ```

2. Modify IDs to match:
   ```python
   import pandas as pd
   
   # Add species prefix if needed
   og = pd.read_csv("Orthogroups.tsv", sep="\t")
   og["Species1_transcript"] = "Species1_" + og["Species1_transcript"]
   ```

3. Use gene-level aggregation:
   - Map transcripts → genes → orthogroups
   - Aggregate transcript counts to genes before cstmm

### Issue: Too few single-copy orthologs

```
Warning: Only 50 single-copy orthologs found
```

**Causes**:
- Low-quality orthogroup prediction
- High gene duplication in study species
- Stringent filtering

**Solutions**:
1. Re-run OrthoFinder with relaxed parameters
2. Check BUSCO completeness scores (should be >90%)
3. Accept multi-copy orthogroups (modify cstmm logic)

### Issue: Empty cstmm output

**Diagnosis**:
```bash
wc -l output/work/cstmm/*/*_cstmm.tc.tsv
# Shows: 1 (header only, no data)
```

**Solutions**:
1. Verify merged matrices exist:
   ```bash
   ls output/work/merge/*/*_tc.tsv
   ```

2. Check metadata contains all species:
   ```bash
   cut -f2 output/work/metadata/metadata.tsv | sort -u
   ```

3. Verify orthogroup table has matching IDs

## Best Practices

### 1. Verify Orthogroup Quality

```bash
# Check BUSCO completeness before orthogroup analysis
for species in Species1 Species2 Species3; do
    busco_summary=$(grep "C:" ${species}_full_table.tsv | head -1)
    echo "$species: $busco_summary"
done

# Should show >90% complete BUSCOs for reliable orthogroups
```

### 2. Validate ID Mapping

```python
# Verify transcript IDs match orthogroup table
import pandas as pd

# Load orthogroups
og = pd.read_csv("Orthogroups.tsv", sep="\t")

# Load expression matrix
expr = pd.read_csv("Species1_tc.tsv", sep="\t", index_col=0)

# Check overlap
og_ids = set(og["Species1_transcript"].dropna())
expr_ids = set(expr.index)

overlap = og_ids & expr_ids
print(f"Overlap: {len(overlap)} / {len(og_ids)} = {len(overlap)/len(og_ids)*100:.1f}%")

# Should be >80% for good results
```

### 3. Inspect Normalization Factors

```bash
# Check cstmm summary
cat output/work/cstmm/*/cstmm_summary.txt

# Look for:
# - Number of single-copy orthologs used
# - Normalization factors (should be ~1.0)
# - Sample counts per species
```

### 4. Validate Cross-Species Comparability

```python
# After cstmm, verify expression values are comparable
import pandas as pd

cstmm = pd.read_csv("output/work/cstmm/multi_species/multi_species_cstmm.tc.tsv",
                    sep="\t", index_col=0)

# Check expression distributions across species
for species in ["Species1", "Species2", "Species3"]:
    species_cols = [c for c in cstmm.columns if species in c]
    print(f"{species}: mean={cstmm[species_cols].mean().mean():.2f}, "
          f"median={cstmm[species_cols].median().median():.2f}")

# Values should be similar across species (normalized)
```

## Real-World Examples

### Example 1: Four Ant Species Comparison

```bash
# Merge each species
for species in pbarbatus cfloridanus amellifera linepithema; do
    amalgkit merge \
      --out_dir output/amalgkit/${species}/work \
      --metadata output/amalgkit/${species}/work/metadata/pivot_qualified.tsv
done

# Cross-species TMM normalization
amalgkit cstmm \
  --out_dir output/amalgkit/comparative \
  --orthogroup_table output/orthofinder/Orthogroups.tsv \
  --dir_count output/amalgkit/comparative/merge
```

**Result**: TMM-normalized matrix with samples from all species in the merge directory

### Example 2: Insect Brain Multi-Species

```bash
amalgkit cstmm \
  --out_dir output/insect_brain_comparative \
  --orthogroup_table output/insect_orthogroups/Orthogroups.tsv \
  --dir_busco output/busco \
  --metadata output/insect_brain/metadata/brain_samples.tsv
```

**Result**: Normalized brain expression data across multiple insect species

## Integration with METAINFORMANT Workflow

### Multi-Species Workflow

```python
from metainformant.rna.workflow import execute_workflow, load_workflow_config

# Run workflow for each species
species_list = ["amellifera", "pbarbatus", "cfloridanus"]
for species in species_list:
    cfg = load_workflow_config(f"config/amalgkit_{species}.yaml")
    execute_workflow(cfg)  # Runs: metadata → select → getfastq → quant → merge

# Then run cross-species normalization
# (cstmm is separate, as it requires all species to be merged first)
```

### Configuration

```yaml
# config/amalgkit_comparative.yaml
steps:
  cstmm:
    orthogroup_table: output/orthofinder/Orthogroups.tsv
    dir_busco: output/busco
    dir_count: output/amalgkit/comparative/merge
```

## References

- **TMM Normalization**: Robinson & Oshlack (2010) "A scaling normalization method for differential expression analysis of RNA-seq data"
- **OrthoFinder**: https://github.com/davidemms/OrthoFinder
- **edgeR**: https://bioconductor.org/packages/release/bioc/html/edgeR.html
- **METAINFORMANT Workflow**: `docs/rna/workflow.md`

## See Also

- **Previous Step**: [`07_merge.md`](07_merge.md) - Merging expression matrices per species
- **Next Step**: [`09_curate.md`](09_curate.md) - Quality control on normalized data
- **Next Step**: [`10_csca.md`](10_csca.md) - Cross-species correlation analysis
- **Workflow Overview**: [`../amalgkit.md`](../amalgkit.md)

---

**Last Updated**: October 29, 2025  
**AMALGKIT Version**: 0.12.19  
**Status**: ✅ Production-ready, requires orthogroup tables from OrthoFinder


