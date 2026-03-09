# Tissue Patching System

The Tissue Patching system in MetaInformAnt provides a robust mechanism for correcting and automating tissue metadata assignment across thousands of samples. This ensures that downstream analyses like `csca` (Cross-Species Correlation Analysis) use accurate, canonical tissue labels even when NCBI metadata is missing or ambiguous.

## Overview

NCBI SRA/ENA metadata often contains:
1.  **Missing Values**: The `tissue` field is empty.
2.  **Ambiguous Terms**: "Nervous system", "head", "central nervous system".
3.  **Experimental Details**: "Brain of forager kept in a group for 2 days".

The MetaInformAnt system uses a Two-Tier normalization strategy:
-   **Synonym Mapping**: Mapping varied strings to canonical names (e.g., "Nervous system" → `brain`).
-   **Patching**: Force-assigning tissues to specific Runs, BioProjects, or BioSamples based on manual research or study titles.

## Configuration Files

The system is controlled by two YAML files in `config/amalgkit/`:

### 1. `tissue_mapping.yaml`
Maps canonical tissue names to a list of synonyms.
```yaml
brain:
  - nerve
  - nervous system
  - head
  - brain
```

### 2. `tissue_patches.yaml`
High-priority overrides for specific accessions.
```yaml
samples:
  SRR12345678: brain
bioprojects:
  PRJNA339620: mushroom_body
biosamples:
  SAMN00849801: brain
```

## Normalization Logic

The `StreamingPipelineOrchestrator` and `scripts/rna/normalize_tissue_metadata.py` apply normalization in the following priority order:

1.  **Sample Patch**: If the Run Accession (SRR) exists in `samples:`.
2.  **BioSample Patch**: If the BioSample Accession exists in `biosamples:`.
3.  **BioProject Patch**: If the BioProject Accession exists in `bioprojects:`.
4.  **Synonym Match**: If the raw metadata value matches a synonym in `tissue_mapping.yaml`.
5.  **Prefix Match**: If the raw metadata starts with a known synonym (e.g., "Brain tissue..." matches `brain`).

## Tools

### Metadata Normalization Script
**Path**: `scripts/rna/normalize_tissue_metadata.py`

Used to batch-process a `metadata.tsv` file and add a `tissue_normalized` column.

```bash
python3 scripts/rna/normalize_tissue_metadata.py \
    --input output/amalgkit/pbarbatus/work/metadata/metadata.tsv \
    --output output/amalgkit/pbarbatus/work/metadata/metadata_normalized.tsv
```

### Coverage Verification
**Path**: `scripts/rna/verify_tissue_coverage.py`

Reports how many samples are mapped vs. unmapped in a species workflow.

```bash
python3 scripts/rna/verify_tissue_coverage.py --species pbarbatus
```

## Validation

The tissue patching system is validated via:
-   **Automated Tests**: [test_rna_tissue_normalization.py](../../../tests/test_rna_tissue_normalization.py) ensures patching priority and synonym matching work correctly.
-   **Production Audit**: Current coverage for the honeybee dataset is **99.9%** (7,265/7,270 samples) using these patches.
