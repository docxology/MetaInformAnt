# Pharmacogenomics

## Overview

Clinical pharmacogenomic analysis module for METAINFORMANT. Covers star allele calling, metabolizer phenotyping, CPIC guideline lookups, drug interaction prediction, and report generation.

## Contents

- **alleles/** - Star allele calling, diplotype determination, phenotype classification
- **annotations/** - CPIC guideline lookups, PharmGKB queries, FDA drug label parsing
- **metabolism/** - CPIC activity score computation, metabolizer phenotype prediction
- **interaction/** - Drug-drug interaction prediction, polypharmacy risk, CYP profiling
- **clinical/** - ACMG variant classification, drug-gene analysis, clinical reports
- **visualization/** - Metabolizer distributions, allele frequencies, drug response heatmaps

## Architecture

```mermaid
graph TD
    subgraph "Pharmacogenomics Module"
        AL[alleles/] --> |star_allele.py| SA[Star Allele Calling]
        AL --> |diplotype.py| DI[Diplotype Determination]
        AL --> |phenotype.py| PH[Phenotype Classification]

        AN[annotations/] --> |cpic.py| CP[CPIC Guidelines]
        AN --> |pharmgkb.py| PG[PharmGKB Queries]
        AN --> |drug_labels.py| DL[FDA Drug Labels]

        ME[metabolism/] --> |metabolizer_status.py| MS[Activity Score + Phenotype]

        IN[interaction/] --> |drug_interactions.py| DDI[Drug-Drug Interactions]

        CL[clinical/] --> |reporting.py| RP[Clinical Reports]
        CL --> |pathogenicity.py| PA[ACMG Classification]
        CL --> |drug_interaction.py| CI[Drug-Gene Interactions]

        VZ[visualization/] --> |plots.py| PL[Publication Plots]
    end

    SA --> MS
    CP --> RP
    CI --> RP
```

## Usage

```python
from metainformant.pharmacogenomics.alleles import star_allele, diplotype, phenotype
from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
from metainformant.pharmacogenomics.metabolism import metabolizer_status
from metainformant.pharmacogenomics.interaction import drug_interactions
from metainformant.pharmacogenomics.clinical import reporting, pathogenicity
```
