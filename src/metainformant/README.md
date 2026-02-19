# METAINFORMANT

## Overview

METAINFORMANT: Comprehensive Bioinformatics Toolkit for Multi-Omic Analysis.

## 📦 Contents

- **[core/](core/)**
- **[dna/](dna/)**
- **[ecology/](ecology/)**
- **[epigenome/](epigenome/)**
- **[gwas/](gwas/)**
- **[information/](information/)**
- **[life_events/](life_events/)**
- **[longread/](longread/)**
- **[math/](math/)**
- **[menu/](menu/)**
- **[metabolomics/](metabolomics/)**
- **[metagenomics/](metagenomics/)**
- **[ml/](ml/)**
- **[multiomics/](multiomics/)**
- **[networks/](networks/)**
- **[ontology/](ontology/)**
- **[pharmacogenomics/](pharmacogenomics/)**
- **[phenotype/](phenotype/)**
- **[protein/](protein/)**
- **[quality/](quality/)**
- **[rna/](rna/)**
- **[simulation/](simulation/)**
- **[singlecell/](singlecell/)**
- **[spatial/](spatial/)**
- **[structural_variants/](structural_variants/)**
- **[visualization/](visualization/)**
- `[__init__.py](__init__.py)`
- `[__main__.py](__main__.py)`

## 📊 Structure

```mermaid
graph TD
    metainformant[metainformant]
    style metainformant fill:#f9f,stroke:#333,stroke-width:2px

    subgraph "Core & Utils"
        core[core]
        menu[menu]
        visualization[visualization]
        quality[quality]
    end

    subgraph "Omics"
        dna[dna]
        rna[rna]
        protein[protein]
        epigenome[epigenome]
        metagenomics[metagenomics]
        metabolomics[metabolomics]
        pharmacogenomics[pharmacogenomics]
        structural_variants[structural_variants]
        longread[longread]
        singlecell[singlecell]
        spatial[spatial]
        multiomics[multiomics]
    end

    subgraph "Analysis"
        gwas[gwas]
        ml[ml]
        networks[networks]
        simulation[simulation]
        math[math]
        information[information]
        phenotype[phenotype]
        ecology[ecology]
        life_events[life_events]
        ontology[ontology]
    end

    metainformant --> core
    metainformant --> dna
    metainformant --> rna
    metainformant --> protein
```

## Usage

Import modules directly:

```python
from metainformant import dna, rna, protein
from metainformant.gwas.analysis import association
from metainformant.longread.workflow.orchestrator import LongReadOrchestrator
```
