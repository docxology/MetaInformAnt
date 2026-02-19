# Agent Directives: metainformant

## 🤖 Role

Specialized agent context for the `metainformant` component.

## 🛠️ Tools & Capabilities

- **Context**: METAINFORMANT: Comprehensive Bioinformatics Toolkit for Multi-Omic Analysis (27 modules).
- **Core & Utils**: `core` (I/O, logging, config), `menu` (interactive CLI), `quality` (QC, batch effects), `visualization` (plotting)
- **Molecular Omics**: `dna`, `rna`, `protein`, `epigenome`, `longread`, `structural_variants`
- **Higher-order Omics**: `singlecell`, `spatial`, `multiomics`, `metabolomics`, `metagenomics`, `pharmacogenomics`
- **Analysis & Methods**: `gwas`, `ml` (ML + deep learning + LLM), `networks`, `simulation`, `math`, `information`
- **Annotation & Ecology**: `ontology`, `phenotype`, `ecology`, `life_events`
- **Pattern**: Source Code Pattern

## ⚠️ Rules & Constraints

- **Imports**: Prefer absolute imports from `metainformant`.
- **I/O**: Use `metainformant.core.io` for all file operations.
- **Logging**: Use `metainformant.core.utils.logging`.
