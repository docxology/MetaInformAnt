# Agent Directives: metagenomics

## Role
Documentation agent for the metagenomics module covering amplicon and shotgun analysis.

## Module Scope
- 16S/ITS amplicon profiling (OTU clustering, ASV denoising)
- Taxonomic classification and tree building
- Shotgun metagenome assembly and binning
- Community profiling and relative abundance
- Functional annotation (ORF prediction, gene families)
- Metabolic pathway reconstruction and completeness
- Microbiome visualization (Krona, rarefaction, ordination)

## Key Source Files
- `src/metainformant/metagenomics/amplicon/` - OTU/ASV processing and taxonomy
- `src/metainformant/metagenomics/shotgun/` - Assembly, binning, profiling
- `src/metainformant/metagenomics/functional/` - Annotation and pathways
- `src/metainformant/metagenomics/visualization/` - Microbiome plots

## External Dependencies
- QIIME2 for amplicon workflows
- MetaSPAdes/MEGAHIT for assembly
- MetaBAT2/CONCOCT for binning
- DIAMOND/HMMER for functional annotation
