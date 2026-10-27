# Agent Directives: metagenomics

## Role
Documentation agent for the metagenomics module covering amplicon and shotgun analysis.

## Module Scope
- 16S/ITS amplicon profiling (OTU clustering, ASV denoising)
- Taxonomic classification and tree building
- Shotgun metagenome assembly and binning, with bin quality computed from marker gene annotations
- Community profiling and relative abundance
- Functional annotation (ORF prediction, gene families)
- Metabolic pathway reconstruction and completeness
- Alpha/beta diversity, PERMANOVA, ordination and differential abundance
- Microbiome visualization (Krona, rarefaction, ordination)

## Key Source Files
- `src/metainformant/metagenomics/amplicon/` - OTU/ASV processing and taxonomy
- `src/metainformant/metagenomics/shotgun/` - Assembly, binning, profiling
- `src/metainformant/metagenomics/functional/` - Annotation and pathways
- `src/metainformant/metagenomics/diversity/` - Diversity metrics, rarefaction, PERMANOVA, ordination
- `src/metainformant/metagenomics/comparative/` - Differential abundance and biomarkers
- `src/metainformant/metagenomics/visualization/` - Microbiome plots

## External Tools (not wrapped)
The package wraps none of these tools; run them yourself and load their results into the plain dicts and lists the functions take:
- QIIME2 for amplicon preprocessing
- MetaSPAdes/MEGAHIT for production assembly
- MetaBAT2/CONCOCT for production binning
- HMMER/CheckM for the single-copy marker gene annotations (`marker_genes`) that `assess_bin_quality` and `refine_bins` need
- DIAMOND/HMMER for functional annotation at scale (`annotate_genes` only scores against the position-specific matrices you supply)
