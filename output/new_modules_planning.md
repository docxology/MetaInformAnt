# New Modules Planning

This document outlines potential new modules for METAINFORMANT, their structure, integration points, and implementation priorities.

## Overview

The following modules are proposed to expand METAINFORMANT's capabilities in specific biological domains that are currently partially covered or missing:

1. **Metabolomics** - Comprehensive metabolite analysis and pathway mapping
2. **Metagenomics** - Taxonomic classification and functional annotation of microbial communities
3. **Chromatin/3D Structure** - Hi-C and chromatin conformation capture analysis
4. **Comparative Genomics** - Synteny analysis, gene family evolution, genome alignment
5. **Structural Biology** - Comprehensive protein and nucleic acid structure analysis

## 1. Metabolomics Module

### Rationale
Multiomics covers metabolomics partially, but a dedicated module would provide:
- Specialized metabolite identification and quantification
- Pathway mapping and enrichment analysis
- Integration with metabolic databases (KEGG, MetaboLights, HMDB)
- Statistical analysis for metabolomics data

### Proposed Structure
```
metabolomics/
├── __init__.py           # Module exports
├── README.md             # Module documentation
├── AGENTS.md             # AI contribution documentation
├── identification.py     # Metabolite identification algorithms
├── quantification.py     # Peak detection and quantification
├── pathways.py          # Pathway mapping and enrichment
├── integration.py       # Integration with multiomics and other modules
├── visualization.py     # Metabolomics-specific plots
└── workflows.py         # End-to-end metabolomics workflows
```

### Key Functions
- `identify_metabolites()` - MS/MS spectrum matching and identification
- `quantify_peaks()` - Peak detection and area quantification
- `map_to_pathways()` - KEGG pathway mapping
- `pathway_enrichment()` - Metabolic pathway enrichment analysis
- `integrate_with_rna()` - Correlation with transcriptomics data
- `integrate_with_protein()` - Integration with proteomics data

### Integration Points
- **Multiomics**: Add metabolomics layer to multi-omics integration
- **Networks**: Metabolic pathway network analysis
- **Information Theory**: Information content of metabolic profiles
- **Visualization**: Volcano plots, pathway diagrams, heatmaps

### Dependencies
- `core` (all utilities)
- `multiomics` (integration)
- `networks` (pathway networks)
- External: PyMS, mzML parsing libraries, KEGG API

### CLI Commands
```bash
metainformant metabolomics identify --input data/ms_data.mzML --output output/metabolomics
metainformant metabolomics pathway --metabolites data/metabolites.csv --output output/metabolomics
metainformant metabolomics integrate --metabolomics data/metabolomics.csv --transcriptomics data/rna.csv
```

## 2. Metagenomics Module

### Rationale
Distinct from ecology (which focuses on community diversity), metagenomics provides:
- Taxonomic classification from sequencing data
- Functional annotation of microbial communities
- Comparative metagenomics across samples
- Integration with single-cell and ecology modules

### Proposed Structure
```
metagenomics/
├── __init__.py           # Module exports
├── README.md             # Module documentation
├── AGENTS.md             # AI contribution documentation
├── taxonomy.py          # Taxonomic classification (Kraken, MetaPhlAn)
├── functional.py        # Functional annotation (HUMAnN, MG-RAST)
├── diversity.py         # Diversity analysis (distinct from ecology)
├── comparison.py        # Comparative metagenomics
├── integration.py        # Integration with single-cell and ecology
└── workflows.py         # End-to-end metagenomics workflows
```

### Key Functions
- `classify_taxonomy()` - Taxonomic classification from reads
- `annotate_function()` - Functional annotation of microbial genes
- `compare_samples()` - Comparative metagenomics analysis
- `integrate_singlecell()` - Integration with single-cell data
- `diversity_analysis()` - Alpha and beta diversity metrics

### Integration Points
- **Ecology**: Community structure analysis
- **Single-Cell**: Integration with single-cell microbial data
- **Quality**: Sequence quality for metagenomic reads
- **Information Theory**: Information content of microbial communities

### Dependencies
- `core` (all utilities)
- `ecology` (diversity metrics)
- `singlecell` (integration)
- External: Kraken, MetaPhlAn, HUMAnN (optional)

### CLI Commands
```bash
metainformant metagenomics classify --input data/reads.fastq --output output/metagenomics
metainformant metagenomics functional --input data/contigs.fasta --output output/metagenomics
metainformant metagenomics compare --samples data/samples/ --output output/metagenomics
```

## 3. Chromatin/3D Structure Module

### Rationale
Epigenome only covers methylation/tracks. This module would provide:
- Hi-C data processing and normalization
- Chromatin conformation capture (3C, 4C, 5C, Hi-C)
- 3D genome structure reconstruction
- Topologically associating domains (TADs) and loop calling
- Integration with epigenome and visualization

### Proposed Structure
```
chromatin/
├── __init__.py           # Module exports
├── README.md             # Module documentation
├── AGENTS.md             # AI contribution documentation
├── hic.py               # Hi-C data processing and normalization
├── domains.py           # TAD calling and domain analysis
├── loops.py             # Loop calling and interaction analysis
├── structure.py         # 3D structure reconstruction
├── integration.py        # Integration with epigenome
└── visualization.py     # Hi-C contact maps, 3D visualizations
```

### Key Functions
- `process_hic()` - Hi-C data processing and normalization
- `call_tads()` - Topologically associating domain calling
- `call_loops()` - Chromatin loop detection
- `reconstruct_3d()` - 3D genome structure reconstruction
- `integrate_methylation()` - Integration with DNA methylation data

### Integration Points
- **Epigenome**: Integration with methylation and chromatin marks
- **DNA**: Genomic coordinates and sequence analysis
- **Visualization**: Hi-C contact maps, 3D structure visualization
- **Networks**: Chromatin interaction networks

### Dependencies
- `core` (all utilities)
- `epigenome` (methylation integration)
- `dna` (genomic coordinates)
- External: Cooler, HiCExplorer (optional)

### CLI Commands
```bash
metainformant chromatin process-hic --input data/hic_data.hic --output output/chromatin
metainformant chromatin call-tads --hic data/normalized.hic --output output/chromatin
metainformant chromatin call-loops --hic data/normalized.hic --output output/chromatin
metainformant chromatin reconstruct-3d --hic data/normalized.hic --output output/chromatin
```

## 4. Comparative Genomics Module

### Rationale
Phylogeny exists but no comparative genomics workflows. This module would provide:
- Synteny analysis and block detection
- Gene family evolution and ortholog/paralog identification
- Whole-genome alignment
- Duplication and loss detection
- Integration with phylogeny and DNA modules

### Proposed Structure
```
comparative/
├── __init__.py           # Module exports
├── README.md             # Module documentation
├── AGENTS.md             # AI contribution documentation
├── synteny.py           # Synteny block detection and analysis
├── families.py          # Gene family evolution
├── orthologs.py         # Ortholog and paralog identification
├── alignment.py         # Whole-genome alignment
├── evolution.py         # Duplication and loss analysis
└── workflows.py         # End-to-end comparative genomics workflows
```

### Key Functions
- `detect_synteny()` - Synteny block detection
- `identify_gene_families()` - Gene family clustering
- `find_orthologs()` - Ortholog identification
- `align_genomes()` - Whole-genome alignment
- `analyze_duplications()` - Gene duplication and loss analysis

### Integration Points
- **DNA**: Sequence data and genome assemblies
- **Phylogeny**: Integration with phylogenetic trees
- **Visualization**: Synteny plots, gene family trees
- **Networks**: Gene family networks

### Dependencies
- `core` (all utilities)
- `dna` (sequences, alignment, phylogeny)
- External: LAST, MUMmer, OrthoFinder (optional)

### CLI Commands
```bash
metainformant comparative synteny --genome1 data/genome1.fasta --genome2 data/genome2.fasta --output output/comparative
metainformant comparative families --genomes data/genomes/ --output output/comparative
metainformant comparative orthologs --genomes data/genomes/ --output output/comparative
```

## 5. Structural Biology Module

### Rationale
Protein has basics but no full structure analysis. This module would provide:
- Structure prediction integration (AlphaFold, RoseTTAFold)
- Fold classification and domain architecture
- Structure comparison and alignment
- Structural phylogenetics
- Integration with protein and visualization modules

### Proposed Structure
```
structural/
├── __init__.py           # Module exports
├── README.md             # Module documentation
├── AGENTS.md             # AI contribution documentation
├── prediction.py        # Structure prediction integration
├── classification.py    # Fold classification (SCOP, CATH)
├── domains.py           # Domain architecture analysis
├── comparison.py        # Structure comparison and alignment
├── phylogenetics.py     # Structural phylogenetics
└── visualization.py     # 3D structure visualization
```

### Key Functions
- `predict_structure()` - Structure prediction integration
- `classify_fold()` - Fold classification (SCOP, CATH)
- `analyze_domains()` - Domain architecture analysis
- `compare_structures()` - Structure comparison and alignment
- `structural_phylogeny()` - Phylogenetic analysis from structures

### Integration Points
- **Protein**: Sequence-structure relationships
- **Phylogeny**: Structural phylogenetics
- **Visualization**: 3D structure visualization
- **Networks**: Structure-based protein networks

### Dependencies
- `core` (all utilities)
- `protein` (sequence analysis)
- `dna.phylogeny` (phylogenetic methods)
- External: BioPython, PyMOL, ChimeraX (optional)

### CLI Commands
```bash
metainformant structural predict --sequence data/protein.fasta --output output/structural
metainformant structural classify --structure data/structure.pdb --output output/structural
metainformant structural compare --structures data/structures/ --output output/structural
```

## Implementation Priority

### High Priority
1. **Metabolomics** - High demand for multi-omics integration
2. **Metagenomics** - Distinct from ecology, widely used

### Medium Priority
3. **Comparative Genomics** - Natural extension of phylogeny module
4. **Chromatin/3D Structure** - Growing field, complements epigenome

### Lower Priority
5. **Structural Biology** - Protein module provides basics, can be enhanced incrementally

## Integration Strategy

### Phase 1: Core Functionality
- Implement basic functions in each module
- Establish integration points with existing modules
- Create CLI commands

### Phase 2: Workflow Integration
- Create workflow functions
- Integrate with multiomics module
- Add visualization support

### Phase 3: Advanced Features
- Add advanced algorithms
- Integrate with external tools
- Comprehensive documentation

## Dependencies and Requirements

### External Tools
Many modules will benefit from optional external tools but should work without them:
- Metabolomics: PyMS, mzML parsers
- Metagenomics: Kraken, MetaPhlAn, HUMAnN
- Chromatin: Cooler, HiCExplorer
- Comparative: LAST, MUMmer, OrthoFinder
- Structural: BioPython, PyMOL, ChimeraX

### Core Module Usage
All new modules should:
- Use `core.io` for file I/O
- Use `core.config` for configuration
- Use `core.logging` for logging
- Use `core.paths` for path validation
- Write outputs to `output/` by default

## Testing Strategy

### Real Implementation Policy
All modules must follow the no-mocking policy:
- Use real data structures and algorithms
- Test with real (small) datasets
- Skip gracefully when external tools unavailable
- Document dependencies clearly

## Documentation Requirements

Each new module must include:
- Comprehensive README.md
- AGENTS.md documenting AI contributions
- API documentation in docstrings
- Usage examples
- Integration examples with other modules
- CLI documentation

## Next Steps

1. Review and prioritize modules based on user needs
2. Create module structure and basic scaffolding
3. Implement core functions following existing patterns
4. Add CLI integration
5. Create comprehensive documentation
6. Add integration tests with existing modules

