# Metagenomics Module Rules

## Purpose
Metagenomic analysis including amplicon profiling (16S/ITS), shotgun metagenomics (assembly, binning, profiling), functional annotation, and microbial ecology visualization.

## Source Structure
```
src/metainformant/metagenomics/
├── amplicon/
│   ├── otu_clustering.py  # OTU clustering, identity calculation, chimera filtering
│   ├── asv_denoising.py   # ASV denoising, error rate estimation, paired-read merging
│   └── taxonomy.py        # Taxonomic classification, tree building, confidence scoring
├── shotgun/
│   ├── assembly.py        # Metagenome assembly, scaffolding, assembly statistics
│   ├── binning.py         # Contig binning, tetranucleotide frequency, bin QC
│   └── profiling.py       # Community profiling, relative abundance, k-mer indexing
├── functional/
│   ├── annotation.py      # Gene annotation, ORF prediction, gene family classification
│   └── pathways.py        # Pathway reconstruction, completeness, profile comparison
└── visualization/
    └── plots.py           # Krona charts, stacked bars, rarefaction, ordination, heatmaps
```

## Dependencies
- **Required**: numpy, pandas, scipy
- **Optional**: QIIME2, MetaSPAdes/MEGAHIT (assembly), MetaBAT2/CONCOCT (binning), DIAMOND/HMMER (annotation)

## Import Patterns
```python
from metainformant.metagenomics.amplicon import otu_clustering, asv_denoising, taxonomy
from metainformant.metagenomics.shotgun import assembly, binning, profiling
from metainformant.metagenomics.functional import annotation, pathways
```

## Configuration
- Environment prefix: `META_` (e.g., `META_THREADS`, `META_WORK_DIR`)
- Output path: `output/metagenomics/<analysis_type>/` (e.g., `amplicon/`, `assembly/`, `binning/`, `functional/`)

## Integration
- **Metagenomics → Ecology**: Community diversity from amplicon/shotgun data
- **Metagenomics → Networks**: Microbial co-occurrence networks
- **Metagenomics → Ontology**: Functional annotation via GO/KEGG

## Testing
- Use `@pytest.mark.external_tool` for tests requiring QIIME2, MetaSPAdes
- Generate test amplicon/shotgun data programmatically
- All test outputs to `tmp_path`
