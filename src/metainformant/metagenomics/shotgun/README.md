# Shotgun

Shotgun metagenomics analysis providing de Bruijn graph assembly, tetranucleotide-based contig binning into MAGs, and k-mer-based community profiling with LCA classification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports assembly, binning, profiling submodules |
| `assembly.py` | De Bruijn graph assembler with multi-k, scaffolding, and N50 stats |
| `binning.py` | Composition/coverage-based binning and CheckM-style quality assessment |
| `profiling.py` | K-mer taxonomic profiling (Kraken-style) with LCA classification |

## Key Functions

| Function | Description |
|----------|-------------|
| `assembly.assemble_contigs()` | Assemble reads into contigs using de Bruijn graph approach |
| `assembly.scaffold_contigs()` | Scaffold contigs using paired-end information |
| `assembly.calculate_assembly_stats()` | Compute N50, L50, GC content, and other assembly metrics |
| `binning.bin_contigs()` | Bin contigs into genome bins using TNF and coverage features |
| `binning.refine_bins()` | Refine bins using single-copy marker gene assessment |
| `binning.assess_bin_quality()` | Assess bin completeness and contamination |
| `profiling.build_kmer_index()` | Build k-mer index from reference genomes |
| `profiling.profile_community()` | Classify reads and generate community profiles |
| `profiling.calculate_relative_abundance()` | Calculate taxon relative abundances |

## Usage

```python
from metainformant.metagenomics.shotgun import assembly, binning, profiling

contigs = assembly.assemble_contigs(reads, k=31)
stats = assembly.calculate_assembly_stats(contigs)
bins = binning.bin_contigs(contigs, coverage_data)
profile = profiling.profile_community(reads, kmer_index)
```
