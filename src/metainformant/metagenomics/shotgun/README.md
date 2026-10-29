# Shotgun

Shotgun metagenomics analysis providing de Bruijn graph assembly, tetranucleotide-based contig binning into MAGs, and k-mer-based community profiling with LCA classification.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports assembly, binning, profiling submodules |
| `assembly.py` | De Bruijn graph assembler with multi-k, scaffolding, and N50 stats |
| `binning.py` | Composition/coverage-based binning and CheckM-style quality assessment from supplied marker gene annotations |
| `profiling.py` | K-mer taxonomic profiling (Kraken-style) with LCA classification |

## Key Functions

| Function | Description |
|----------|-------------|
| `assembly.assemble_contigs()` | Assemble reads into contigs using de Bruijn graph approach |
| `assembly.scaffold_contigs()` | Scaffold contigs from paired-end reads (FR libraries): order, orientation (`Scaffold.orientations`) and gap sizes are inferred from the strands and positions of the mates |
| `assembly.calculate_assembly_stats()` | Compute N50, L50, GC content, and other assembly metrics |
| `binning.bin_contigs()` | Bin contigs into genome bins using TNF and coverage features |
| `binning.refine_bins()` | Refine bins using single-copy marker gene assessment (requires `marker_genes`); keeps only bins and split sub-bins that meet both the completeness and the contamination threshold |
| `binning.assess_bin_quality()` | Compute completeness, contamination and quality score from per-contig marker gene annotations (`None` without them) |
| `binning.count_high_quality_bins()` | Count assessed bins with completeness >= 0.9 and contamination <= 0.05 |
| `profiling.build_kmer_index()` | Build a canonical (strand-independent) k-mer index from reference genomes |
| `profiling.profile_community()` | Classify reads and generate community profiles |
| `profiling.calculate_relative_abundance()` | Calculate taxon relative abundances |

## Usage

```python
from metainformant.metagenomics.shotgun import assembly, binning, profiling

contigs = assembly.assemble_contigs(reads, k_range=[21, 33])
stats = assembly.calculate_assembly_stats(contigs)
bins = binning.bin_contigs(contig_sequences, coverage_data)  # bin quality is None without marker_genes
# marker_genes: contig ID -> single-copy marker names found on it (from a marker search such as HMMER)
assessed = binning.assess_bin_quality(bins.bins, marker_genes=marker_genes)
index = profiling.build_kmer_index(refs, taxonomy, k=31)
profile = profiling.profile_community(reads, database=index)
```
