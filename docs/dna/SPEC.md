# Specification: dna

## Scope
Comprehensive documentation for the DNA sequence analysis domain in MetaInformAnt, covering sequence processing, alignment, phylogenetics, population genetics, variant analysis, and NCBI genome retrieval.

## Architecture
- **Dependency Level**: Domain (relies on core.io, core.utils.config, math.population)
- **Component Type**: Analysis Module
- **Location**: `docs/dna/`
- **Source Root**: `src/metainformant/dna/`

## Data Structures
- **Documentation Files**:
  - `README.md`: Module overview and installation
  - `index.md`: Entry point with quick-start links
  - `AGENTS.md`: This file's AI attribution & directives
  - `SPEC.md`: Technical specification (this file)
  - `sequence.md`: FASTA/FASTQ processing, translation, reverse complement
  - `alignment.md`: Pairwise and multiple alignment algorithms
  - `phylogeny.md`: Tree building (UPGMA, NJ, ML), bootstrap
  - `population.md`: Population genetics (Tajima's D, Fst, heterozygosity, PCA)
  - `variants.md`: Variant calling, filtering, annotation (VEP, SnpEff)
  - `example_workflows.md`: End-to-end analysis pipelines

## Integration
- **Source**: `src/metainformant/dna/sequence/`, `src/metainformant/dna/alignment/`, `src/metainformant/dna/population/`, `src/metainformant/dna/phylogeny/`, `src/metainformant/dna/variants/`
- **Scripts**: `scripts/dna/` — CLI wrappers for common tasks
- **Tests**: `tests/dna/test_dna_sequence.py`, `tests/dna/test_dna_alignment.py`, `tests/dna/test_dna_population.py`, `tests/dna/test_dna_variants.py`
- **Dependencies**:
  - Core: `metainformant.core.io`, `metainformant.core.utils.config`
  - Math: `metainformant.math.coalescent`, `metainformant.math.fst`
  - Visualization: `metainformant.visualization.phylogeny`

## Key Functions (by submodule)

| Submodule | Key Functions | Purpose |
|-----------|---------------|---------|
| sequence  | `read_fasta()`, `reverse_complement()`, `translate()` | Basic sequence I/O and transformations |
| alignment | `global_align()`, `local_align()`, `multiple_align()`, `pairwise_distances()` | Sequence alignment algorithms |
| phylogeny | `upgma()`, `neighbor_joining()`, `maximum_likelihood()`, `bootstrap()` | Tree construction |
| population | `tajima_d()`, `weir_cockerham_fst()`, `heterozygosity()`, `pca()` | Population genetic statistics |
| variants  | `call_snps()`, `filter_vcf()`, `annotate_vep()` | Variant discovery and annotation |

## Testing Policy
- **Real Implementation**: All tests must use real implementations with actual genomic data. Synthetic data allowed only for edge cases.
- **Data**: Test datasets in `tests/data/dna/` (small FASTA, VCF, BAM fixtures).
- **Benchmarks**: Performance tests in `tests/benchmark/test_dna_speed.py` using `pytest-benchmark`.

## Public API
```python-snippet
# All submodules importable from top-level
from metainformant.dna import (
    alignment,      # pairwise & MSA
    phylogeny,      # tree building
    population,     # statistics
)
from metainformant.dna.sequence import composition, core as sequences
from metainformant.dna.variation import variants
```

## Examples
```python
from metainformant.dna import alignment, population
from metainformant.dna.sequence.core import gc_content, read_fasta

# 1. Read FASTA, compute GC%
seqs = read_fasta("genome.fasta")
for name, seq in seqs.items():
    print(f"{name}: GC={gc_content(seq):.1%}")

# 2. Align two sequences (global)
aligned1, aligned2 = alignment.global_align(seq1, seq2, match=2, mismatch=-1, gap=-2)

# 3. Compute Fst between two populations
fst = population.weir_cockerham_fst(pop1_vcf, pop2_vcf)
print(f"Mean Fst = {fst.mean():.3f}")
```

## Related Specifications
- **Core I/O**: [../core/SPEC.md](../core/SPEC.md) — File reading/writing infrastructure
- **Visualization**: [../visualization/SPEC.md](../visualization/SPEC.md) — Phylogenetic tree plotting
- **GWAS**: [../gwas/SPEC.md](../gwas/SPEC.md) — Variant association downstream
- **Task Reference**: [../../tasks/analyze_dna.md](../tasks/analyze_dna.md) — Quick reference guide
