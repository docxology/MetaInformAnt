# Specification: dna

## Scope
Comprehensive documentation for the DNA sequence analysis domain in MetaInformAnt, covering sequence processing, alignment, phylogenetics, population genetics, variant analysis, and NCBI genome retrieval.

## Architecture
- **Dependency Level**: Domain (relies on `core.utils` for logging, errors and Newick helpers, plus `core.io`, `core.sequence` and `core.ncbi`; imports nothing from `math` or `visualization`)
- **Component Type**: Analysis Module
- **Location**: `docs/dna/`
- **Source Root**: `src/metainformant/dna/`

## Data Structures
- **Documentation Files**:
  - `README.md`: Module overview and installation
  - `index.md`: Entry point with quick-start links
  - `AGENTS.md`: AI attribution & directives for this folder
  - `SPEC.md`: Technical specification (this file)
  - `sequences.md`: FASTA reading, reverse complement, GC content, k-mers
  - `alignment.md`: Pairwise alignment (`msa.md` covers multiple alignment)
  - `phylogeny.md`: Tree building (UPGMA, NJ), Newick export, bootstrap support
  - `population.md`: Population genetics (π, Tajima's D and other neutrality tests, Hudson's Fst, heterozygosity)
  - `variants.md`: VCF parsing (`calling.md` covers variant calling and filtering)
  - `population_workflows.md`: Population genetics workflows (summary statistics, population comparison, neutrality tests)

## Integration
- **Source**: `src/metainformant/dna/sequence/`, `src/metainformant/dna/expression/`, `src/metainformant/dna/alignment/`, `src/metainformant/dna/population/`, `src/metainformant/dna/phylogeny/`, `src/metainformant/dna/variation/`, `src/metainformant/dna/external/`
- **Scripts**: `scripts/dna/run_dna_analysis.py` — command-line wrapper for DNA analysis
- **Tests**: `tests/dna/test_dna_sequences.py`, `tests/dna/test_dna_alignment.py`, `tests/dna/test_dna_population_genetics.py`, `tests/dna/test_dna_variants_vcf.py`
- **Dependencies**:
  - Core: `metainformant.core.utils` (logging, errors, Newick helpers), `metainformant.core.io`, `metainformant.core.sequence`, `metainformant.core.ncbi`
  - Math: none imported by `dna`. `metainformant.math.population_genetics.coalescent` and `metainformant.math.population_genetics.fst` hold the theory-side counterparts (coalescent expectations, Weir and Cockerham's Fst from haplotype counts)
  - Visualization: population genetics plots live in `metainformant.dna.population.visualization`; phylogenetic tree plots come from `metainformant.visualization.genomics.trees` (`plot_phylo_tree`)

## Key Functions (by submodule)

| Submodule | Key Functions | Purpose |
|-----------|---------------|---------|
| sequence  | `read_fasta()`, `write_fasta()`, `reverse_complement()`, `gc_content()` (in `sequence.core`) | Basic sequence I/O and transformations |
| expression | `transcribe_dna_to_rna()` (in `dna.transcription`), `translate_dna()` (in `expression.translation`) | Transcription and translation |
| alignment | `pairwise.global_align()`, `pairwise.local_align()`, `msa.align_msa()`, `distances.distance_matrix()` | Sequence alignment algorithms and distances |
| phylogeny | `upgma_tree()`, `neighbor_joining_tree()`, `bootstrap_support()`, `to_newick()` | Distance-based tree construction, bootstrap support and Newick export |
| population | `tajimas_d()`, `hudson_fst()`, `observed_heterozygosity()`, `nucleotide_diversity()` | Population genetic statistics |
| variation | `parse_vcf()`, `call_variants_pileup()`, `filter_variants_by_quality()`, `annotate_variants()` | Variant parsing, calling, filtering and annotation from a supplied annotation table |

## Testing Policy
- **Real Implementation**: All tests must use real implementations with actual genomic data. Synthetic data allowed only for edge cases.
- **Data**: The one stored fixture is `tests/dna/data/dna/toy.fasta`; other test inputs are built inline or written under `tmp_path`.
- **Benchmarks**: There are no DNA-specific benchmarks. Repository benchmark scripts live in `tests/perf/`, and `pytest-benchmark` is a development dependency.

## Public API
```python-snippet
# All subpackages are importable from the top level
from metainformant.dna import (
    alignment,      # pairwise, msa and distances modules
    phylogeny,      # tree building
    population,     # statistics
)
from metainformant.dna.sequence import composition, core as sequences
from metainformant.dna.variation import variants
```

## Examples
```python
from pathlib import Path

from metainformant.dna.alignment import pairwise
from metainformant.dna.population import hudson_fst
from metainformant.dna.sequence.core import gc_content, read_fasta

# 1. Read FASTA, compute GC% (a tiny file keeps the example self-contained)
Path("genome.fasta").write_text(">seq1\nATGCGCGTAT\n>seq2\nATGCGTAT\n")
seqs = read_fasta("genome.fasta")
for name, seq in seqs.items():
    print(f"{name}: GC={gc_content(seq):.1%}")

# 2. Align two sequences (global)
result = pairwise.global_align(seqs["seq1"], seqs["seq2"], match=2, mismatch=-1, gap=-2)
aligned1, aligned2 = result.aligned_seq1, result.aligned_seq2

# 3. Compute Fst between two populations (lists of aligned sequences)
pop1 = ["ATGCGC", "ATGCGC", "ATGCGT", "ATGCGC"]
pop2 = ["ATACGC", "ATACGT", "ATACGT", "ATACGC"]
fst = hudson_fst(pop1, pop2)
print(f"Fst = {fst:.3f}")
```

## Related Specifications
- **Core I/O**: [../core/SPEC.md](../core/SPEC.md) — File reading/writing infrastructure
- **Visualization**: [../visualization/SPEC.md](../visualization/SPEC.md) — Phylogenetic tree plotting
- **GWAS**: [../gwas/SPEC.md](../gwas/SPEC.md) — Variant association downstream
- **Task Reference**: [../tasks/analyze_dna.md](../tasks/analyze_dna.md) — Quick reference guide
