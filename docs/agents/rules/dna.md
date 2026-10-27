# DNA Module Rules (`metainformant.dna`)

## Purpose
DNA sequence analysis, genomics, population genetics, and variant calling.

## Dependencies
- **Required**: `core`
- **Optional**: `biopython`, `ncbi-datasets-pylib`, `pysam`

## Source Structure
```
src/metainformant/dna/
 alignment/
 distances.py, msa.py, pairwise.py
 expression/
 codon.py, transcription.py, translation.py
 external/
 entrez.py, genomes.py, ncbi.py
 integration/
 rna.py
 io/
 fasta.py, fastq.py
 phylogeny/
 tree.py
 population/
 analysis.py, core.py, visualization.py
 sequence/
 composition.py, consensus.py, core.py, kmer.py, motifs.py, restriction.py
 variation/
 mutations.py, variants.py
```

## Package Management
- **ALWAYS use `uv`** for all Python package management and environment operations
- Install optional dependencies: `uv add biopython`, `uv add ncbi-datasets-pylib`, `uv add pysam`
- Use `uv run` to execute commands: `uv run pytest`, `uv run python scripts/dna/example.py`
- **NEVER use `pip` directly**

## Key Submodules

### Sequences (`sequences`)
- FASTA/FASTQ reading and writing
- Basic sequence operations
- Sequence validation

**Patterns**:
```python
from metainformant.dna.sequence import core as sequence_core

seqs = sequence_core.read_fasta("data/sequences.fasta")  # {sequence_id: sequence}
for seq_id, sequence in seqs.items():
    # Process sequence
    pass
```

### Alignment (`alignment`)
- Pairwise sequence alignment
- Global and local alignment algorithms
- Needleman-Wunsch and Smith-Waterman
- Alignment scoring

### Multiple Sequence Alignment (`msa`)
- Progressive alignment algorithm
- External tool integration (MUSCLE, MAFFT, ClustalW)
- Consensus generation from alignments
- Gap handling and optimization

### Phylogeny (`phylogeny`)
- Phylogenetic tree construction (Neighbor-Joining)
- Tree format: Newick
- Output to `output/dna/phylogeny/`

**Patterns**:
```python
from metainformant.dna.phylogeny import tree

tree_result = tree.neighbor_joining_tree(seqs)  # seqs: {sequence_id: aligned sequence}
# Returns a tree dictionary; tree.to_newick(tree_result) gives the Newick string
```

### Population (`population`)
- Population genetics statistics
- π (pi), Tajima's D, Fst
- `core` functions take aligned sequences and return plain numbers; `analysis.calculate_summary_statistics` returns a dictionary of named statistics

**Patterns**:
```python
from metainformant.dna.population import analysis, core as population

aligned_seqs = list(seqs.values())  # equal-length aligned sequences

pi = population.nucleotide_diversity(aligned_seqs)
# Returns a float, e.g. 0.001

stats = analysis.calculate_summary_statistics(sequences=aligned_seqs)
# Returns: {"nucleotide_diversity": 0.001, "segregating_sites": 42, "tajima_d": ..., ...}
```

### Population Analysis (`population_analysis`)
- Advanced population genetics analysis
- Extended population statistics
- Complex population genetic calculations

### Population Visualization (`population_viz`)
- Population genetics visualization
- Statistical plot generation
- Population structure visualization

### Variants (`variants`)
- Variant calling and analysis
- VCF file support
- Variant filtering and quality control

### Genomes (`genomes`)
- Genome data management
- Genome accession validation
- Genome metadata handling

### NCBI (`ncbi`, `genomes`)
- `ncbi`: E-utilities client (`NCBIClient`) for nucleotide/protein search, sequence fetch and taxonomy lookups; needs a contact email (`email=` or `NCBI_EMAIL`)
- `genomes`: accession validation and genome download through the NCBI Datasets API
- API integration with error handling

**Patterns**:
```python
from metainformant.dna.external import ncbi, genomes

# Validate accession (returns a bool)
if genomes.validate_accession("GCF_000001405.40"):
    # Download the genome package; returns the directory it was extracted to
    package_dir = genomes.download_genome_package(
        "GCF_000001405.40",
        output_dir="output/dna/genomes/",
    )

# Search nucleotide records through the E-utilities client
client = ncbi.NCBIClient(email="you@example.org")
records = client.search_nucleotide("Apis mellifera mitochondrion", max_results=5)
# Returns: list of dicts with id, accession, title, organism, length, moltype, created, updated
```

### Entrez (`entrez`)
- Entrez query construction and execution
- Sequence search and retrieval
- Database query utilities

### FastQ (`fastq`)
- FASTQ-specific processing
- Quality score analysis
- Read filtering

### Motifs (`motifs`)
- Pattern discovery
- Position Weight Matrix (PWM) analysis
- Motif enrichment

### Restriction (`restriction`)
- Restriction enzyme analysis
- Cutting site prediction
- Fragment size calculation

### Transcription (`transcription`)
- DNA to RNA transcription
- Transcription utilities
- RNA sequence generation

### Translation (`translation`)
- RNA to protein translation
- Genetic code handling
- Protein sequence generation

### Mutations (`mutations`)
- Mutation detection and analysis
- Mutation type classification
- Mutation spectrum analysis
- Mutation rate estimation

### Codon (`codon`)
- Codon usage analysis
- Codon frequency calculations
- Genetic code utilities

### Composition (`composition`)
- Sequence composition analysis
- GC content calculation
- Dinucleotide and trinucleotide frequencies
- Compositional bias detection

### Consensus (`consensus`)
- Consensus sequence generation
- Multiple sequence consensus
- Consensus quality scoring
- Ambiguity code handling

### Distances (`distances`)
- Evolutionary distance estimation
- Jukes-Cantor distance
- Kimura 2-parameter distance
- Tamura-Nei distance
- Distance matrix construction

### RNA Integration (`rna_integration`)
- Integration with RNA module
- Cross-domain analysis utilities

## Patterns

### I/O Operations
- **CRITICAL**: Always use `metainformant.core.io` for JSON/CSV/TSV operations
- Never use direct `import json` or `import csv`

**Pattern**:
```python
from metainformant.core import io

# Save results
io.dump_json(results, "output/dna/population/stats.json")
io.write_csv(df, "output/dna/variants/variants.csv")
```

### Lazy Imports
- Use lazy imports for optional dependencies
- Avoid import failures during unrelated commands

**Pattern**:
```python
def _lazy_import_biopython():
    try:
        from Bio import SeqIO
        return SeqIO
    except ImportError:
        return None
```

### File Path Handling
- Functions accept `Path` or `str` for file paths
- Always use `metainformant.core.io.paths` utilities for path operations
- Always resolve and validate paths using `paths.expand_and_resolve()` and `paths.is_within()`
- Default outputs to `output/dna/<analysis_type>/`

**Pattern**:
```python
from metainformant.core.io import paths
from pathlib import Path

# Expand and resolve user paths
path = paths.expand_and_resolve("~/data/file.txt")

# Security check before using user-provided paths
if paths.is_within(user_path, base_dir):
    # Safe to use
    pass
```

### Type Hints
- Comprehensive type hints throughout
- Use `from __future__ import annotations` for forward references
- Union types: `str | Path | None`

### Return Formats
- Population genetics stats: dictionaries with named keys
- Phylogenetic trees: Newick format strings
- Sequence data: iterators or lists of (id, sequence) tuples

### Error Handling
- Handle missing optional dependencies gracefully
- Provide clear error messages for missing data files
- Validate input sequences before processing

## Configuration

- No module-specific config (uses core config)
- Genome data paths can be specified via config
- Output paths default to `output/dna/`
- Environment prefix: `DNA_` (e.g., `DNA_THREADS`, `DNA_WORK_DIR`)

## Output Paths

- Phylogenetic trees: `output/dna/phylogeny/`
- Population genetics: `output/dna/population/`
- Variants: `output/dna/variants/`
- Genomes: `output/dna/genomes/`
- Alignments: `output/dna/alignments/`
- Multiple alignments: `output/dna/msa/`
- Mutations: `output/dna/mutations/`
- Composition: `output/dna/composition/`
- Consensus: `output/dna/consensus/`
- Distances: `output/dna/distances/`

## Integration

- **Used by**: `gwas` (variants, population genetics), `rna` (genomic coordinates), `information` (sequence analysis)
- **Integrates with**: `visualization` (phylogenetic trees), `math.population_genetics` (theoretical models)

## Reference

See Global Project Rules (root `.cursorrules`) for:
- Common directory structure and path handling
- Configuration patterns with env overrides
- Testing policy (real-implementation policy)
- Import patterns and code style
- Documentation guidelines

## Testing

- **STRICTLY REAL IMPLEMENTATION**: Test real implementations only (see Global Project Rules (root `.cursorrules`) real-implementation policy)
- Test with real FASTA/FASTQ files from `data/`
- Optional dependency tests should skip gracefully with `@pytest.mark.external_tool`
- Write test outputs to `output/dna/test/` using `tmp_path` fixture
- Real implementations only - use real file operations


