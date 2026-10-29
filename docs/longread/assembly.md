# Long-Read Assembly

The assembly module provides tools for de novo genome assembly from long reads, including minimizer-based overlap finding, consensus sequence generation, iterative polishing, and hybrid assembly with short reads.

## Key Concepts

### Minimizer Sketching

Minimizers are a compact representation of k-mer content that enables efficient seed-based overlap detection. The algorithm selects the lexicographically smallest k-mer in each window of w consecutive k-mers, creating a sparse but representative sketch of each read.

### Overlap-Layout-Consensus (OLC)

The classical OLC assembly approach:
1. **Overlap**: Find all pairwise read overlaps using minimizer index
2. **Layout**: Build and traverse the overlap graph
3. **Consensus**: Generate a consensus sequence from overlapping reads

### Consensus Generation

The implementation aligns reads independently to a backbone and counts bases at each backbone position:

1. Select the supplied backbone or the longest read
2. Align reads by semi-global dynamic programming (a full matrix up to 10^8 cells, then a k-mer anchored approximation)
3. Build a position-weight matrix
4. Call the most frequent base at each backbone position
5. Optionally polish with iterative re-alignment

Insertions relative to the backbone and extensions past its ends are not added to the consensus, whose length is capped by that backbone.

### Hybrid Assembly

Combines long reads (for scaffolding and spanning repeats) with short reads (for error correction). Long reads provide structural continuity while short reads improve per-base accuracy.

## Data Structures

### Minimizer

```python
@dataclass
class Minimizer:
    hash_value: int     # Integer hash of the k-mer
    position: int       # Position in sequence (0-based)
    is_reverse: bool    # From reverse complement
    kmer: str           # Actual k-mer string
```

### Overlap

```python
@dataclass
class Overlap:
    query_name: str; query_length: int
    query_start: int; query_end: int
    target_name: str; target_length: int
    target_start: int; target_end: int
    strand: str         # "+" or "-"
    n_matches: int; alignment_score: float
    identity: float
```

### ConsensusResult

```python
@dataclass
class ConsensusResult:
    sequence: str           # Consensus sequence
    quality: list[float]    # Per-base Phred quality scores
    coverage: list[int]     # Per-base read coverage depths
    num_reads: int          # Reads used in consensus
    length: int
    mean_quality: float
    mean_coverage: float
```

## Function Reference

### Overlap Finding

```python-snippet
def minimizer_sketch(
    sequence: str, k: int = 15, w: int = 10,
) -> list[Minimizer]
```

Compute minimizer sketch for a sequence. `k` is k-mer size, `w` is window size.

```python-snippet
def find_overlaps(
    reads: dict[str, str],
    k: int = 15, w: int = 10, min_overlap: int = 500,
) -> list[Overlap]
```

Find all pairwise overlaps using minimizer index. `reads` maps read_name to sequence.

```python-snippet
def compute_overlap_graph(overlaps: list[Overlap]) -> dict[str, Any]
def filter_contained_reads(overlaps: list[Overlap]) -> list[Overlap]
```

Build and clean the overlap graph. `filter_contained_reads` removes reads fully contained within larger reads.

### Consensus

```python-snippet
def generate_consensus(reads: Sequence[str | dict[str, Any]], backbone: str | None=None, match_score: int=2, mismatch_penalty: int=-4, gap_open: int=-4, gap_extend: int=-2) -> ConsensusResult
def polish_consensus(consensus: str | ConsensusResult, reads: Sequence[str | dict[str, Any]], iterations: int=2) -> ConsensusResult
def multiple_sequence_alignment(sequences: Sequence[str], gap_open: int=-4, gap_extend: int=-2) -> MSAResult
def calculate_consensus_quality(consensus: str | ConsensusResult, reads: Sequence[str | dict[str, Any]]) -> dict[str, Any]
```

`generate_consensus` aligns reads to a backbone and calls a consensus; a single read has coverage 1 and quality 0.0 because it provides no independent error evidence. `polish_consensus` iteratively realigns the reads. `multiple_sequence_alignment` aligns each read to the backbone independently: its rows are not guaranteed to be column-consistent when insertions differ, and a warning identifies this limitation. Use a full multiple-sequence aligner when shared insertion columns are required. `calculate_consensus_quality` reports quality metrics.

### Hybrid Assembly

```python-snippet
def hybrid_assemble(long_reads: Sequence[str | dict[str, Any]], short_reads: Sequence[str | dict[str, Any]], min_long_read_length: int=2000, kmer_size: int=21, min_overlap: int=1000) -> HybridAssemblyResult
def correct_with_short_reads(long_reads: Sequence[str | dict[str, Any]], short_reads: Sequence[str | dict[str, Any]], kmer_size: int=21, min_kmer_support: int=3, window_size: int=50) -> list[CorrectedRead]
def scaffold_with_long_reads(contigs: Sequence[str | dict[str, Any]], long_reads: Sequence[str | dict[str, Any]], min_alignment_length: int=1000, min_bridge_reads: int=2) -> list[Scaffold]
```

`hybrid_assemble` combines long and short reads. `correct_with_short_reads` polishes a single long read using short read alignments. `scaffold_with_long_reads` orders contigs using long-read spanning information; it does not infer contig orientation.

## Usage Examples

```python-snippet
from metainformant.longread import assembly

# Compute minimizer sketches
sketch = assembly.minimizer_sketch("ACGTACGT" * 100, k=15, w=10)

# Find overlaps between reads
reads = {"read1": seq1, "read2": seq2, "read3": seq3, ...}
overlaps = assembly.find_overlaps(reads, k=15, w=10, min_overlap=500)
print(f"Found {len(overlaps)} overlaps")

# Build overlap graph
graph = assembly.compute_overlap_graph(overlaps)
clean_overlaps = assembly.filter_contained_reads(overlaps)

# Generate consensus
read_sequences = [reads[name] for name in selected_reads]
result = assembly.generate_consensus(read_sequences)
print(f"Consensus: {result.length} bp, mean Q{result.mean_quality:.1f}")

# Polish consensus
polished = assembly.polish_consensus(result.sequence, read_sequences, iterations=3)

# Hybrid assembly
hybrid_result = assembly.hybrid_assemble(long_reads, short_reads)

# Scaffold contigs with long reads
scaffolds = assembly.scaffold_with_long_reads(contigs, long_reads)
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: numpy
- Default k-mer size (k=15) balances sensitivity and specificity for typical error rates
- Window size (w=10) controls sketch density
- Minimum overlap default (500bp) filters spurious short overlaps

## Related Modules

- `longread.io` -- Read FAST5/BAM files for assembly input
- `longread.quality` -- Quality filtering before assembly
- `longread.utils` -- Batch processing and assembly summaries
- `longread.visualization` -- `plot_dotplot` for visualizing overlaps
