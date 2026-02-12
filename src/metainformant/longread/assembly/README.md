# Assembly

Long-read assembly module providing minimizer-based overlap computation, partial-order alignment consensus generation, and hybrid assembly combining long and short reads.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports overlap, consensus, hybrid submodules |
| `overlap.py` | Minimizer sketching and all-vs-all overlap detection |
| `consensus.py` | POA-style consensus generation with iterative polishing |
| `hybrid.py` | Hybrid assembly combining long and short reads |

## Key Functions

| Function | Description |
|----------|-------------|
| `overlap.minimizer_sketch()` | Compute minimizer sketches for read sequences |
| `overlap.find_overlaps()` | Find all-vs-all overlaps using minimizer index |
| `overlap.compute_overlap_graph()` | Build overlap graph from detected overlaps |
| `consensus.generate_consensus()` | Generate consensus sequence from read pile-up |
| `consensus.polish_consensus()` | Iteratively polish consensus with re-alignment |
| `consensus.multiple_sequence_alignment()` | Align multiple reads for consensus building |
| `hybrid.hybrid_assemble()` | End-to-end hybrid assembly pipeline |
| `hybrid.correct_with_short_reads()` | Error-correct long reads using short-read k-mers |
| `hybrid.scaffold_with_long_reads()` | Scaffold short-read contigs with long reads |

## Usage

```python
from metainformant.longread.assembly import overlap, consensus, hybrid

sketch = overlap.minimizer_sketch(sequence, k=15, w=10)
overlaps = overlap.find_overlaps(reads)
result = consensus.generate_consensus(reads)
```
