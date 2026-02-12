# Chromatin State Learning

ChromHMM-style chromatin state discovery and annotation using Expectation-Maximization with Gaussian Mixture Models. Segments the genome into functionally distinct chromatin states from histone modification data.

## Key Concepts

**Chromatin states** are recurring combinations of histone modifications that correspond to functional genomic elements (active promoters, enhancers, heterochromatin, etc.). States are discovered in an unsupervised manner from multi-mark ChIP-seq or similar data.

**EM/GMM learning** fits a Gaussian Mixture Model to the epigenomic signal matrix, where each component represents a chromatin state characterized by its mean signal levels across marks and a covariance structure.

**Viterbi decoding** assigns each genomic bin to its most likely chromatin state using the learned model parameters, producing a genome-wide segmentation.

**Biological interpretation** maps learned states to known chromatin annotations (e.g., "Active TSS", "Strong Enhancer", "Heterochromatin") based on the mark emission patterns of each state.

## Function Reference

### `learn_chromatin_states(signal_matrix, n_states=10, max_iter=200, tol=1e-6, seed=42) -> Dict`

Discover chromatin states from a signal matrix (bins x marks) using EM-based GMM fitting.

**Parameters:**
- `signal_matrix`: 2D list of signal values (genomic bins x histone marks).
- `n_states`: Number of states to learn (default 10).
- `max_iter`: Maximum EM iterations (default 200).
- `tol`: Convergence threshold for log-likelihood change.

**Returns** dict with:
- `means`: Per-state mean signal levels (states x marks).
- `covariances`: Per-state covariance matrices.
- `weights`: Prior probability of each state.
- `log_likelihood`: Final log-likelihood.
- `n_iterations`: Number of EM iterations completed.
- `converged`: Whether the model converged.

### `assign_states(signal_matrix, model) -> List[int]`

Assign each genomic bin to its most likely state via Viterbi decoding (MAP assignment from the learned GMM posterior).

### `interpret_states(model, mark_names=None) -> List[Dict]`

Generate biological interpretations of each learned state based on emission patterns. Labels states as "Active TSS", "Enhancer", "Transcription", "Heterochromatin", etc.

Returns per-state dicts with `state_id`, `label`, `dominant_marks`, and `description`.

### `compute_state_enrichment(state_assignments, genomic_features, n_bins) -> Dict`

Test enrichment of chromatin states in genomic features (promoters, gene bodies, etc.) using Fisher's exact test or hypergeometric approximation.

### `segment_genome(state_assignments, bin_size=200) -> List[Dict]`

Convert per-bin state assignments into contiguous genomic segments. Adjacent bins with the same state are merged.

Returns segments with `chrom`, `start`, `end`, `state`, and `length`.

### `compare_chromatin_states(model_a, model_b, mark_names=None) -> Dict`

Compare chromatin state models between two conditions or cell types. Computes state correspondence by correlating emission profiles.

## Usage Examples

```python
from metainformant.epigenome import (
    learn_chromatin_states,
    assign_states,
    interpret_states,
    compute_state_enrichment,
    segment_genome,
)

# Signal matrix: 1000 genomic bins x 5 histone marks
# Marks: H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K27me3
signal = [[2.5, 0.1, 3.0, 0.2, 0.0],   # Active promoter-like
          [0.1, 2.0, 2.5, 0.1, 0.0],   # Enhancer-like
          [0.0, 0.0, 0.0, 0.0, 3.0],   # Repressed-like
          # ... more bins
         ]

# Learn 8 chromatin states
model = learn_chromatin_states(signal, n_states=8, max_iter=200)
print(f"Converged: {model['converged']}, iterations: {model['n_iterations']}")

# Assign states to genomic bins
assignments = assign_states(signal, model)

# Interpret state biology
mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3", "H3K27me3"]
labels = interpret_states(model, mark_names=mark_names)
for label in labels:
    print(f"State {label['state_id']}: {label['label']} ({label['description']})")

# Segment genome into contiguous state blocks
segments = segment_genome(assignments, bin_size=200)

# Enrichment in promoters
features = {"promoter": [0, 5, 10, 50]}  # bin indices
enrichment = compute_state_enrichment(assignments, features, n_bins=len(signal))
```

## Configuration

Environment variable prefix: `EPI_`

All algorithms are pure Python with numpy for matrix operations.

## Related Modules

- `metainformant.epigenome.chipseq` -- input histone mark peaks
- `metainformant.epigenome.atacseq` -- chromatin accessibility data
- `metainformant.epigenome.peak_calling` -- signal-based peak detection
- `metainformant.epigenome.workflow` -- end-to-end pipelines
