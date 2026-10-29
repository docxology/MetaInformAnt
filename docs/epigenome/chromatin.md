# Chromatin State Learning

Simplified ChromHMM-style chromatin state discovery and annotation, implemented in `metainformant.epigenome.chromatin_state.state_learning`, using Expectation-Maximization on a Gaussian Mixture Model. Segments the genome into functionally distinct chromatin states from histone modification data.

## Key Concepts

**Chromatin states** are recurring combinations of histone modifications that correspond to functional genomic elements (active promoters, enhancers, heterochromatin, etc.). States are discovered in an unsupervised manner from multi-mark ChIP-seq or similar data.

**EM/GMM learning** fits a Gaussian Mixture Model to the epigenomic signal matrix, where each component represents a chromatin state characterized by its mean signal and variance for each mark (a diagonal covariance: marks are treated as independent). A state-to-state transition matrix is estimated from the responsibilities of consecutive bins and is used for decoding.

**Viterbi decoding** assigns each genomic bin to a chromatin state using the learned emission parameters, transition matrix and state weights, producing the most likely state sequence along the genome.

**Biological interpretation** maps learned states to named chromatin categories (for example "Active Promoter", "Active Enhancer", "Polycomb Repressed", "Heterochromatin") based on which histone marks are high in each state's emission means.

## Function Reference

### `learn_chromatin_states(signal_matrix, n_states=15, max_iter=200, tol=1e-4, random_state=42) -> Dict`

Discover chromatin states from a signal matrix (bins x marks) using EM-based GMM fitting with diagonal covariance. Centroids are initialized with k-means++ using a private generator. The same input and integer `random_state` reproduce the model; `random_state=None` draws a fresh seed without reading or changing Python's global random state.

**Parameters:**
- `signal_matrix`: 2D list or numpy array of signal values (genomic bins x histone marks). `ValueError` if it is an empty list, a numpy array that is not 2-D, has rows of different lengths, or has fewer bins than `n_states`.
- `n_states`: Number of states to learn (default 15).
- `max_iter`: Maximum EM iterations (default 200).
- `tol`: Positive, finite convergence threshold for the log-likelihood change (default 1e-4).
- `random_state`: Integer initialization seed (default 42), or `None` for fresh initialization. Non-finite signals, invalid dimensions/state counts and non-positive iteration limits raise `ValueError`.

**Returns** dict with:
- `states`: Number of states.
- `emission_params`: dict with `means` (states x marks) and `variances` (states x marks, the diagonal of each state's covariance).
- `transition_matrix`: State-to-state transition probabilities (states x states).
- `assignments`: Most likely state of each bin (list of int): the state with the highest responsibility, without the transition-based decoding of `assign_states`.
- `weights`: Mixing weight (prior probability) of each state.
- `log_likelihood`: Final log-likelihood.
- `n_iterations`: Number of EM iterations performed.
- `converged`: Whether the log-likelihood change fell below `tol`.
- `state_counts`: Number of bins assigned to each state, as a dict keyed by state index (states with no bins are absent).

### `assign_states(signal_matrix, model) -> List[int]`

Assign each genomic bin to a state by Viterbi decoding with the emission parameters, transition matrix and weights of `model` (the dict returned by `learn_chromatin_states`), so neighboring bins influence each other. `ValueError` if the matrix has a different number of marks than the model.

### `interpret_states(emission_params, mark_names) -> List[Dict]`

Label each state from its emission means. `emission_params` is `model["emission_params"]` (not the whole model) and `mark_names` is required, in the column order of the signal matrix. Marks are recognized by substring, case-insensitively (`K4me3`, `K4me1`, `K27ac`, `K27me3`, `K9me3`, `K36me3`), and a mark counts as high in a state when its mean is strictly above the mean across states plus half a standard deviation. Labels are assigned in this order: "Bivalent Promoter" (K4me3 and K27me3), "Active Promoter" (K4me3 and K27ac), "Promoter" (K4me3), "Active Enhancer" (K4me1 and K27ac), "Poised Enhancer" (K4me1 and K27me3), "Weak Enhancer" (K4me1), "Polycomb Repressed" (K27me3), "Heterochromatin" (K9me3), "Transcribed" (K36me3), "Active Regulatory" (K27ac), "Quiescent" (no mark high), otherwise "Unknown".

Returns per-state dicts with `state_id`, `label`, `category`, `dominant_marks` (the marks that are high in the state), `mean_signal`, `description` and `emission_means`. Non-finite emission values or a mark-name count that differs from the number of emission columns raise `ValueError`.

### `compute_state_enrichment(state_assignments, annotations, n_states) -> Dict`

Test enrichment of chromatin states in annotations. `annotations` maps an annotation name (promoters, gene bodies, ...) to the indices of the bins it covers, and `n_states` is the total number of states (not the number of bins). For every state and annotation it computes the observed and expected overlap, the fold enrichment and a one-sided Fisher's exact test p-value (an exact hypergeometric tail when SciPy is unavailable).

Returns a dict with `enrichment_matrix` (keyed `"state_<id>__<annotation>"`; each value has `state_id`, `annotation`, `observed_overlap`, `expected_overlap`, `fold_enrichment`, `p_value`, `n_state_bins` and `n_annotation_bins`), `summary` (per state: `state_id`, `n_bins`, `most_enriched_annotation`, `best_fold_enrichment`, `best_p_value`), `n_bins`, `n_states` and `n_annotations`. Without annotations, the best-annotation and associated enrichment/p-value summary fields are `None`.

### `segment_genome(state_assignments, bin_size=200, min_segment=1) -> List[Dict]`

Convert per-bin state assignments into contiguous segments. Adjacent bins with the same state are merged, and segments shorter than `min_segment` bins are dropped. `ValueError` for an empty list.

Returns segments with `state`, `start_bin`, `end_bin` (exclusive), `start_bp`, `end_bp`, `length_bins` and `length_bp`. Bins are numbered from 0 and there is no chromosome field.

### `compare_chromatin_states(states_a, states_b) -> Dict`

Compare two models returned by `learn_chromatin_states` (for example two conditions or cell types binned the same way). The comparison is bin by bin on their `assignments`, which must have equal lengths, and it assumes state numbers are comparable between the two models: emission profiles are not correlated or matched.

Returns `n_bins`, `concordance`, `discordance`, `concordant_bins`, `discordant_bins`, `state_switches` (counts keyed `"<a>-><b>"`), `switch_matrix`, `differential_segments`, `n_differential_segments` and `summary`.

## Usage Examples

```python
import numpy as np

from metainformant.epigenome.chromatin_state.state_learning import (
    assign_states,
    compute_state_enrichment,
    interpret_states,
    learn_chromatin_states,
    segment_genome,
)

# Signal matrix: genomic bins x 5 histone marks, built from three block types
mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3", "H3K27me3"]
block_types = np.array(
    [
        [2.5, 0.1, 3.0, 0.2, 0.0],  # Active promoter-like
        [0.1, 2.0, 2.5, 0.1, 0.0],  # Enhancer-like
        [0.0, 0.0, 0.0, 0.0, 3.0],  # Repressed-like
    ]
)
rng = np.random.default_rng(42)
block_order = np.repeat([0, 0, 1, 2, 2, 1, 0, 2], 25)  # 25-bin blocks
signal = block_types[block_order] + rng.normal(0, 0.15, (len(block_order), 5))

# Learn 3 chromatin states
model = learn_chromatin_states(signal, n_states=3, max_iter=200)
print(f"Converged: {model['converged']}, iterations: {model['n_iterations']}")

# Assign states to genomic bins (Viterbi decoding)
assignments = assign_states(signal, model)

# Interpret state biology from the emission parameters
labels = interpret_states(model["emission_params"], mark_names)
for label in labels:
    print(f"State {label['state_id']}: {label['label']} ({label['description']})")

# Segment genome into contiguous state blocks
segments = segment_genome(assignments, bin_size=200)

# Enrichment of states in an annotation given as bin indices
annotations = {"promoter": [0, 5, 10, 50]}
enrichment = compute_state_enrichment(assignments, annotations, n_states=model["states"])
```

## Configuration

No dedicated configuration file or environment prefix is used by the chromatin state module; call the functions above directly from Python.

The algorithms are pure Python. numpy arrays are accepted as input and converted to nested lists, and scipy is used for Fisher's exact test when it is installed.

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- histone mark ChIP-seq peaks
- `metainformant.epigenome.assays.atacseq` -- chromatin accessibility data
- `metainformant.epigenome.peak_calling` -- signal-based peak detection
- `metainformant.epigenome.workflow` -- end-to-end pipelines
