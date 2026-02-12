# Chromatin State Learning

ChromHMM-style chromatin state discovery from histone modification data using Gaussian mixture models and EM algorithm, with state assignment, biological interpretation, and cross-condition comparison.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `state_learning` module |
| `state_learning.py` | State discovery, assignment, interpretation, enrichment, segmentation |

## Key Functions

| Function | Description |
|----------|-------------|
| `learn_chromatin_states()` | Discover chromatin states from histone modification signal matrix via GMM/EM |
| `assign_states()` | Assign learned states to genomic regions using Viterbi-like decoding |
| `interpret_states()` | Biologically interpret states based on known histone mark patterns |
| `compute_state_enrichment()` | Test state enrichment against genomic annotations |
| `segment_genome()` | Segment genome into contiguous state blocks |
| `compare_chromatin_states()` | Compare chromatin states across conditions or cell types |
| `compute_state_transition_rates()` | Calculate state-to-state transition rate matrix |

## Usage

```python
from metainformant.epigenome.chromatin_state import state_learning

states = state_learning.learn_chromatin_states(signal_matrix, n_states=10)
assignments = state_learning.assign_states(signal_matrix, states)
interpretation = state_learning.interpret_states(states)
enrichment = state_learning.compute_state_enrichment(assignments, annotations)
```
