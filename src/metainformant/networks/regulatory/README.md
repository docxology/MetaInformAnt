# Regulatory Network Analysis

Gene regulatory network inference from expression data and transcription factor binding motif analysis, including correlation-based, mutual information-based, and regression-based approaches.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `grn_inference` and `motif_analysis` |
| `grn_inference.py` | GRN inference (correlation, MI/ARACNE, regression), motif detection, validation |
| `motif_analysis.py` | PWM construction, motif scoring, TF binding site scanning |

## Key Functions

| Function | Description |
|----------|-------------|
| `infer_grn_correlation()` | Infer GRN using Pearson/Spearman correlation |
| `infer_grn_mutual_info()` | ARACNE-like MI-based network inference |
| `infer_grn_regression()` | Regression-based GRN inference |
| `score_regulators()` | Score transcription factor regulatory activity |
| `compute_network_motifs()` | Detect network motifs (feed-forward loops, etc.) |
| `validate_grn()` | Validate inferred network against known interactions |
| `build_pwm()` | Build position weight matrix from aligned sequences |
| `score_motif_match()` | Score a sequence against a PWM |
| `find_tf_binding_motifs()` | Discover TF binding motifs via k-mer overrepresentation |
| `scan_sequence_for_motifs()` | Scan a sequence for matches to known motif library |

## Usage

```python
from metainformant.networks.regulatory import grn_inference, motif_analysis

grn = grn_inference.infer_grn_correlation(expression_matrix, gene_names)
grn_mi = grn_inference.infer_grn_mutual_info(expression_matrix, gene_names)
motifs = grn_inference.compute_network_motifs(grn, motif_size=3)
pwm = motif_analysis.build_pwm(aligned_sequences)
hits = motif_analysis.scan_sequence_for_motifs(sequence, motif_library)
```
