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
| `compute_network_motifs()` | Count motifs of a directed edge list (a feed-forward loop is A -> B, B -> C and the direct A -> C) with z-scores against 100 randomised networks; takes `random_state` |
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
motifs = grn_inference.compute_network_motifs(grn["edges"], motif_size=3, random_state=0)
pwm = motif_analysis.build_pwm(aligned_sequences)
hits = motif_analysis.scan_sequence_for_motifs(sequence, motif_library)
```

`compute_network_motifs` draws its randomised networks from a generator of its own: `random_state` is an integer seed or a `random.Random`, the same seed gives the same z-scores, and the global `random` state is neither seeded nor consumed (`None` draws fresh entropy).

## Behaviour notes

- `tf_list` of `infer_grn_correlation`, `infer_grn_mutual_info` and `infer_grn_regression` restricts the regulators. Names
  that are not genes of the matrix are ignored and repeated names count once, but a list that matches no gene raises
  `ValueError`; it never widens silently to every gene.
- `score_regulators()` reports `enrichment_p`, the exact one-sided Fisher p-value (hypergeometric upper tail, summed in
  log space, with or without SciPy). Its population is the set of distinct targets in the GRN and its successes are the
  gene-set members among them; gene-set entries the GRN never targets, and duplicates, are not counted. The p-value is
  not rounded.
- `compute_network_motifs()` counts motifs of two or three distinct nodes: autoregulation (self-loops) and repeated edges
  are ignored, and a z-score is `+inf` / `-inf` when the randomised counts never vary and the observed count differs.
- `validate_grn()` evaluates its curves once per distinct score, so `auroc` (ties count half, as in the Mann-Whitney
  statistic) and `auprc` (step-wise average precision) do not depend on the order of tied edges or on `PYTHONHASHSEED`.
- `find_tf_binding_motifs()` ranks k-mers by the exact binomial upper tail of their count under the expected frequency.

## Related inference functions

`infer_grn()` and `infer_regulatory_network_from_expression()` live in `metainformant.networks.interaction.regulatory_analysis`, next to `GeneRegulatoryNetwork` (see [`../interaction/`](../interaction/)):

- `infer_grn()` supports only the methods `'correlation'` and `'mutual_info'`; any other value raises `ValueError`. For Granger causality on time series use `metainformant.information.network_info.information_flow.information_flow_network(series, method="granger_causality")`.
- `infer_regulatory_network_from_expression()` supports the same two methods. `'grnboost2'` is not implemented (it needs the `arboreto` package) and raises `NotImplementedError`.
- Both take `random_state` for the `'mutual_info'` method (scikit-learn's mutual-information estimator adds a random jitter); one generator is shared across all gene pairs of a run.
- A feed-forward loop in `regulatory_motifs()` and `GeneRegulatoryNetwork.regulatory_motifs()` needs all three edges, including the direct TF1 -> target edge; a chain without it is a cascade.
