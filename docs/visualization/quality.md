# Quality Control Plots

Quality control visualization functions for quality control metrics including single-cell QC panels, sequencing QC summaries, per-base quality plots, adapter content plots, and read length distributions. The implementations live in `metainformant.visualization.analysis.quality_omics` and `metainformant.visualization.analysis.quality_sequencing`; `metainformant.visualization.analysis.quality` re-exports them, together with `quality_assessment`, under one namespace. Inputs are type-checked, so pass the `dict`, `list` or `numpy.ndarray` described below. Every function accepts `output_path=` to save the figure at 300 dpi.

## Functions

### `plot_singlecell_qc_metrics(qc_metrics, *, ax=None, output_path=None, figsize=(14, 10), **kwargs)`

Plot multiple single-cell quality control metrics. `qc_metrics` is a dict of per-cell arrays; the
recognized keys are `n_counts`, `n_genes`, `percent_mito` and `doublet_score`. The function always draws its own
2 x 3 panel figure (`ax` is not used) and returns the first axes of that figure.

**Example:**
```python
from metainformant.visualization.analysis.quality_omics import plot_singlecell_qc_metrics
import numpy as np

metrics = {
    'n_counts': np.random.poisson(1000, 100),
    'n_genes': np.random.poisson(2000, 100),
    'percent_mito': np.random.uniform(0, 10, 100)
}
ax = plot_singlecell_qc_metrics(metrics)
```

### `plot_quality_metrics(qc_data, *, ax=None, output_path=None, **kwargs)`

Plot a 2 x 2 summary of sequencing QC metrics and return the first axes (`ax` is not used; `figsize` is the only keyword argument read). `qc_data` is a non-empty dict; each panel is drawn only when its key is present:

- `per_base_quality`: `{'positions': [...], 'mean_qualities': [...]}`
- `gc_content_distribution`: `{'bins': [...], 'counts': [...]}`, with one more bin edge than counts
- `sequence_length_distribution`: `{'lengths': [...], 'counts': [...]}`
- `basic_statistics`: a dict with `num_reads`, `total_bases`, `min_length`, `max_length` and `mean_length`

**Example:**
```python
from metainformant.visualization.analysis.quality_sequencing import plot_quality_metrics
import numpy as np

positions = np.arange(1, 101)
qc_data = {
    'per_base_quality': {'positions': positions, 'mean_qualities': 38 - 0.05 * positions},
    'gc_content_distribution': {
        'bins': np.arange(0, 101, 10),
        'counts': [2, 5, 20, 60, 140, 180, 90, 30, 8, 2],
    },
    'sequence_length_distribution': {'lengths': [148, 149, 150], 'counts': [10, 40, 950]},
    'basic_statistics': {
        'num_reads': 1000,
        'total_bases': 149_900,
        'min_length': 148,
        'max_length': 150,
        'mean_length': 149.9,
    },
}
ax = plot_quality_metrics(qc_data)
```

### `plot_per_base_quality_boxplot(per_base_qualities, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)`

Plot a FastQC-style per-base quality boxplot. `per_base_qualities` maps a read position (`'7'`) or position range (`'1-5'`) to a list of quality scores; each box is drawn at the range midpoint and colored on a red-to-green scale. Extra keyword arguments are passed to `Axes.boxplot` (`positions`, `widths` and `patch_artist` are fixed by the function).

**Example:**
```python
from metainformant.visualization.analysis.quality_sequencing import plot_per_base_quality_boxplot
import numpy as np

rng = np.random.default_rng(0)
per_base = {
    f'{start}-{start + 4}': rng.normal(36 - start / 10, 2, 50).tolist()
    for start in range(1, 51, 5)
}
ax = plot_per_base_quality_boxplot(per_base)
```

### `plot_adapter_content(adapter_data, *, ax=None, output_path=None, **kwargs)`

Plot adapter contamination as bars. `adapter_data` maps an adapter name to a list of content percentages (for example one per read position); the bar for each adapter shows the maximum of its list, sorted from highest to lowest and annotated with the value. Extra keyword arguments are passed to `Axes.bar` (`alpha` is fixed by the function).

**Example:**
```python
from metainformant.visualization.analysis.quality_sequencing import plot_adapter_content

adapter_data = {
    'Illumina Universal': [0.0, 0.1, 0.4, 1.2, 2.6],
    'Nextera Transposase': [0.0, 0.0, 0.1, 0.2, 0.3],
    'Small RNA 3p': [0.0, 0.0, 0.0, 0.1, 0.1],
}
ax = plot_adapter_content(adapter_data)
```

### `plot_length_distribution(length_data, *, ax=None, output_path=None, **kwargs)`

Plot a read length histogram with the mean, median and range written on the axes. `length_data` is a non-empty list or NumPy array of read lengths. Extra keyword arguments are passed to `Axes.hist` (`bins` defaults to `'auto'`; `alpha` and `edgecolor` are fixed by the function).

**Example:**
```python
from metainformant.visualization.analysis.quality_sequencing import plot_length_distribution
import numpy as np

lengths = np.random.default_rng(0).normal(150, 10, 1000).astype(int)
ax = plot_length_distribution(lengths)
```

## Other QC plots

The same modules provide further plots. Each returns a matplotlib `Axes` (the first axes for the multi-panel summaries, which draw their own figure and do not use `ax`); check the function's source for the keys it reads from its dict argument.

- `plot_gc_distribution(gc_data, *, ax=None, output_path=None, **kwargs)`: GC content histogram (list or array of percentages) with mean, median and standard deviation.
- `plot_sequence_duplication_levels(duplication_levels, *, ax=None, output_path=None, figsize=(8, 6), **kwargs)`: bars for a dict of duplication level to percentage.
- `plot_overrepresented_sequences(overrepresented_seqs, *, ax=None, output_path=None, figsize=(10, 6), **kwargs)`: horizontal bars for a list of dicts with `sequence`, `count` and `percentage`.
- `plot_kmer_profiles(kmer_counts, top_n=20, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)`: bars for the most frequent k-mers of a dict of k-mer to count.
- `plot_vcf_quality_metrics(vcf_qc_data, *, ax=None, output_path=None, figsize=(14, 10), **kwargs)`: 2 x 3 panel summary of VCF quality, depth, allele frequency, variant types, Ti/Tv ratio and summary statistics.
- `plot_protein_structure_quality(structure_quality, *, ax=None, output_path=None, figsize=(12, 8), **kwargs)`: B-factor, Ramachandran, clash score and overall quality panels.
- `plot_multiomics_quality_overview(quality_reports, *, ax=None, output_path=None, figsize=(14, 8), **kwargs)`: one quality score per omics layer.
- `plot_coverage_uniformity(coverage_data, positions=None, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)`: coverage line with mean and median reference lines.
- `plot_error_profiles(error_profiles, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)`: error rate by read position for each error type, on a log scale.
- `plot_batch_effects_qc(batch_qc_data, *, ax=None, output_path=None, figsize=(14, 10), **kwargs)`: up to six batch-effect panels.
- `plot_data_integrity_metrics(integrity_metrics, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)`: bars for the numeric entries of a metrics dict.
