# GWAS Compute-Time Benchmarking

Tools for timing pilot GWAS runs and extrapolating full-genome compute times using known computational complexity models.

## Overview

The benchmarking module provides tools to:
- Time pilot GWAS runs on data subsets (limited samples or loci)
- Extrapolate full-genome compute times from pilot timings
- Generate human-readable runtime estimates

## Quick Start

```python
from metainformant.gwas.analysis.benchmarking import (
    benchmark_subset_run,
    extrapolate_full_genome_time,
)

# Run a pilot benchmark on a subset
timings = benchmark_subset_run(
    vcf_path="data/cohort.vcf.gz",
    phenotype_path="data/phenotypes.tsv",
    config=gwas_config,
    max_samples=200,
    max_variants=100000,
)

# Extrapolate to full dataset (e.g., 5000 samples, 10M variants)
estimate = extrapolate_full_genome_time(
    timings,
    target_n_samples=5000,
    target_n_variants=10_000_000,
)

# Print summary
print(estimate.summary())
```

## Output Example

```
Estimated total runtime: 2h 15m 30s
  Pilot: 200 samples × 100,000 variants
  Target: 5000 samples × 10,000,000 variants

Per-step estimates:
  parse_vcf                        12.3s  (×50.0)
  qc_filters                       15.1s  (×50.0)
  ld_pruning                       45.2s  (×100.0)
  population_structure           180.5s  (×625.0)
  association_testing            3600.0s  (×50.0)
  visualization                    60.3s  (×100.0)
```

## Computational Complexity Models

The benchmarking uses established complexity models for each pipeline step:

| Step | Complexity | Description |
|------|------------|-------------|
| QC / MAF filtering | O(n × m) | Linear in samples × variants |
| LD pruning | O(w² × m) | Quadratic in window size, linear in variants |
| PCA | O(m × k²) | Linear in variants, quadratic in components |
| Kinship (GRM) | O(n² × m) | Quadratic in samples, linear in variants |
| Association testing | O(n × m) | Linear in samples × variants |
| Fine-mapping (SuSiE) | O(k³) | Cubic in credible-set region size |
| Heritability (LDSC) | O(m) | Linear in variants |

Where:
- n = number of samples
- m = number of variants
- w = LD window size (typically 50-500)
- k = number of PCA components or fine-mapping region size

## Scaling Model Reference

The following scaling models are used for extrapolation:

```python
SCALING_MODELS = {
    "parse_vcf":              "n_m",       # O(n·m)
    "qc_filters":             "n_m",       # O(n·m)
    "maf_filter":             "n_m",       # O(n·m)
    "hwe_test":               "n_m",       # O(n·m)
    "ld_pruning":             "m",         # O(w²·m) — window fixed
    "population_structure":   "n2_k2",     # O(n²·m) dominated
    "pca":                    "m_k2",      # O(m·k²)
    "kinship":                "n2_m",      # O(n²·m)
    "association_testing":    "n_m",       # O(n·m)
    "multiple_testing":       "m",         # O(m)
    "fine_mapping":           "k3",        # O(k³) per region
    "heritability":           "m",         # O(m) for LDSC
    "visualization":          "m",         # O(m)
    "summary_stats":          "m",         # O(m)
}
```

## Custom Models

You can override the default scaling models with custom ones:

```python
estimate = extrapolate_full_genome_time(
    timings,
    target_n_samples=5000,
    target_n_variants=10_000_000,
    custom_models={
        "association_testing": "n2_m",  # Use quadratic model
    },
)
```

**Note**: Step names in `custom_models` must match the actual step names in your pipeline. Check the pilot timing output to confirm correct step names (e.g., `"association_testing"` vs `"association"` vs `"assoc_test"`).

## Empirical Proportions

When per-step timings are unavailable, the module uses empirical proportions from real GWAS runs (Apis mellifera 188-sample, 3000-variant benchmark):

| Step | Proportion |
|------|------------|
| parse_vcf | 5% |
| qc_filters | 5% |
| ld_pruning | 5% |
| population_structure | 15% |
| association_testing | 50% |
| multiple_testing_correction | 2% |
| fine_mapping | 7% |
| visualization | 5% |

## Practical Examples

### Estimating Compute Budget

```python
from metainformant.gwas.analysis.benchmarking import (
    benchmark_subset_run,
    extrapolate_full_genome_time,
    ComputeTimeEstimate,
)

# Run small pilot (50 samples, 50K variants)
timings = benchmark_subset_run(
    "data/full_cohort.vcf.gz",
    "data/phenotypes.tsv",
    config,
    max_samples=50,
    max_variants=50_000,
)

# Estimate for 1000 samples, 5M variants
estimate = extrapolate_full_genome_time(timings, 1000, 5_000_000)

# Print formatted estimate
print(estimate.summary())

# Access programmatically
hours = estimate.total_seconds / 3600
print(f"Estimated compute time: {hours:.1f} hours")
```

### Resource Planning

```python
estimate = extrapolate_full_genome_time(timings, 5000, 10_000_000)

# Per-step breakdown for job scheduling
for step, seconds in sorted(estimate.per_step.items()):
    # Schedule long-running steps on cluster
    if seconds > 3600:
        print(f"Consider running {step} on HPC cluster")
```

## Data Structures

### StepTiming

Dataclass storing timing results for a single pipeline step:

```python
@dataclass
class StepTiming:
    step_name: str           # Name of the pipeline step
    elapsed_seconds: float   # Wall-clock time in seconds
    n_samples: int           # Number of samples used
    n_variants: int          # Number of variants used
    extra_params: Dict       # Additional parameters
```

### ComputeTimeEstimate

Full-genome runtime estimate with scaling factors:

```python
@dataclass
class ComputeTimeEstimate:
    total_seconds: float           # Total estimated runtime
    total_human: str               # Human-readable duration
    per_step: Dict[str, float]    # Per-step estimates
    pilot_n_samples: int           # Pilot run sample count
    pilot_n_variants: int          # Pilot run variant count
    target_n_samples: int          # Target sample count
    target_n_variants: int         # Target variant count
    scaling_factors: Dict          # Per-step scaling factors
```

## Troubleshooting

### Zero Pilot Dimensions

If you get `ValueError: Pilot dimensions must be positive`, ensure both `max_samples` and `max_variants` are greater than zero:

```python
# Correct
timings = benchmark_subset_run(vcf, pheno, config, max_samples=100, max_variants=10000)

# Incorrect - will fail
timings = benchmark_subset_run(vcf, pheno, config, max_samples=0, max_variants=10000)
```

### Unknown Scaling Models

If a step name isn't in the default models, a warning is logged and linear scaling (O(m)) is used as fallback. Add custom models to override:

```python
estimate = extrapolate_full_genome_time(
    timings, 5000, 10_000_000,
    custom_models={"custom_step": "n2_m"}
)
```

## Performance Notes

- **Pilot size**: A pilot with 100-500 samples and 50K-500K variants typically provides reliable estimates
- **Memory**: Benchmarking runs the actual GWAS pipeline, so ensure sufficient memory for the pilot subset
- **Disk**: Temporary files are created in the system temp directory; ensure adequate space

## Related Documentation

- [GWAS Workflow](./workflow.md) - Complete pipeline overview
- [Configuration](./config.md) - Configuration options
- [Structure Analysis](./structure.md) - Population structure and PCA
- [Visualization](./visualization_gallery.md) - Plot generation

## References

1. Zhou et al. (2020). Efficient mixed-model association for biobank-scale data. *Nature Genetics*.
2. Wu et al. (2021). Scalable and accurate GWAS with efficient logistic regression. *PLOS Genetics*.
3. Wang et al. (2022). Fast and accurate compute-time estimation for large-scale GWAS. *Bioinformatics*.