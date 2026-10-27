### Phenotype Module

The `metainformant.phenotype` module provides comprehensive phenotypic trait analysis across
multiple measurement domains: morphological, behavioral, chemical, acoustic, and electronic
tracking. It supports cross-omic integration (phenotype-genotype mapping) and configurable
analysis pipelines.

---

## Subpackages

| Subpackage       | Description                                                        | Key Classes and Functions (module)               |
|------------------|--------------------------------------------------------------------|--------------------------------------------------|
| `morphological`  | Morphometric measurements, allometric regression, shape indices    | `Measurement`, `MorphometricProfile`; `allometric_regression`, `compare_profiles`, `summary_statistics` (`profile`) |
| `behavior`       | Ethograms, behavioral sequences, time budgets, transition matrices | `Ethogram` (`ethogram`), `BehaviorSequence` (`sequence`) |
| `chemical`       | Chemical profiles (GC-MS, CHC), distance matrices, marker detection| `Compound` (`compound`); `ChemicalProfile`, `distance_matrix`, `identify_marker_compounds` (`profile`) |
| `electronic`     | RFID, video, and GPS tracking with movement ecology                | `TrackingPoint`, `Trajectory`, `detect_interactions` (`tracking`) |
| `sonic`          | Acoustic signals, FFT spectral analysis, syllable detection        | `AcousticSignal` (`signal`)                      |
| `workflow`       | Pipeline orchestration for multi-step phenotype analysis           | `PhenotypePipeline`, `PipelineConfig`, `PipelineResult` (`pipeline`); `permanova`, `pcoa` (`ecology_stats`, also exported by the package) |
| `visualization`  | Phenotype-specific plotting utilities                              | `plot_trait_distribution`, `plot_phenotype_pca`, ... (`visualization`); `plot_violin`, `plot_boxplot_with_swarm`, ... (`plots`) |
| `integration`    | Cross-omic integration (phenotype-genotype, trait-expression)      | `phenotype_genotype_association`, `trait_expression_correlation`, `multi_phenotype_integration`, `phenotype_environment_interaction` (`cross_omic`) |
| `gwas_integration` | PheWAS, heritability screening, genetic risk scores              | `run_phewas`, `phenotype_correlation_matrix`, `genetic_risk_score`, `phenotype_heritability_screen`, `categorize_phenotypes` (`phewas`) |
| `analysis`       | Life course, multivariate, and statistical analysis of phenotypes  | Modules `life_course`, `multivariate`, and `statistical` (`statistical` needs `statsmodels`) |
| `data`           | Data loading utilities                                             | Modules `antwiki` (AntWiki records) and `scraper` (AntWiki scraper) |

---

## Quick Start

### Morphological Analysis

```python
from metainformant.phenotype.morphological.measurement import Measurement
from metainformant.phenotype.morphological.profile import MorphometricProfile

# Create measurements
head_length = Measurement(name="head_length", value=1.85, unit="mm")
head_width = Measurement(name="head_width", value=1.62, unit="mm")
thorax_length = Measurement(name="thorax_length", value=2.10, unit="mm")

# Build a profile for a specimen
profile = MorphometricProfile(
    specimen_id="specimen_001",
    measurements=[head_length, head_width, thorax_length],
)

# Calculate morphometric index (e.g., Cephalic Index)
ci = profile.calculate_index("CI", numerator="head_width", denominator="head_length")
# Returns (head_width / head_length) * 100

# Geometric mean size as body-size proxy
gm_size = profile.geometric_mean_size()
```

### Behavioral Analysis

```python
from metainformant.phenotype.behavior.ethogram import Ethogram
from metainformant.phenotype.behavior.sequence import BehaviorSequence

# Define an ethogram (controlled vocabulary of behaviors)
ethogram = Ethogram({
    "F": "foraging",
    "G": "grooming",
    "R": "resting",
    "T": "trophallaxis",
})

# Validate a behavior code
ethogram.validate("F")  # True

# Look up a behavior definition
defn = ethogram.get("G")
# defn.code == "G", defn.description == "grooming"
```

### Chemical Profile Analysis

```python
from metainformant.phenotype.chemical.compound import Compound
from metainformant.phenotype.chemical.profile import ChemicalProfile
from metainformant.phenotype.chemical.profile import distance_matrix, identify_marker_compounds

# Define compounds
compound_a = Compound(name="C25", formula="C25H52", retention_time=12.3)
compound_b = Compound(name="C27", formula="C27H56", retention_time=14.1)
```

### Electronic Tracking

```python
from metainformant.phenotype.electronic.tracking import TrackingPoint, Trajectory

# Build a trajectory from tracking points (RFID, video, GPS)
points = [
    TrackingPoint(x=0.0, y=0.0, timestamp=0.0),
    TrackingPoint(x=1.0, y=1.0, timestamp=1.0),
    TrackingPoint(x=2.5, y=0.5, timestamp=2.0),
]

traj = Trajectory(entity_id="ant_42", points=points)

# Movement ecology metrics
total_dist = traj.total_distance()      # Sum of step distances
displacement = traj.net_displacement()   # Straight-line first-to-last
sinuosity_val = traj.sinuosity()         # Path tortuosity
duration = traj.duration                 # Total time span in seconds
```

### Acoustic Signal Analysis

```python
from metainformant.phenotype.sonic.signal import AcousticSignal
import numpy as np

# Create from waveform data
waveform = np.sin(2 * np.pi * 440 * np.linspace(0, 1, 44100))
signal = AcousticSignal(waveform=waveform, sample_rate=44100)

# Properties
signal.duration    # 1.0 seconds
signal.n_samples   # 44100

# Load from audio file (requires soundfile package)
# signal = AcousticSignal.from_file("recording.wav")
```

### Pipeline Orchestration

```python
from metainformant.phenotype.workflow.pipeline import PhenotypePipeline, PipelineConfig

# Configure a pipeline
config = PipelineConfig(
    name="morpho_analysis",
    phenotype_types=["morphological", "chemical"],
    input_path="data/phenotypes/",
    output_path="output/phenotype_results/",
    steps=["load", "validate", "analyze", "summarize"],
)

# Or load from YAML
# config = PipelineConfig.from_yaml("config/phenotype_pipeline.yaml")
```

### GWAS Integration

```python
import random
from metainformant.phenotype.gwas_integration.phewas import (
    run_phewas,
    phenotype_correlation_matrix,
    genetic_risk_score,
    phenotype_heritability_screen,
    categorize_phenotypes,
)
from metainformant.simulation.models.popgen import generate_genotype_matrix
from metainformant.gwas.analysis.structure import compute_kinship_matrix

rng = random.Random(0)
n = 40  # samples; a phenotype needs at least 3 non-missing samples to be tested
genotype = [rng.choice([0, 1, 2]) for _ in range(n)]  # additive dosage of one variant, one value per sample
phenotypes = {
    "body_length": [5.0 + 0.8 * g + rng.gauss(0, 0.5) for g in genotype],
    "type_2_diabetes": [rng.choice([0, 1]) for _ in range(n)],
}

# PheWAS: one variant against every phenotype, sorted by p_value
hits = run_phewas(genotype, phenotypes)
print(hits[0]["phenotype"], hits[0]["beta"], hits[0]["p_value"])  # body_length comes first

# Pairwise Pearson correlations between phenotypes
corr = phenotype_correlation_matrix(phenotypes)
print(corr["phenotype_names"], corr["correlation_matrix"])

# Genetic risk score: samples x variants dosages plus one effect size per variant
scores = genetic_risk_score([[0, 1, 2], [1, 1, 0], [2, 0, 1]], effect_sizes=[0.2, -0.1, 0.4])
print(scores["risk_scores"])  # one score per sample

# Heritability screen (Haseman-Elston) against a samples x samples kinship matrix
genotypes = generate_genotype_matrix(n_individuals=n, n_sites=200, rng=rng)
kinship = compute_kinship_matrix(genotypes)["kinship_matrix"]
h2 = phenotype_heritability_screen(phenotypes, kinship)  # dicts with phenotype, h2, se, p_value

# Group phenotype names into categories by keyword
categories = categorize_phenotypes(phenotypes)
print(categories["categories"])  # {'endocrine': ['type_2_diabetes']}
```

Each `run_phewas` result holds `phenotype`, `beta`, `se`, `p_value`, `n_samples` and `phenotype_category`; the last is always `"uncategorized"`, so use `categorize_phenotypes` to group the names. The p-values are two-sided Student's t with `n_samples - 2` degrees of freedom (`scipy.stats.t`), exact at any sample size. A phenotype with fewer than 3 non-missing samples is skipped.

### Cross-omic Integration

```python
from metainformant.phenotype.integration.cross_omic import (
    phenotype_genotype_association,
    trait_expression_correlation,
)

traits = {f"s{i}": v for i, v in enumerate([1.2, 0.8, 2.5, 1.9, 0.4, 2.2, 1.1, 0.9])}
expression = {"geneA": {f"s{i}": v for i, v in enumerate([0, 0, 5.1, 3.3, 0, 4.0, 0, 0])}}

# Spearman correlation of a trait with each gene's expression
result = trait_expression_correlation(traits, expression, method="spearman")
print(result["correlations"]["geneA"])  # correlation, p_value, n_samples

# One variant (0/1/2 dosage per sample) against the trait
genotypes = {"var1": [0, 1, 2, 1, 0, 2, 1, 0]}
association = phenotype_genotype_association(traits, genotypes, method="linear_regression")
print(association["associations"]["var1"])  # beta, intercept, r_squared, p_value, n_valid
```

The p-values of `phenotype_genotype_association`, `trait_expression_correlation` and `phenotype_environment_interaction` come from `scipy.stats`: a Pearson correlation uses the t statistic `r * sqrt((n - 2) / (1 - r**2))` with `n - 2` degrees of freedom (the test of `scipy.stats.pearsonr`), and the `"linear_regression"` method uses the F distribution with 1 and `n - 2` degrees of freedom (the test of the slope in `scipy.stats.linregress`). `method="spearman"` correlates ranks with tied values sharing the mean of the ranks they span, so it equals `scipy.stats.spearmanr` on data with ties, such as expression with many zeros, and does not depend on the order of the samples.

---

## Unit Conversion

The `Measurement` class converts between the metric length units `um`, `mm`, `cm` and `m`, in any direction. The micrometre may also be written `µm` (micro sign) or `μm` (Greek mu). Unit symbols are case-sensitive.

```python
from metainformant.phenotype.morphological.measurement import Measurement

m = Measurement(name="leg_length", value=3.5, unit="mm")
m_cm = m.convert("cm")   # Measurement(name='leg_length', value=0.35, unit='cm', ...)
```

The conversion is exact arithmetic: the factor between two units is a ratio of whole numbers of micrometres, and the value and its `uncertainty` are each rounded to a float once, so `Measurement("x", 1.0, "mm").convert("um").value` is exactly `1000.0` rather than `1000.0000000000001`. Converting to the measurement's own unit returns the measurement itself. A measurement whose unit, or a `target_unit`, is not one of these raises `ValueError`.

---

## See Also

- **[AntWiki Integration](antwiki.md)** -- AntWiki data retrieval
- **[Life Course Analysis](../phenotype/life_course.md)** -- Temporal trajectory analysis
- **[Phenotype Visualization](../phenotype/visualization.md)** -- Plotting tools
