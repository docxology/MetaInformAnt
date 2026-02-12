### Phenotype Module

The `metainformant.phenotype` module provides comprehensive phenotypic trait analysis across
multiple measurement domains: morphological, behavioral, chemical, acoustic, and electronic
tracking. It supports cross-omic integration (phenotype-genotype mapping) and configurable
analysis pipelines.

---

## Subpackages

| Subpackage       | Description                                                        | Key Exports                                      |
|------------------|--------------------------------------------------------------------|--------------------------------------------------|
| `morphological`  | Morphometric measurements, allometric regression, shape indices    | `Measurement`, `MorphometricProfile`             |
| `behavior`       | Ethograms, behavioral sequences, time budgets, transition matrices | `Ethogram`, `BehaviorSequence`                   |
| `chemical`       | Chemical profiles (GC-MS, CHC), distance matrices, marker detection| `Compound`, `ChemicalProfile`                    |
| `electronic`     | RFID, video, and GPS tracking with movement ecology                | `TrackingPoint`, `Trajectory`                    |
| `sonic`          | Acoustic signals, FFT spectral analysis, syllable detection        | `AcousticSignal`                                 |
| `workflow`       | Pipeline orchestration for multi-step phenotype analysis           | `PhenotypePipeline`, `PipelineConfig`            |
| `visualization`  | Phenotype-specific plotting utilities                              | (domain-specific plot functions)                 |
| `integration`    | Cross-omic integration (phenotype-genotype, trait-expression)      | (integration functions)                          |
| `gwas_integration` | PheWAS, heritability screening, genetic risk scores              | `run_phewas`, `genetic_risk_score`               |
| `analysis`       | Life course and temporal trajectory analysis                       | (temporal analysis functions)                    |
| `data`           | Data loading utilities                                             | (loader functions)                               |

---

## Quick Start

### Morphological Analysis

```python
from metainformant.phenotype import Measurement, MorphometricProfile

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
from metainformant.phenotype import Ethogram, BehaviorSequence

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
from metainformant.phenotype import Compound, ChemicalProfile
from metainformant.phenotype.chemical import distance_matrix, identify_marker_compounds

# Define compounds
compound_a = Compound(name="C25", formula="C25H52", retention_time=12.3)
compound_b = Compound(name="C27", formula="C27H56", retention_time=14.1)
```

### Electronic Tracking

```python
from metainformant.phenotype import TrackingPoint, Trajectory

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
from metainformant.phenotype import AcousticSignal
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
from metainformant.phenotype import PhenotypePipeline, PipelineConfig

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
from metainformant.phenotype import (
    run_phewas,
    phenotype_correlation_matrix,
    genetic_risk_score,
    phenotype_heritability_screen,
    categorize_phenotypes,
)
```

---

## Unit Conversion

The `Measurement` class supports automatic unit conversion for metric length units
(um, mm, cm, m):

```python
m = Measurement(name="leg_length", value=3.5, unit="mm")
m_cm = m.convert("cm")   # Measurement(name="leg_length", value=0.35, unit="cm")
```

---

## See Also

- **[AntWiki Integration](antwiki.md)** -- AntWiki data retrieval
- **[Life Course Analysis](../phenotype/life_course.md)** -- Temporal trajectory analysis
- **[Phenotype Visualization](../phenotype/visualization.md)** -- Plotting tools
