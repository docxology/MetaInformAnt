# Metabolomics Module

Metabolomics analysis: mass spectrometry data processing, metabolite identification, pathway mapping, and metabolite-gene integration.

## Overview

Metabolomics analysis: mass spectrometry data processing, metabolite identification, pathway mapping, and metabolite-gene integration.


## Table of Contents

- [Architecture](#architecture)
- [Submodules](#submodules)
- [Quick Start](#quick-start)
- [Related](#related)

## Architecture

```mermaid
graph TD
    subgraph "Metabolomics Module"
        IO[io/] --> |formats.py| FMT[CSV / MGF I/O]

        AN[analysis/] --> |identification.py| ID[Metabolite Identification & Quantification]

        PW[pathways/] --> |enrichment.py| PE[Metabolite Set Enrichment]

        VZ[visualization/] --> |plots.py| PL[Volcano / PCA / Peaks]
    end

    IO --> AN
    AN --> PW
    AN --> VZ
```

## Submodules

| Module | Purpose |
|--------|---------|
| [`io/`](io/) | Intensity-matrix CSV I/O and MGF (Mascot Generic Format) spectrum reading/writing, spectra filtering, extracted-ion chromatograms |
| [`analysis/`](analysis/) | Metabolite identification (`identify_metabolites`, `identify_with_adducts`), normalization, fold change, differential abundance, imputation, spectral similarity |
| [`pathways/`](pathways/) | Over-representation analysis (`metabolite_set_enrichment`), Benjamini-Hochberg FDR, pathway activity scoring |
| [`visualization/`](visualization/) | Volcano plot data, PCA ordination, heatmap data, chromatographic peak detection, spectrum plot prep, RT alignment |

## Quick Start

```python
import numpy as np
from metainformant.metabolomics.io import formats
from metainformant.metabolomics.analysis import identification
from metainformant.metabolomics.pathways import enrichment

# Load an intensity matrix (first column = metabolite IDs, first row = samples)
data = formats.read_csv("metabolomics_data.csv")

# Identify metabolites by matching observed m/z against a reference of exact masses
database = {"glucose": 180.0634, "citrate": 192.0270, "lactate": 90.0317}
identified = identification.identify_metabolites(
    np.array([180.0641, 192.0279]), database, ppm_tolerance=10.0
)

# Over-representation analysis against a pathway database
pathway_db = {"glycolysis": ["glucose", "lactate"], "TCA_cycle": ["citrate"]}
enriched = enrichment.metabolite_set_enrichment(["glucose", "lactate"], pathway_db)
```

## Input handling

The readers and matchers report bad input instead of repairing it silently.

- `formats.read_csv`: a blank cell or `NA` / `N/A` / `#N/A` / `NaN` / `null` / `None` (any case) is returned as
  `NaN` (never `0.0`); any other non-numeric cell, including a decimal comma, raises `ValueError` naming the line,
  metabolite and sample, and a row whose width differs from the header raises `ValueError` naming the line.
  Columns before `data_start_col` other than `metabolite_col` fill `MetabolomicsDataset.metadata` as
  `{metabolite: {column: value}}`.
  `missing_value_imputation` treats `NaN` and zero as missing, so the two compose. Impute before
  `normalize_intensities(method="total_ion_count" | "median")`, which raises on `NaN` because one missing value
  would distort the scale of its whole sample; `method="log2"` keeps a `NaN` in place and `fold_change` returns
  `NaN` (not a floored value) for a metabolite with a missing value in either group.
- `formats.read_mgf`: peak lines are read whatever their first character (`.5 99`, `+200.0 20`), comment lines
  (`#`, `;`, `!`, `/`) are skipped, and an empty or non-numeric `PEPMASS`, `RTINSECONDS` or `SCANS` value or a
  malformed peak line raises `ValueError` naming the spectrum. A `BEGIN IONS` block that is never closed by
  `END IONS` (for example the last spectrum of a truncated file) is dropped without an error.
- `identification.identify_metabolites` / `identify_with_adducts`: `ppm_tolerance` must be a positive finite
  number; in `identify_with_adducts`, `ion_mode` is `"positive"` or `"negative"` in any case and anything else
  raises `ValueError`, even when `adducts` is given.
- `enrichment.metabolite_set_enrichment`: a `background_size` smaller than the query or a pathway raises
  `ValueError`, because it cannot be a universe that contains them.
- `visualization.plots.detect_chromatographic_peaks` integrates peaks with the trapezoidal rule and runs on
  every supported NumPy (>= 1.24).

## Related
- [API Reference](SPEC.md) — Type signatures, error codes, data structures

- [metainformant.multiomics](../multiomics/) - Multi-omic integration
- [metainformant.protein](../protein/) - Protein-metabolite interactions
- [metainformant.visualization](../visualization/) - General plotting
