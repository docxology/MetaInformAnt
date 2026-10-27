# GWAS Visualization Gallery

This document showcases the main visualization types in the METAINFORMANT GWAS module, organized by category. Interactive, composite, fine-mapping, phenotype and map plots are listed under [Other Plot Modules](#other-plot-modules).

## Conventions

- **Inputs**: association results are passed as Python objects, not file paths. Most functions take a list of dictionaries with `CHROM`, `POS` and `p_value` keys (plus `BETA`, `SE`, `MAF` and similar where a plot needs them), and some also accept a pandas DataFrame with `CHR`, `BP` and `P` columns. Each entry below lists the keys it reads. Only `generate_all_plots` reads association results from a `.csv` or `.tsv` file.
- **Return value**: the functions documented below return a dictionary with a `status` key (`"success"`, `"failed"` with an `error` message, or `"skipped"`) and counts such as `n_variants`. Problems with the data are usually reported there instead of raised, so check `result["status"]`. The package-level `manhattan_plot` and `qq_plot` (see the notes under those entries) return a matplotlib figure.
- **Output paths**: the keyword is `output_path` in some functions and `output_file` in others, and the file extension selects the format. Some functions create the missing parent directories of the output file themselves: the entry points of `gwas.visualization.general` (the package-level `manhattan_plot` and `qq_plot` among them), `regional_plot`, `gene_annotation_plot` and `recombination_rate_plot` in `genomic/regional.py`, `kinship_dendrogram` and `kinship_clustermap`, the interactive, composite, fine-mapping and phenotype plots, and `generate_all_plots` (in `general.py` and in `interactive/suite.py`, which create their `output_dir`). Most other functions call `plt.savefig` directly, which does not create missing directories, so create the output directory first: depending on the function, a missing directory either returns `status: "failed"` or raises `FileNotFoundError`.

## Quick Start

```python
from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.gwas.visualization.interactive.suite import generate_all_plots

# Synthetic association table (use your own results in practice)
rng = np.random.default_rng(1)
table = pd.DataFrame(
    [
        {
            "ID": f"rs{chrom}_{pos}",
            "CHROM": str(chrom),
            "POS": pos,
            "p_value": float(rng.uniform()),
            "BETA": float(rng.normal(0, 0.3)),
        }
        for chrom in (1, 2, 3)
        for pos in range(100_000, 5_000_000, 50_000)
    ]
)
Path("output/gwas/results").mkdir(parents=True, exist_ok=True)
table.to_csv("output/gwas/results/association_results.tsv", sep="\t", index=False)

# Generate all visualizations at once
results = generate_all_plots(
    association_results="output/gwas/results/association_results.tsv",
    output_dir="output/gwas/plots",
    significance_threshold=5e-8
)

print(f"Generated {results['num_plots_generated']} plots")
```

`association_results` is a `.csv` or `.tsv` file (any other suffix is rejected) or a list of dictionaries. It needs a p-value column (`p_value`, `P`, `pvalue`, ...), a chromosome column (`CHROM`, `CHR`, ...) and a position column (`POS`, `BP`, ...); a `BETA` column adds the effect-size plots. Population-structure plots need the in-memory inputs described under [Visualization Suite](#visualization-suite).

## Visualization Categories

### 1. Genome-Wide Visualizations

#### Manhattan Plot
**Module**: `visualization/genomic/genome.py`  
**Function**: `manhattan_plot()`

Standard genome-wide association plot: -log10(p) against position, one point per variant.

**Input**: `results`, a list of dictionaries with `CHROM`, `POS` and `p_value` keys (`chr`, `pos` and `P` are also read); keywords `output_path`, `significance_threshold` (5e-8), `suggestiveness_threshold` (off by default) and `figsize`.

**Features**:
- Alternating chromosome colors
- Genome-wide significance line, plus an optional suggestive line
- Returns `n_variants`, `n_chromosomes` and the list of chromosomes

**Use Cases**:
- Primary GWAS results visualization
- Identifying genome-wide significant loci
- Publication figures

**Note**: this version draws one scatter artist per variant, so it suits small result sets (see [Performance Notes](#performance-notes)). The package-level `manhattan_plot` (`from metainformant.gwas.visualization import manhattan_plot`) is a separate implementation for large studies: it takes dictionaries with lowercase `chrom`, `pos` and `p_value` keys, draws all points in one call, returns a matplotlib figure, and also supports `highlight_regions`, `gene_annotations`, `label_top_n`, `colors` and a `PlotStyle` (`style=`).

**Example**:
```python-snippet
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization.genomic.genome import manhattan_plot

# Synthetic results: one dictionary per variant
rng = np.random.default_rng(1)
results = [
    {"CHROM": str(chrom), "POS": pos, "p_value": float(rng.uniform())}
    for chrom in (1, 2, 3)
    for pos in range(100_000, 5_000_000, 100_000)
]
results[10]["p_value"] = 3e-9  # one genome-wide significant hit

Path("plots").mkdir(exist_ok=True)
result = manhattan_plot(
    results,
    output_path="plots/manhattan.png",
    significance_threshold=5e-8,
    suggestiveness_threshold=1e-5
)
print(result["status"], result["n_variants"])  # success 147
```

#### Circular Manhattan Plot
**Function**: `circular_manhattan_plot()`

Circular genome-wide view with radial p-values: each chromosome gets an equal sector, the radius grows with -log10(p), and rings mark the genome-wide and suggestive thresholds.

**Input**: `results_df`, a DataFrame with `CHR`, `BP` and `P` columns or a list of dictionaries with `CHROM`, `POS` and `p_value` keys; keywords `output_file`, `significance_threshold` (5e-8), `suggestive_threshold` (1e-5, `None` removes the ring), `title` and `label_top_n` (0) to label the top hits.

**Use Cases**:
- Compact genome-wide overview
- Identifying genome-wide patterns
- Alternative to linear Manhattan plot

#### Chromosome Ideogram
**Function**: `chromosome_ideogram()`

Chromosome bars scaled by length, with the regions you choose shaded.

**Input**: `chromosome_data`, a `{chromosome: length}` dictionary or a list of result dictionaries (lengths are taken from the largest `POS` per chromosome); `highlighted_regions`, a list of `{"chrom": "1", "start": 1000000, "end": 1500000}` dictionaries (optional `color` and `alpha`) marking the loci to show; `output_file` and `title`. It does not detect significant loci itself.

**Use Cases**:
- Visualizing distribution of hits (pass them as highlighted regions)
- Quick overview of significant regions
- Identifying clustered associations

---

### 2. Statistical Diagnostics

#### Q-Q Plot
**Module**: `visualization/statistical/statistical.py`  
**Function**: `qq_plot()`

Quantile-quantile plot with 95% CI and genomic inflation factor.

**Input**: `results`, a list of p-values or a list of dictionaries with a `p_value` key (every p-value must be greater than 0; an exact 0 makes the plot fail with `math domain error`); keywords `output_path` and `figsize`.

**Features**:
- 95% confidence band from Beta-distribution quantiles (with SciPy; a rough approximation without it)
- Lambda GC calculation, shown in the title and returned as `lambda_gc`
- A warning mark in the title when lambda GC exceeds 1.1

**Interpretation**:
- Points on diagonal: well-calibrated
- Deviation at right: true associations
- Global inflation: population stratification or relatedness

**Note**: the package-level `qq_plot` (`from metainformant.gwas.visualization import qq_plot`) is a separate implementation that takes a list of p-values and returns a matplotlib figure. It labels the inflation (nominal to elevated, with a small-sample note when `n_samples` is given) and adds the mean chi-square, a Kolmogorov-Smirnov test and a GC-corrected overlay.

**Example**:
```python-snippet
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization.statistical.statistical import qq_plot

rng = np.random.default_rng(1)
p_values = [float(p) for p in rng.uniform(size=2000)]

Path("plots").mkdir(exist_ok=True)
result = qq_plot(
    p_values,
    output_path="plots/qq_plot.png"
)
print(f"lambda GC = {result['lambda_gc']:.2f}")  # close to 1 for uniform p-values
```

#### Stratified Q-Q Plot
**Function**: `qq_plot_stratified()`

Q-Q plot separated by MAF bins.

**Input**: `results`, a list of dictionaries with `p_value` and `MAF` keys, together with `maf_bins`, a list of half-open `(min, max)` tuples such as `[(0.0, 0.1), (0.1, 0.3), (0.3, 0.5)]`; or a list of p-values together with `strata`, one label per p-value. Each stratum's lambda GC is returned under `strata_results`.

**Use Cases**:
- Checking if inflation varies by frequency
- Identifying frequency-specific issues
- Quality control by allele frequency

#### Lambda GC by Chromosome
**Function**: `lambda_gc_plot()`

Genomic inflation factor per chromosome.

**Input**: `results`, a list of dictionaries with `p_value` and `CHROM` keys (a chromosome needs at least two p-values), or a list of lambda values you computed yourself; keywords `output_file` and `title`. Returns `overall_lambda_gc` and, for dictionary input, `lambda_gc_by_chrom`.

**Use Cases**:
- Identifying chromosome-specific issues
- Quality control across genome
- Detecting localized problems

#### Volcano Plot
**Function**: `volcano_plot()`

Effect size vs significance scatter plot.

**Input**: `results`, a list of dictionaries with `p_value` and lowercase `beta` keys (a missing `beta` counts as 0); keywords `output_path`, `significance_threshold` (5e-8) and `effect_size_threshold` (0.5). Returns the counts `n_significant`, `n_large_effect` and `n_significant_large_effect`.

**Use Cases**:
- Identifying large-effect variants
- Combining significance with biological relevance
- Publication figures

#### Power Plot
**Function**: `power_plot()`

Statistical power of a two-sample t-test as a contour map over sample size and effect size, with contours at 80% and 90% power.

**Input**: `sample_sizes` (per group) and `effect_sizes` (standardized, Cohen's d) as lists, `alpha` (0.05), and `output_file` (or its alias `output_path`). Requires SciPy.

**Use Cases**:
- Study design and planning
- Post-hoc power analysis
- Sample size estimation

---

### 3. Regional/Locus-Specific

#### Regional Association Plot
**Module**: `visualization/genomic/regional.py`  
**Function**: `regional_plot()`

Detailed association view for genomic window: -log10(p) for the variants in the window with a SNP-density histogram below.

**Input**: `results`, a list of dictionaries or a DataFrame (`CHROM`/`CHR`, `POS`/`BP`, `p_value`/`P`); keywords `output_file`, `chrom`, `start`, `end` (needed unless you pass `title`), `significance_threshold`, `figsize`, `title` and `recombination_rate` (a DataFrame with `CHR`, `BP` and `RATE` columns, drawn as an overlay). `chrom` must equal the stored values exactly (`"chr1"` does not match `"1"` or the integer 1); with no match the function returns `status: "success"` with `n_variants: 0` and draws nothing.

**Use Cases**:
- Fine-mapping significant loci
- Identifying lead SNPs
- Examining local association structure (for linkage disequilibrium use the regional LD plot below)

**Example**:
```python-snippet
import numpy as np
import pandas as pd

from metainformant.gwas.visualization.genomic.regional import regional_plot

# Synthetic results with CHROM, POS and p_value columns
rng = np.random.default_rng(1)
results = pd.DataFrame({"CHROM": "chr1", "POS": range(100_000, 3_000_000, 20_000)})
results["p_value"] = rng.uniform(size=len(results))

result = regional_plot(
    results,
    output_file="plots/region_chr1.png",
    chrom="chr1",
    start=1000000,
    end=2000000
)
print(result["status"], result["n_variants"])  # success 51
```

#### Regional LD Plot
**Function**: `regional_ld_plot()`

Linkage disequilibrium around a locus: the SNP positions (with the lead SNP marked) above an r² heatmap of the window.

**Input**: `vcf_path`, an uncompressed VCF file; keywords `output_file`, `chrom`, `start`, `end`, `lead_snp_pos` and `ld_threshold` (0.8). The function computes pairwise r² from the genotypes itself, so no external LD tool is needed; it returns `status: "skipped"` when the file is missing or the window holds fewer than two SNPs.

#### Gene Annotation Plot
**Function**: `gene_annotation_plot()`

Regional association plot with a gene track below it.

**Input**: the same `results`, `chrom`, `start` and `end` as `regional_plot`, plus `gene_data`, a dictionary such as `{"genes": [{"name": "geneA", "start": 1100000, "end": 1300000}]}`, and `output_file`. Without `gene_data` the lower panel reads "Gene annotations not provided".

**Use Cases**:
- Identifying variants in genes
- Functional annotation
- Publication figures

#### Recombination Rate Plot
**Function**: `recombination_rate_plot()`

Recombination rate profile along each chromosome, with a smoothed line.

**Input**: keywords `recombination_rates` (`{chromosome: [rate in cM/Mb, ...]}`) and `positions` (`{chromosome: [bp, ...]}`) plus `output_file`. Without them the function returns `status: "skipped"` (recombination map data not available).

**Note**: requires recombination map data. To overlay a rate on an association plot, pass `recombination_rate` to `regional_plot`.

---

### 4. Population Structure

#### PCA Plot
**Module**: `visualization/population/population_pca.py`  
**Function**: `pca_plot()`

Principal component scatter plot.

**Input**: `pca_data`, a dictionary with `pcs` (samples x components array), `explained_variance` (ratios) and optional `sample_names` and `populations`; or the path of a TSV file (sample ID in the first column, then `PC1`, `PC2`, ...; explained variance is then estimated from the variance of each PC column). Keywords `output_file`, `title`, `pc1` and `pc2` (0-based component indices).

**Features**:
- Optional population coloring (the `populations` key)
- Any pair of components; the 3D view is `pca_3d()` in the same module
- Sample names annotated when there are 50 samples or fewer

**Use Cases**:
- Visualizing population stratification
- Quality control (outlier detection)
- Adjusting for ancestry

**Example**:
```python-snippet
from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.gwas.visualization.population.population_pca import pca_plot

# Synthetic components: 40 samples x 5 PCs, two groups separated on PC1
rng = np.random.default_rng(2)
pcs = rng.normal(size=(40, 5))
pcs[:20, 0] += 3

Path("structure").mkdir(exist_ok=True)
components = pd.DataFrame(pcs, columns=[f"PC{i}" for i in range(1, 6)])
components.insert(0, "sample_id", [f"S{i}" for i in range(40)])
components.to_csv("structure/pca_components.tsv", sep="\t", index=False)

Path("plots").mkdir(exist_ok=True)
result = pca_plot(
    "structure/pca_components.tsv",
    output_file="plots/pca.png",
    pc1=0,  # 0-based component indices (PC1 vs PC2)
    pc2=1
)
print(result["status"], result["n_samples"])  # success 40
```

#### PCA Scree Plot
**Function**: `pca_scree_plot()`

Variance explained by each PC, as bars with a cumulative line.

**Input**: `explained_variance`, a list of ratios between 0 and 1 or the path of a text file with one value per line; keywords `output_file` and `title`.

**Use Cases**:
- Determining number of PCs to adjust for
- Understanding population structure complexity

#### Kinship Heatmap
**Module**: `visualization/population/population_admixture.py`  
**Function**: `kinship_heatmap()`

Pairwise genetic relationships, as a hierarchically clustered heatmap (seaborn `clustermap`, Ward linkage).

**Input**: `kinship_matrix` (array, list of lists, or the path of a TSV with an optional header row of sample IDs followed by numeric rows and no row-label column), then the output path as the second argument (or as `output_path_or_labels=`); `population_labels` adds color bars. Requires seaborn. `kinship_clustermap()` and `kinship_dendrogram()` in the same module are variants of this plot.

**Use Cases**:
- Identifying related samples
- Mixed model adjustment
- Quality control

#### Admixture Plot
**Function**: `admixture_plot()`

Ancestry proportions (stacked bar plot).

**Input**: `admixture_data`, a dictionary with `admixture_proportions` (samples x K) and optional `population_labels` and `ancestry_names`, or the path of an ADMIXTURE-style `.Q` file (whitespace-delimited, one row of K proportions per sample, no header); keywords `output_path`, `figsize` and `k`.

**Note**: ADMIXTURE is not run by this function; supply its output or your own ancestry proportions.

#### Kinship Dendrogram
**Function**: `kinship_dendrogram()`

Hierarchical clustering dendrogram of the samples in a kinship matrix (distance = 1 - kinship).

**Input**: `kinship_matrix` (square array or list of lists), optional `sample_labels`, `output_file`, the scipy linkage `method` (`"ward"`), and `color_by` with a `metadata` dictionary (`{sample_id: {field: value}}`) to color the labels by group. Requires SciPy.

**Use Cases**:
- Visualizing population relationships
- Identifying genetic clusters

---

### 5. Variant Properties

#### MAF Distribution
**Module**: `visualization/genomic/variants.py`  
**Function**: `maf_distribution()`

Allele frequency spectrum histogram.

**Input**: `maf_values`, a list of MAF values between 0 and 0.5 (others are clipped) or a list of dictionaries with a `MAF` key; keywords `output_file` and `title`. Returns `n_variants` and `mean_maf`.

**Features**:
- 50-bin histogram with the mean marked
- Counts per frequency class: rare (<1%), low (1-5%), common (5-10%), high (10-20%) and very high (20-50%)
- Summary statistics, including the number of rare (<1%) and common (>=5%) variants

**Use Cases**:
- Quality control
- Understanding variant frequency distribution
- Assessing QC filter effects

**Example**:
```python-snippet
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization.genomic.variants import maf_distribution

rng = np.random.default_rng(3)
maf_values = [float(m) for m in rng.beta(1, 5, size=1000) * 0.5]

Path("plots").mkdir(exist_ok=True)
result = maf_distribution(
    maf_values,
    output_file="plots/maf_dist.png"
)
print(result["n_variants"], round(result["mean_maf"], 2))  # 1000 0.08
```

#### Variant Density Plot
**Function**: `variant_density_plot()`

SNP density across genome, with the chromosomes laid end to end.

**Input**: `variants`, a list of dictionaries with `CHROM` and `POS` keys or a plain list of positions (treated as one chromosome); keywords `output_file`, `title` and `window_size` (1 Mb), which only pads each chromosome's length (it is not the density bin width).

**Use Cases**:
- Identifying high/low density regions
- Quality control
- Visualizing variant distribution

#### HWE Deviation Plot
**Function**: `hwe_deviation_plot()`

Hardy-Weinberg equilibrium p-value distribution: a histogram of -log10(p) with the p = 0.05 line, next to a Q-Q panel (its expected values are random uniform draws, so that panel changes between runs).

**Input**: `hwe_p_values`, a list of p-values or a list of dictionaries with an `HWE_P` key; keywords `output_file` and `title`.

**Use Cases**:
- Quality control
- Identifying genotyping errors
- Detecting population substructure

#### Missingness Plot
**Function**: `missingness_plot()`

Missing data patterns: per-sample and per-variant missingness histograms.

**Input**: `vcf_path`, an uncompressed VCF file (missing calls such as `./.` are counted), or a dictionary with `genotypes` (rows are variants, columns are samples, `None` or -1 for missing); `output_file` is required (omitting it raises `ValueError`). The `by_sample` keyword is accepted but not used.

**Note**: a basic VCF reader is built in; gzipped VCF files are not supported.

#### Transition/Transversion Ratio
**Function**: `transition_transversion_plot()`

Ts/Tv ratio visualization and interpretation.

**Input**: `variants`, a list of dictionaries with single-base `REF` and `ALT` keys (other variants are ignored), or a dictionary with `transitions` and `transversions` arrays; keywords `output_file` and `title`. The plot shows transition and transversion counts and a substitution matrix, with the ratio in the title; the result holds `ts_tv_ratio`, `transitions` and `transversions`.

**Interpretation**:
- Expected Ts/Tv: 2.0-2.1 (whole genome)
- Low Ts/Tv: possible quality issues

---

### 6. Effect Sizes and Functional Annotation

#### Forest Plot
**Module**: `visualization/statistical/effects.py`  
**Function**: `effect_size_forest_plot()`

Effect estimates with 95% confidence intervals.

**Input**: `effect_data`, a list of dictionaries with `CHROM`, `POS`, `BETA`, `SE` and `p_value` keys (lowercase `beta` and `se`, and `P`, are also read; the interval is the effect +/- 1.96 SE, and a missing `SE` defaults to 0.3 times the absolute effect), or a dictionary with `variants`, `effects`, `lower_ci` and `upper_ci` lists; keywords `output_file`, `title` and `top_n`.

**Features**:
- Top N most significant variants (list input is sorted by p-value)
- Confidence intervals

**Use Cases**:
- Reporting significant associations
- Comparing effect sizes
- Publication figures

**Example**:
```python-snippet
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization.statistical.effects import effect_size_forest_plot

rng = np.random.default_rng(4)
results = [
    {
        "CHROM": "1",
        "POS": 100_000 * i,
        "BETA": float(rng.normal(0, 0.3)),
        "SE": 0.1,
        "p_value": float(rng.uniform() ** 3),
    }
    for i in range(1, 41)
]

Path("plots").mkdir(exist_ok=True)
result = effect_size_forest_plot(
    results,
    output_file="plots/forest.png",
    top_n=20
)
print(result["n_variants"])  # 20
```

#### Effect Direction Plot
**Function**: `effect_direction_plot()`

Effect sizes of all variants along the genome, colored by sign (blue positive, red negative).

**Input**: `results`, a list of dictionaries with `CHROM`, `POS` and `BETA` keys (lowercase `chr`, `pos` and `beta` are also read); keywords `output_path` and `figsize`. Returns `positive_effects` and `negative_effects`.

**Use Cases**:
- Checking for directional bias
- Quality control
- Understanding effect patterns

**Note**: `visualization/genomic/regional.py` has a different function of the same name for a single region, `effect_direction_plot(effect_sizes, positions, chrom, ...)`, which returns a matplotlib figure.

#### Functional Enrichment Plot
**Function**: `functional_enrichment_plot()`

Enrichment of variant categories (functional annotation classes).

**Input**: `enrichment_results`, a dictionary of the form `{"missense": {"enrichment_ratio": 3.1, "p_value": 1e-5}, ...}` computed from your annotation; categories with a ratio of 1 or less are dropped. Keywords `output_path` and `top_n` (20). The plot shows the enrichment ratios and enrichment against -log10(p).

**Note**: requires variant annotation (ANNOVAR, VEP, SnpEff) summarized into that dictionary. The function does not read annotation files: given a list of result dictionaries instead (where `significance_threshold` applies), it does not annotate anything and sorts the significant variants into eight placeholder categories by `POS % 8`, so use the dictionary form for real data.

#### Allelic Series Plot
**Function**: `allelic_series_plot()`

Effect size and significance of each variant in a locus, ordered by position.

**Input**: `results`, a list of dictionaries with `CHROM`, `POS`, `BETA` and `p_value` keys (optional `ID` for the labels); keywords `output_path` and `gene_region`, a `(chromosome, start, end)` tuple that restricts the variants shown.

**Note**: for complex loci with multiple functional alleles.

---

### 7. Multi-Trait and Cross-Cohort Comparisons

#### Miami Plot
**Module**: `visualization/statistical/comparison.py`  
**Function**: `miami_plot()`

Back-to-back Manhattan plots for two traits (the second is mirrored below the first).

**Input**: `trait1_results` and `trait2_results`, each a list of dictionaries or a DataFrame (`CHROM`/`CHR`, `POS`/`BP`, `p_value`/`P`; a DataFrame you pass in has its columns renamed to `CHR`, `BP` and `P` in place); keywords `output_file`, `trait1_name`, `trait2_name` and `significance_threshold`.

**Use Cases**:
- Identifying pleiotropic loci
- Comparing two phenotypes
- Publication figures for trait pairs

**Example**:
```python-snippet
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization.statistical.comparison import miami_plot

rng = np.random.default_rng(5)
height_results = [
    {"CHROM": str(chrom), "POS": pos, "p_value": float(rng.uniform())}
    for chrom in (1, 2, 3)
    for pos in range(100_000, 3_000_000, 100_000)
]
bmi_results = [dict(row, p_value=float(rng.uniform())) for row in height_results]

Path("plots").mkdir(exist_ok=True)
result = miami_plot(
    height_results,
    bmi_results,
    output_file="plots/miami.png",
    trait1_name="Height",
    trait2_name="BMI"
)
print(result["status"])  # success
```

#### Multi-Trait Manhattan
**Function**: `multi_trait_manhattan()`

Manhattan plots for multiple phenotypes in a grid (up to three per row, shared y-axis).

**Input**: `trait_results`, a dictionary mapping each trait name to its results (list of dictionaries or DataFrame), with at least two traits; keywords `output_file`, `significance_threshold` and `title`.

**Use Cases**:
- Comparing many traits simultaneously
- Identifying shared associations
- Comprehensive trait analysis

#### Cross-Cohort Forest Plot
**Function**: `cross_cohort_forest()`

Forest plot of one variant's effect across cohorts.

**Input**: `cohort_results`, a dictionary of the form `{"Discovery": {"effect_size": 0.40, "ci_lower": 0.25, "ci_upper": 0.55}, ...}`; keywords `output_file`, `variant_id`, `trait_name` and `title`. Raw variant lists are not meta-analyzed: they return `status: "skipped"`, and the "overall effect" line is the unweighted mean of the cohort effect sizes.

**Note**: for comparing variant effects across studies.

#### Concordance Plot
**Function**: `concordance_plot()`

Discovery vs replication scatter plot of -log10(p) for the variants the two studies share.

**Input**: `study1_results` and `study2_results`, lists of dictionaries with an `ID` key (or `CHROM` and `POS`), `p_value` or `P`, and optionally `beta` or `BETA`; keywords `output_path`, `figsize` and `title`.

**Features**:
- Correlation calculation (returned as `correlation`; computed from the effect sizes when both studies have them, otherwise from -log10(p))
- Diagonal line of perfect concordance
- Concordance statistics: `n_overlapping_snps`, `study1_snps` and `study2_snps`

**Use Cases**:
- Validating replication
- Meta-analysis quality control
- Assessing consistency

---

## Visualization Suite

Generate all plots at once with `generate_all_plots`. Which plots appear depends on the inputs you give it:

- `association_results` alone: Manhattan and Q-Q plots, regional plots for the five most significant variants (a 1 Mb window around each), a summary panel (its PCA and kinship quadrants read "No data available", because the suite does not pass `pca_data` or `kinship_matrix` on to it), a credible-set plot and an interactive Manhattan HTML file; with a `BETA` column also an effect-size plot and an interactive volcano HTML file.
- `kinship_matrix`: kinship dendrogram and clustermap; together with `pca_data`, a population-structure panel.
- `pca_data`, a dictionary with `pcs` (samples x components), `explained_variance_ratio` and `sample_ids`: an interactive 3D PCA HTML file when there are at least three components (pass `pcs` as a list of lists; a NumPy array fails there).
- `metadata`, a `{sample_id: {"population": ..., "latitude": ..., "longitude": ...}}` dictionary: sample map and population count map.
- `genotype_matrix` (variants x samples, in the same order as the association results): LD decay and an LD heatmap of the top-hit region; with `phenotype_values`, genotype-phenotype plots for the top hits.
- `phenotype_values` with `pca_data`: phenotype-PC correlation plot; `phenotypes_dict` with two or more traits: phenotype correlation matrix.
- `finemapping_results` and `h2_data` (the dictionary returned by `partition_heritability_by_chromosome`): credible-set plot (otherwise built from the association results) and heritability bar chart.

`pca_file` and `kinship_file` are read with `numpy.load`, so a text or TSV file is skipped with a logged warning; pass `pca_data` and `kinship_matrix` instead. `vcf_file` is accepted but not used by this function.

```python
import numpy as np

from metainformant.gwas.visualization.interactive.suite import generate_all_plots

# Synthetic inputs: association results plus per-sample structure data
rng = np.random.default_rng(2)
variants = [
    {"CHROM": str(chrom), "POS": pos, "p_value": float(rng.uniform()), "BETA": float(rng.normal(0, 0.3))}
    for chrom in (1, 2, 3)
    for pos in range(100_000, 5_000_000, 100_000)
]
samples = [f"S{i}" for i in range(30)]
pca_data = {
    "pcs": rng.normal(size=(30, 5)).tolist(),  # samples x components
    "explained_variance_ratio": [0.30, 0.20, 0.10, 0.05, 0.03],
    "sample_ids": samples,
}
kinship_matrix = np.cov(rng.normal(size=(30, 200)))
metadata = {
    sample: {"population": "A" if i < 15 else "B", "latitude": 38.0 + rng.normal(), "longitude": -120.0 + rng.normal()}
    for i, sample in enumerate(samples)
}

results = generate_all_plots(
    association_results=variants,
    output_dir="output/gwas/plots",
    pca_data=pca_data,
    kinship_matrix=kinship_matrix,
    metadata=metadata,
    significance_threshold=5e-8
)

# Check results
print(f"Status: {results['status']}")
print(f"Generated: {results['num_plots_generated']}")

# List generated plots (plot name -> file path)
for plot_name, plot_path in results['plots'].items():
    print(f" {plot_name}: {plot_path}")
```

## Workflow Integration

The GWAS workflow generates its plots itself; there is no configuration key that switches the visualization step on or off (the workflow reads neither `comprehensive_plots` nor `create_plots`; the latter only affects the runtime estimate). `execute_gwas_workflow` writes its plots into the results directory (`output.results_dir`) and returns their paths in `results["outputs"]["plots"]`: Manhattan, Q-Q, PCA, kinship, effect-size and top-hit regional plots, plus strain-aware plots (strain PCA, dendrogram, clustered kinship heatmap, allele-frequency heatmap) when sample IDs allow it. `run_gwas` calls the visualization suite described above instead, forwarding the top-level `significance_threshold` key (default 5e-8), and then draws the summary and population-structure panels from its own PCA and kinship results.

```yaml
# config/gwas/gwas_amellifera.yaml
output:
  results_dir: output/gwas/amellifera/results  # results and plots are written here
```

Or programmatically:

```python
from metainformant.gwas import execute_gwas_workflow
from metainformant.gwas.data.config import load_gwas_config

config = load_gwas_config("config/gwas/gwas_amellifera.yaml")

# Validate the configuration without running anything
print(execute_gwas_workflow(config, check=True)["status"])  # validated

# Run the workflow (needs the VCF and phenotype files named in the config)
results = execute_gwas_workflow(config)
for plot_name, plot_path in results["outputs"].get("plots", {}).items():
    print(plot_name, plot_path)
```

## Performance Notes

### Large-Scale Data

None of the functions thin or downsample their input, so the cost grows with the number of variants. Timings below come from one test machine:

1. **Manhattan Plot**: `genomic.genome.manhattan_plot` draws one scatter artist per variant (about 2 s for 5,000 variants, over 20 s for 20,000), so keep it to small result sets. The package-level `manhattan_plot` draws all points in one call and rasterizes them above 100,000 variants (about 1 s for 500,000).
2. **Q-Q Plot**: both versions evaluate a Beta quantile for every point to draw the 95% band: `statistical.qq_plot` took about 4 s for 50,000 p-values, and the package-level `qq_plot` about 7 s for 100,000 and about a minute for 1,000,000. For very large studies plot a random subsample of the p-values.
3. **Circular Manhattan**: draws one scatter per chromosome (about 0.3 s for 500,000 variants), so building the input list is the main cost.
4. **Kinship Heatmap**: draws every sample (seaborn `clustermap` with Ward linkage, 800 samples in under a second) and labels the axes only when there are 50 samples or fewer; there is no automatic downsampling.

### Optimization Tips

- Filter before plotting: there is no `max_points_per_chrom` or other thinning option, so pass only the variants you need (for example everything below a p-value cutoff plus a random sample of the rest)
- Generate regional plots only for significant loci
- Use summary statistics for large cohorts
- Consider rasterization for publication figures: the package-level `manhattan_plot` rasterizes automatically above 100,000 variants, and `scatter_kwargs={"s": 8, "alpha": 0.7, "rasterized": True}` forces it

## Customization

Each function takes its own parameters (the main ones are listed under its entry); there is no shared set of matplotlib keyword arguments (`genomic.genome.manhattan_plot`, for instance, accepts only `significance_threshold`, `suggestiveness_threshold` and `figsize` besides the data and the output path). The package-level `manhattan_plot` is the most configurable: besides `colors`, `highlight_regions`, `gene_annotations` and `label_top_n` it takes a `PlotStyle`, built with `get_style` from one of the themes `"publication"`, `"presentation"` and `"poster"` plus field overrides.

```python
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization import manhattan_plot
from metainformant.gwas.visualization.config import get_style

# Results with lowercase chrom, pos and p_value keys
rng = np.random.default_rng(7)
results = [
    {"chrom": str(chrom), "pos": pos, "p_value": float(rng.uniform())}
    for chrom in (1, 2, 3)
    for pos in range(100_000, 5_000_000, 50_000)
]

Path("plots").mkdir(exist_ok=True)

# Custom colors, with the three strongest hits labeled
fig = manhattan_plot(
    results,
    output_path="plots/manhattan.png",
    colors=["#FF6B6B", "#4ECDC4"],
    label_top_n=3
)

# A preset theme with an override
fig = manhattan_plot(
    results,
    output_path="plots/manhattan_slides.png",
    style=get_style("presentation", dpi=150)
)
```

## Exporting for Publication

All plots are saved with `savefig`, so the file extension selects the format:

- Default DPI: 300 (publication quality); the package-level functions follow their `PlotStyle` (300 for `"publication"`, 150 for `"presentation"`)
- Formats: any format matplotlib supports, such as PNG, SVG and PDF
- Vector graphics for scalability

```python
from pathlib import Path

import numpy as np

from metainformant.gwas.visualization import manhattan_plot

rng = np.random.default_rng(7)
results = [
    {"chrom": str(chrom), "pos": pos, "p_value": float(rng.uniform())}
    for chrom in (1, 2, 3)
    for pos in range(100_000, 5_000_000, 50_000)
]

Path("plots").mkdir(exist_ok=True)

# Save as vector graphic: the extension selects the format
manhattan_plot(results, output_path="plots/figure.pdf")

# Functions that return a figure can also be saved by hand
fig = manhattan_plot(results)
fig.savefig("plots/figure.svg", dpi=300, bbox_inches="tight")
```

## Troubleshooting

### Common Issues

**Issue**: `generate_all_plots` returns `status: "failed"` with no plots
- The association file must be `.csv` or `.tsv` (other suffixes are rejected) with a p-value, a chromosome and a position column; the exact column names it accepts are listed under Quick Start
- The reason is logged, not raised: look for "No p-value column found" or "Unsupported file format" in the log

**Issue**: `error` reads "'str' object has no attribute 'get'" or "must be real number, not str"
- A file path was passed where a list of dictionaries (or a DataFrame) is expected; load the table first, for example `pd.read_csv(path, sep="\t").to_dict("records")`

**Issue**: `error` reads "[Errno 2] No such file or directory"
- The function does not create its output directory (see [Conventions](#conventions)); create it before plotting

**Issue**: `regional_plot` returns `n_variants: 0` and writes no file
- `chrom` does not match the values in the chromosome column exactly (`"chr1"` versus `"1"` versus the integer 1), or the window holds no variants

**Issue**: `statistical.qq_plot` fails with "math domain error"
- A p-value of exactly 0 was supplied; replace it with the smallest positive value you can represent (the package-level `qq_plot` ignores p-values outside (0, 1])

**Issue**: "Plot generation skipped"
- Some plots require additional data (PCA, kinship, genotypes, metadata, a VCF, a recombination map); see Visualization Suite and the Input line of each entry
- Check log messages for specific requirements

**Issue**: "Memory error with large datasets"
- Filter or subsample the variants before plotting (there is no thinning option)
- Use regional plots instead of genome-wide
- Increase system memory or use HPC

### Getting Help

For issues or questions:
- Check log files for detailed error messages
- See individual function docstrings for parameters
- Refer to examples in `tests/gwas/` directory

## Other Plot Modules

These live under `src/metainformant/gwas/visualization/` and are not covered entry by entry above:

- `interactive/interactive.py`: `interactive_manhattan`, `interactive_pca`, `interactive_volcano` (HTML output; Plotly when installed, otherwise a styled HTML table of the top results). `interactive_manhattan` reads each result's chromosome from `chromosome`, `chrom`, `chr`, `CHROM` or `CHR` (the first one present), its position from `position`, `pos` or `POS`, and its p-value from `p_value`, `pvalue` or `P`; a result that lacks a chromosome or a position, or a positive finite p-value, is skipped rather than placed on a default chromosome. It returns `status`, `output_path`, `n_variants`, `n_skipped` and `interactive`
- `interactive/composite.py`: `gwas_summary_panel`, `population_structure_panel`, `top_hit_detail_panel`
- `interactive/finemapping.py`: `credible_set_plot`, `conditional_analysis_plot`, `pip_vs_ld_plot`. `credible_set_plot` and the module's `compute_credible_set(assoc_results, credible_level=0.95)` take their posterior inclusion probabilities from Wakefield approximate Bayes factors: when every result has `beta` and a positive `se`, the Z-score is `beta / se` with `V = se²`; otherwise the Z-score is recovered from the two-sided `p_value` (or `P`) and the Z-score scale is used (see [Fine-Mapping & Colocalization](README.md#fine-mapping--colocalization)). `compute_credible_set` returns `pips`, `credible_set_indices` (highest PIP first), `credible_set_size` and `cumulative_probability`
- `interactive/phenotype.py`: `phenotype_distribution`, `phenotype_correlation_matrix`, `genotype_phenotype_boxplot`, `top_hits_genotype_phenotype`, `phenotype_pca_correlation`
- `genomic/ld.py` and `genomic/genome.py`: `ld_decay_plot`, `ld_heatmap_region`, `genome_wide_ld_heatmap`
- `population/geography.py`: `sample_map`, `allele_frequency_map`, `population_count_map`
- `population/population_pca.py`: `pca_multi_panel`, `pca_3d`
- `eqtl_visualization.py`: `plot_eqtl_volcano`, `plot_eqtl_boxplot`, `plot_locuszoom_eqtl`, `plot_eqtl_summary`
- `strain_plots.py` and `structure_plots.py`: strain-aware PCA, dendrogram and Fst plots, and kinship, PCA scatter, scree and loading plots

## Future Enhancements

Planned visualization additions:
- PheWAS plots (the analysis exists as `run_phewas` in `metainformant.phenotype.gwas_integration.phewas`; there is no dedicated plot)
- GWAS catalog integration (study data can be downloaded by `gwas/data/download.py`; there is no plotting integration)

---

## Summary

**Total Visualization Functions**: 30

**Categories**:
- Genome-wide: 3 functions
- Statistical: 5 functions
- Regional: 4 functions
- Population: 5 functions
- Variants: 5 functions
- Effects: 4 functions
- Comparison: 4 functions

**Quick Reference**: Use `generate_all_plots()` to create comprehensive visualization suite automatically.
