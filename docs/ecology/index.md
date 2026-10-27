### Ecology: Overview

Ecological analysis including community composition, biodiversity metrics, ordination, indicator-species tests, functional and phylogenetic diversity, and macroecology. The package is `metainformant.ecology`, with the subpackages `analysis` (the modules `community`, `functional`, `indicators`, `macroecology` and `ordination`), `phylogenetic`, `traits` (an empty placeholder namespace) and `visualization`.

## Key Components

### Community Ecology Analysis
- **Diversity metrics**: Shannon, Simpson, and inverse Simpson indices
- **Species richness**: Observed counts and the Chao1 estimate of total richness
- **Evenness measures**: Pielou's evenness
- **Rarefaction analysis**: Standardization of sampling effort

### Community Composition
- **Similarity measures**: Bray-Curtis, Jaccard, and Sorensen (`community.beta_diversity`), plus Euclidean, Manhattan, and Canberra distances (`ordination.distance_matrix`)
- **Ordination methods**: PCoA, NMDS, CCA, and Procrustes rotation for community visualization
- **Clustering analysis**: Hierarchical clustering of communities (UPGMA, single, or complete linkage)
- **Beta diversity**: Analysis of variation between communities (alpha-beta-gamma partitioning, PERMANOVA, ANOSIM, PERMDISP)

### Environmental Integration
- **Environmental variables**: Constrained ordination (`ordination.cca`) relates community composition to a matrix of measured environmental variables
- **Group tests**: `indicators.indval`, `permanova`, `anosim`, and `simper` compare composition across categorical environments such as habitat or treatment
- **Biotic interactions**: Not modeled in `metainformant.ecology`; `metainformant.simulation.models.agents` simulates predator-prey and competition in agent-based ecosystems

### Trait, Phylogenetic, and Large-Scale Analysis
- **Functional diversity** (`analysis.functional`): FRic, FEve, FDiv, FDis, Rao's Q, and community weighted means
- **Phylogenetic diversity** (`phylogenetic.diversity`): Faith's PD, UniFrac, NRI/NTI, phylogenetic signal, and tree building from a distance matrix
- **Macroecology** (`analysis.macroecology`): Species abundance distribution fits, species-area and distance-decay relationships, occupancy-frequency analysis, metabolic scaling, endemism indices, and Taylor's power law

## Analysis Workflows

### Biodiversity Assessment
```python
from metainformant.ecology.analysis import community

community_data = [12, 7, 3, 1]

# Calculate diversity metrics
shannon_diversity = community.shannon_diversity(community_data)
simpson_diversity = community.simpson_diversity(community_data)
species_richness = community.species_richness(community_data)
metrics = community.community_metrics(community_data)
```

### Community Comparison
```python
from metainformant.ecology.analysis import community, indicators, ordination

# Three sites (rows) by four species (columns)
community_matrix = [[10, 5, 0, 3], [0, 6, 8, 2], [7, 7, 7, 7]]
site1_data, site2_data = community_matrix[0], community_matrix[1]

# Compare communities using different metrics
bray_curtis = community.beta_diversity(site1_data, site2_data, method="bray_curtis")
jaccard = community.beta_diversity(site1_data, site2_data, method="jaccard")

# Ordination analysis (NMDS needs more sites than axes; the default is two axes)
distance = ordination.distance_matrix(community_matrix, method="bray_curtis")
ordination_result = ordination.nmds(distance)

# Clustering of communities (methods: "upgma", "single", "complete")
clusters = indicators.cluster_communities(distance, method="upgma", n_clusters=2)
```

### Environmental Correlations
Constrained ordination relates community composition to measured environmental variables, and the indicator and group tests compare composition across categorical environments.

```python
from metainformant.ecology.analysis import indicators, ordination

community_matrix = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0], [5, 4, 2]]  # sites x species
environment = [[1.0, 20.0], [1.5, 22.0], [3.0, 12.0], [3.5, 11.0], [2.0, 16.0]]  # sites x variables
environmental_categories = ["forest", "forest", "grassland", "grassland", "forest"]  # one label per site

# Species and sites in the space constrained by the environmental variables
constrained = ordination.cca(community_matrix, environment)
print(constrained["variance_explained"])  # one fraction per constrained axis

# Indicator values of each species for the categorical environment
indicator_values = indicators.indval(community_matrix, environmental_categories)
```

## Integration with Other Modules

### With Phenotype Data
```python
from metainformant.ecology.analysis import functional
from metainformant.phenotype.data import antwiki

# AntWiki records carry numeric morphology (antwiki.load_antwiki_json reads them from a JSON file)
raw_records = [
    {"species_name": "rufa", "genus": "Formica", "morphology": {"body_length_mm": 7.5, "head_width_mm": 1.9}},
    {"species_name": "fusca", "genus": "Formica", "morphology": {"body_length_mm": 6.0, "head_width_mm": 1.5}},
    {"species_name": "niger", "genus": "Lasius", "morphology": {"body_length_mm": 4.0, "head_width_mm": 1.1}},
]
records = [antwiki.AntWikiRecord(raw) for raw in raw_records]

# Combine ecological and phenotypic data: a species-by-trait matrix plus one abundance per species
species, trait_names, trait_matrix = antwiki.create_phenotype_matrix(records, ["body_length", "head_width"])
abundances = [20, 12, 5]
trait_summary = functional.functional_diversity_suite(trait_matrix, abundances)
```

Species missing a trait value come back as `None` in the matrix; drop them first (for example with `antwiki.filter_antwiki_records(records, required_phenotypes=[...])`).

### With Statistical Analysis
```python
from metainformant.ecology.analysis import community

abundances = [45, 23, 12, 8, 5]
community1 = [10, 5, 0, 3]
community2 = [0, 6, 8, 2]

# Calculate diversity indices
shannon = community.shannon_diversity(abundances)
simpson = community.simpson_diversity(abundances)

# Community comparison
beta_diversity = community.beta_diversity(community1, community2, method="bray_curtis")
```

## Planned Extensions

Future ecological analysis modules will include:
- **Conservation biology**: Biodiversity assessment and management
- **Landscape ecology**: Spatial analysis and connectivity
- **Ecosystem modeling**: Dynamic ecosystem simulations

Macroecology and functional ecology are already available as `analysis.macroecology` and `analysis.functional`.
