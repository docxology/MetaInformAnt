# Specification: phenotype

## Scope

Phenotype module for MetaInformAnt. Multi-modal phenotyping across morphological,
behavioral, chemical, acoustic, and electronic tracking domains.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: analysis, behavior, chemical, data, electronic, gwas_integration, integration, morphological, sonic, visualization, workflow
- **Key Concepts**: Morphology, behavior, chemical profiles, acoustic, electronic tracking

## API Definition

### Exports — `analysis/life_course.py`

- `extract_phenotypes_from_events` — Extract phenotype features from event sequences (life_events or local)
- `aggregate_temporal_phenotypes` — Time-windowed phenotype aggregation
- `map_events_to_traits` — Map events to trait definitions with counts/timestamps

### Exports — `morphological/`

- `Measurement` — Dataclass for morphometric measurements with unit conversion (mm↔cm↔m↔µm)
- `MorphometricProfile` — Collection of measurements for a specimen

### Exports — `analysis/statistical.py`

- `perform_linear_regression` — Simple OLS regression summary (slope, r², p)
- `perform_multifactor_anova` — Multi-factor ANOVA table via OLS (Type 2, Type 1 fallback)
- `get_comprehensive_pairwise_ttests` — Pairwise Welch's t-tests sorted by p-value
- `calculate_summary_stats` — Descriptive statistics by group
- `perform_anova` / `perform_kruskal` / `perform_ttest` — One-way parametric/non-parametric/two-group tests
- `correlate_phenotypes` — Pearson correlation matrix for numeric phenotype columns

### Exports — `mappings.py`

- `BIOLOGICAL_GROUP_MAP` / `PHENOTYPE_LINK_MAP` — BeeWAS biological-group labels and W/Q phenotype axis
- `STRAIN_FULL_NAMES` / `STRAIN_PALETTE` / `STRAIN_ORDER` — Strain metadata for C/I/M/R strains
- `map_biological_groups` / `map_phenotype_links` — DataFrame column mappers (in-place)

### Exports — `workflow/ecology_stats.py`

- `pcoa` / `permanova` — Re-exports of the ecology domain's ordination and PERMANOVA implementations
