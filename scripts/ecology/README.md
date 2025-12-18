# Ecology Analysis Scripts

Ecological community analysis and biodiversity assessment workflow orchestrators.

## Directory Structure

```
scripts/ecology/
├── run_ecology_analysis.py       # Ecology analysis workflow orchestrator
└── README.md                     # This file
```

## Ecology Analysis Workflow (`run_ecology_analysis.py`)

Comprehensive ecological analysis workflow orchestrator for biodiversity assessment and community ecology.

**Features:**
- Alpha diversity analysis (Shannon, Simpson, species richness)
- Beta diversity analysis (Bray-Curtis, Jaccard, UniFrac)
- Rarefaction analysis and species accumulation curves
- Community composition analysis
- Statistical testing and visualization

**Usage:**
```bash
# Diversity analysis
python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --diversity --output output/ecology/diversity

# Full community analysis
python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --diversity --beta-diversity --rarefaction

# Beta diversity only
python3 scripts/ecology/run_ecology_analysis.py --input abundance.tsv --beta-diversity
```

**Options:**
- `--input`: Input abundance table (CSV/TSV format, species x sites)
- `--output`: Output directory (defaults to output/ecology/)
- `--diversity`: Calculate alpha diversity indices
- `--beta-diversity`: Perform beta diversity analysis
- `--rarefaction`: Generate rarefaction curves
- `--group-column`: Column name for grouping samples
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/ecology/
├── diversity_indices.json        # Alpha diversity metrics
├── beta_diversity.json          # Beta diversity matrices and statistics
├── rarefaction_curves.json      # Rarefaction analysis results
├── community_composition.json   # Species composition summaries
├── diversity_plots/             # Generated visualizations
│   ├── shannon_diversity.png
│   ├── rarefaction_curves.png
│   └── beta_diversity_heatmap.png
└── analysis_report.json         # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.ecology**: Core ecological analysis functionality
- **Core utilities**: I/O, logging, path management
- **Statistical analysis**: Diversity calculations and testing
- **Visualization**: Plot generation for ecological data

## Dependencies

- **metainformant.ecology**: Ecological analysis module
- **scikit-bio/scipy**: Diversity and distance calculations
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Ecology Analysis Documentation](../../docs/ecology/README.md)
- [Core Utilities](../../docs/core/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)

