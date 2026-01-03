# Population Genetics Scripts

Population genetics analysis and synthetic data generation workflow orchestrators.

## Directory Structure

```
scripts/popgen/
├── __init__.py                    # Package initialization
├── analysis.py                    # Main workflow orchestrator ⭐
├── generate_dataset.py            # Dataset generation functions
├── analyze.py                     # Statistical analysis functions
├── visualize.py                   # Visualization functions
└── report.py                      # Reporting functions
```

## Population Genetics Analysis (`analysis.py`)

Comprehensive population genetics workflow orchestrator for synthetic data generation and full statistical analysis pipeline.

**Features:**
- Large-scale synthetic dataset generation
- Multiple demographic scenarios (bottlenecks, expansions, migrations)
- Complete statistical analysis suite
- Neutrality tests and demographic inference
- Comprehensive visualization suite

**Usage:**
```bash
# Run complete population genetics analysis
python3 scripts/popgen/analysis.py
```

**What it generates:**
- Synthetic populations with realistic demographic histories
- Genome-wide SNP data across multiple populations
- Complete statistical analysis (FST, π, θ, Tajima's D)
- Neutrality test suites and demographic comparisons
- 20+ comprehensive visualizations

**Output Structure:**
```
output/popgen_comprehensive/
├── synthetic_data/                # Generated synthetic datasets
│   ├── populations.json
│   ├── genotypes.json
│   └── demographic_scenarios.json
├── summary_statistics/            # Population genetics statistics
│   ├── fst_values.json
│   ├── nucleotide_diversity.json
│   ├── tajimas_d.json
│   └── neutrality_tests.json
├── demographic_analysis/          # Demographic inference results
│   ├── population_size_history.json
│   ├── migration_rates.json
│   └── bottleneck_detection.json
├── comparative_analysis/          # Cross-population comparisons
│   ├── population_differences.json
│   └── statistical_comparisons.json
├── visualizations/                # Generated plots (20+ types)
│   ├── allele_frequency_spectra.png
│   ├── demographic_histories.png
│   ├── fst_matrices.png
│   ├── neutrality_test_distributions.png
│   └── pca_population_structure.png
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.dna.population**: Core population genetics functionality
- **metainformant.dna.population_analysis**: Statistical analysis suite
- **metainformant.dna.population_viz**: Visualization utilities
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.dna**: DNA and population genetics modules
- **NumPy/SciPy**: Statistical computing
- **pandas**: Data manipulation
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Population Genetics Documentation](../../docs/dna/population_genetics.md)
- [Synthetic Data Generation](../../docs/dna/synthetic_data.md)
- [Statistical Analysis](../../docs/dna/statistical_analysis.md)
- [METAINFORMANT CLI](../../docs/cli.md)

