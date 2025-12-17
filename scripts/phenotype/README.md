# Phenotype Analysis Scripts

Phenotypic trait analysis and biological data curation workflow orchestrators.

## Directory Structure

```
scripts/phenotype/
├── run_phenotype_analysis.py      # Main phenotype analysis orchestrator ⭐
├── scrape_antwiki.py              # AntWiki species data scraper
├── load_antwiki_example.py        # AntWiki data loading examples
├── test_scraper_cloudscraper.py   # Scraper testing utilities
└── README.md                      # This file
```

## Phenotype Analysis Workflow (`run_phenotype_analysis.py`)

Comprehensive phenotype analysis workflow orchestrator for trait data analysis, statistics, and correlation studies.

**Features:**
- Phenotype data loading and preprocessing
- Statistical analysis of trait distributions
- Correlation analysis between traits
- Data quality assessment and validation
- Integration with biological databases

**Usage:**
```bash
# Load AntWiki JSON data
python3 scripts/phenotype/run_phenotype_analysis.py --input antwiki_species.json --output output/phenotype/antwiki

# Analyze CSV/TSV phenotype data
python3 scripts/phenotype/run_phenotype_analysis.py --input traits.csv --analyze-statistics --analyze-correlations

# Basic statistics only
python3 scripts/phenotype/run_phenotype_analysis.py --input traits.tsv --analyze-statistics
```

## AntWiki Data Acquisition

### AntWiki Scraper (`scrape_antwiki.py`)

Command-line tool for scraping comprehensive species phenotype data from AntWiki.org.

**Features:**
- Automated species page scraping
- Comprehensive phenotype extraction
- Rate limiting and error handling
- Batch processing capabilities

**Usage:**
```bash
# Scrape specific species
python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus

# Scrape multiple species with delay
python3 scripts/phenotype/scrape_antwiki.py --all --limit 10 --delay 3.0

# Batch processing with custom output
python3 scripts/phenotype/scrape_antwiki.py --all --output output/phenotype/antwiki/
```

### AntWiki Loading Examples (`load_antwiki_example.py`)

Examples and utilities for loading and processing AntWiki phenotype data.

**Usage:**
```bash
# Run loading examples
python3 scripts/phenotype/load_antwiki_example.py
```

### Scraper Testing (`test_scraper_cloudscraper.py`)

Testing utilities for the AntWiki scraper functionality.

**Usage:**
```bash
# Test scraper functionality
python3 scripts/phenotype/test_scraper_cloudscraper.py
```

**Output Structure:**
```
output/phenotype/
├── antwiki_data/                  # Scraped AntWiki data
│   ├── species_data.json
│   ├── phenotype_summaries.json
│   └── scraping_metadata.json
├── statistical_analysis/          # Statistical analysis results
│   ├── trait_distributions.json
│   ├── descriptive_statistics.json
│   └── normality_tests.json
├── correlation_analysis/          # Correlation analysis results
│   ├── correlation_matrix.json
│   ├── correlation_significance.json
│   └── correlation_plots/
├── data_quality/                  # Quality assessment results
│   ├── missing_data_report.json
│   ├── outlier_detection.json
│   └── data_validation.json
├── phenotype_plots/               # Generated visualizations
│   ├── trait_distributions.png
│   ├── correlation_heatmap.png
│   ├── phenotype_boxplots.png
│   └── statistical_summary.png
└── analysis_report.json           # Comprehensive analysis report
```

## Key Features

✅ **Multi-Source Support**: AntWiki, CSV/TSV, JSON phenotype data
✅ **Comprehensive Scraping**: Automated AntWiki species data extraction
✅ **Statistical Analysis**: Trait distributions, correlations, quality metrics
✅ **Data Validation**: Quality assessment and outlier detection
✅ **Rich Visualization**: Statistical plots and phenotype analysis graphics
✅ **Batch Processing**: Efficient handling of multiple species/samples

## Integration

Integrates with:
- **metainformant.phenotype**: Core phenotype analysis functionality
- **Web scraping**: BeautifulSoup, cloudscraper for AntWiki access
- **Statistical analysis**: pandas, numpy, scipy
- **Data validation**: Custom quality assessment utilities
- **Core utilities**: I/O, logging, configuration management

## Dependencies

- **metainformant.phenotype**: Phenotype analysis module
- **beautifulsoup4**: HTML parsing for web scraping
- **cloudscraper**: Anti-bot protection handling
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Phenotype Analysis Documentation](../../docs/phenotype/README.md)
- [AntWiki Integration](../../docs/phenotype/antwiki.md)
- [Data Curation](../../docs/phenotype/data_curation.md)
- [METAINFORMANT CLI](../../docs/cli.md)
