# Epigenome Analysis Scripts

Epigenetic modification analysis and chromatin biology workflow orchestrators.

## Directory Structure

```
scripts/epigenome/
├── run_epigenome_analysis.py     # Epigenome analysis workflow orchestrator
└── README.md                     # This file
```

## Epigenome Analysis Workflow (`run_epigenome_analysis.py`)

Comprehensive epigenome analysis workflow orchestrator for methylation analysis, chromatin track processing, and epigenetic data analysis.

**Features:**
- DNA methylation analysis (CpG beta values, M-values)
- Chromatin track processing (BedGraph, BigWig)
- Genome-wide methylation summaries
- Chromosomal distribution analysis
- Quality control and validation

**Usage:**
```bash
# Methylation analysis
python3 scripts/epigenome/run_epigenome_analysis.py --methylation cpg_table.tsv --compute-beta --summarize-by-chromosome

# Chromatin track analysis
python3 scripts/epigenome/run_epigenome_analysis.py --bedgraph track.bedgraph --output output/epigenome/tracks

# Full methylation pipeline
python3 scripts/epigenome/run_epigenome_analysis.py --methylation cpg.tsv --compute-beta --summarize-by-chromosome --output output/epigenome/methylation
```

**Options:**
- `--methylation`: Input CpG methylation table (TSV format)
- `--bedgraph`: Input chromatin track (BedGraph format)
- `--bigwig`: Input chromatin track (BigWig format)
- `--compute-beta`: Compute beta values from methylation ratios
- `--summarize-by-chromosome`: Generate chromosome-wise summaries
- `--output`: Output directory (defaults to output/epigenome/)
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/epigenome/
├── methylation_summary.json      # Genome-wide methylation statistics
├── beta_values.json             # Computed beta values
├── chromosomal_summary.json     # Per-chromosome methylation profiles
├── chromatin_tracks.json        # Processed chromatin track data
├── qc_metrics.json              # Quality control metrics
├── methylation_plots/           # Generated visualizations
│   ├── beta_value_distribution.png
│   ├── chromosomal_methylation.png
│   └── chromatin_tracks.png
└── analysis_report.json         # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.epigenome**: Core epigenetic analysis functionality
- **Core utilities**: I/O, logging, path management
- **Genomic coordinates**: Integration with DNA sequence data
- **Visualization**: Plot generation for epigenetic data

## Dependencies

- **metainformant.epigenome**: Epigenetic analysis module
- **pyBigWig**: BigWig file handling
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Epigenome Analysis Documentation](../../docs/epigenome/README.md)
- [Core Utilities](../../docs/core/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)
