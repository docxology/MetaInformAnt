# Visualization Scripts

Biological data visualization and plotting workflow orchestrators.

## Directory Structure

```
scripts/visualization/
├── run_visualization.py           # Visualization workflow orchestrator
└── README.md                      # This file
```

## Visualization Workflow (`run_visualization.py`)

Comprehensive visualization workflow orchestrator for generating publication-quality plots and interactive visualizations.

**Features:**
- Line plots and time series
- Heatmaps and correlation plots
- Network and graph visualizations
- Statistical plot generation
- Interactive dashboard creation
- Customizable styling and themes

**Usage:**
```bash
# Generate line plots
python3 scripts/visualization/run_visualization.py --data time_series.csv --plot-type line --output output/visualization/plots

# Create heatmap
python3 scripts/visualization/run_visualization.py --data correlation_matrix.csv --plot-type heatmap --output output/visualization/heatmaps

# Generate multiple plot types
python3 scripts/visualization/run_visualization.py --data analysis_results.json --all-plots --output output/visualization/complete
```

**Options:**
- `--data`: Input data file for visualization
- `--plot-type`: Type of plot to generate (line, heatmap, network, etc.)
- `--output`: Output directory (defaults to output/visualization/)
- `--all-plots`: Generate all applicable plot types
- `--interactive`: Create interactive visualizations
- `--style`: Plot styling theme
- `--dpi`: Resolution for static plots
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/visualization/
├── line_plots/                    # Line and time series plots
│   ├── time_series.png
│   ├── trend_analysis.png
│   └── line_plot_data.json
├── heatmaps/                      # Heatmap visualizations
│   ├── correlation_heatmap.png
│   ├── expression_heatmap.png
│   └── heatmap_data.json
├── network_plots/                 # Network and graph visualizations
│   ├── interaction_network.png
│   ├── pathway_diagram.png
│   └── network_data.json
├── statistical_plots/             # Statistical visualizations
│   ├── box_plots.png
│   ├── violin_plots.png
│   ├── histogram_plots.png
│   └── statistical_data.json
├── interactive/                   # Interactive visualizations
│   ├── dashboard.html
│   ├── interactive_plots.html
│   └── interactive_data.json
└── visualization_report.json      # Comprehensive visualization report
```

## Integration

Integrates with:
- **metainformant.visualization**: Core visualization functionality
- **matplotlib/seaborn**: Static plotting
- **plotly/bokeh**: Interactive visualizations
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.visualization**: Visualization module
- **matplotlib**: Core plotting
- **seaborn**: Statistical visualization
- **plotly**: Interactive plots
- **bokeh**: Web-based visualizations

## Related Documentation

- [Visualization Documentation](../../docs/visualization/README.md)
- [Plot Types](../../docs/visualization/plot_types.md)
- [METAINFORMANT CLI](../../docs/cli.md)

