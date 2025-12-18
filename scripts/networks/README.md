# Network Analysis Scripts

Biological network analysis and graph theory workflow orchestrators.

## Directory Structure

```
scripts/networks/
├── run_network_analysis.py        # Network analysis workflow orchestrator
└── README.md                      # This file
```

## Network Analysis Workflow (`run_network_analysis.py`)

Comprehensive network analysis workflow orchestrator for biological networks, including construction, metrics calculation, community detection, and visualization.

**Features:**
- Network construction from interaction data
- Graph theory metrics and analysis
- Community detection algorithms
- Centrality analysis
- Network visualization and plotting

**Usage:**
```bash
# Basic network construction and metrics
python3 scripts/networks/run_network_analysis.py --input interactions.tsv --output output/networks/basic

# Full analysis with all modules
python3 scripts/networks/run_network_analysis.py --input interactions.tsv --analyze-metrics --detect-communities --analyze-centrality

# Metrics and centrality only
python3 scripts/networks/run_network_analysis.py --input interactions.tsv --analyze-metrics --analyze-centrality
```

**Options:**
- `--input`: Input interaction data (TSV/CSV format: node1, node2, [weight])
- `--output`: Output directory (defaults to output/networks/)
- `--analyze-metrics`: Calculate network topology metrics
- `--detect-communities`: Perform community detection
- `--analyze-centrality`: Compute centrality measures
- `--visualize`: Generate network visualizations
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/networks/
├── network_data/                  # Constructed network files
│   ├── network_graph.pkl
│   ├── adjacency_matrix.csv
│   └── node_attributes.json
├── metrics/                       # Network metrics results
│   ├── global_metrics.json
│   ├── local_metrics.json
│   └── degree_distribution.json
├── communities/                   # Community detection results
│   ├── community_membership.json
│   ├── modularity_scores.json
│   └── community_summary.json
├── centrality/                    # Centrality analysis results
│   ├── degree_centrality.json
│   ├── betweenness_centrality.json
│   ├── closeness_centrality.json
│   └── eigenvector_centrality.json
├── visualizations/                # Generated network plots
│   ├── network_graph.png
│   ├── degree_distribution.png
│   ├── community_structure.png
│   └── centrality_heatmap.png
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.networks**: Core network analysis functionality
- **NetworkX**: Graph algorithms and analysis
- **Core utilities**: I/O, logging, path management
- **Visualization**: Network plotting and graph visualization

## Dependencies

- **metainformant.networks**: Network analysis module
- **NetworkX**: Graph theory algorithms
- **pandas/numpy**: Data manipulation and statistics
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Network Analysis Documentation](../../docs/networks/README.md)
- [Graph Theory Methods](../../docs/networks/graph_theory.md)
- [Community Detection](../../docs/networks/community_detection.md)
- [METAINFORMANT CLI](../../docs/cli.md)

