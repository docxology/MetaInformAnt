# Protein Analysis Scripts

Protein sequence and structure analysis workflow orchestrators.

## Directory Structure

```
scripts/protein/
├── run_protein_analysis.py        # Protein analysis workflow orchestrator
└── README.md                      # This file
```

## Protein Analysis Workflow (`run_protein_analysis.py`)

Comprehensive protein analysis workflow orchestrator for sequence analysis, structure prediction, and functional annotation.

**Features:**
- Protein sequence analysis and validation
- Structure prediction and modeling
- Domain and motif identification
- Functional annotation and classification

**Usage:**
```bash
# Protein sequence analysis
python3 scripts/protein/run_protein_analysis.py --input sequences.fasta --output output/protein/basic

# Full analysis with structure prediction
python3 scripts/protein/run_protein_analysis.py --input sequences.fasta --analyze-sequence --predict-structure --find-domains
```

**Options:**
- `--input`: Input protein sequences (FASTA format)
- `--output`: Output directory (defaults to output/protein/)
- `--analyze-sequence`: Perform sequence analysis
- `--predict-structure`: Run structure prediction
- `--find-domains`: Identify protein domains
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/protein/
├── sequence_analysis/             # Sequence analysis results
│   ├── sequence_properties.json
│   ├── physicochemical_properties.json
│   └── sequence_statistics.json
├── structure_prediction/          # Structure prediction results
│   ├── predicted_structures.pdb
│   ├── structure_quality.json
│   └── structural_features.json
├── domain_analysis/               # Domain identification results
│   ├── identified_domains.json
│   ├── domain_annotations.json
│   └── domain_visualization.png
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.protein**: Core protein analysis functionality
- **Biopython**: Sequence and structure handling
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.protein**: Protein analysis module
- **Biopython**: Biological sequence analysis
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Protein Analysis Documentation](../../docs/protein/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)
