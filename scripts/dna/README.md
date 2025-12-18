# DNA Analysis Scripts

DNA sequence analysis and genomics workflow orchestrators.

## Directory Structure

```
scripts/dna/
├── run_dna_analysis.py          # DNA analysis workflow orchestrator
└── README.md                    # This file
```

## DNA Analysis Workflow (`run_dna_analysis.py`)

Comprehensive DNA analysis workflow orchestrator for sequence processing, quality control, and genomic analysis.

**Features:**
- Sequence quality control and filtering
- Composition analysis (GC content, nucleotide frequencies)
- Population genetics analysis
- Phylogenetic analysis support
- Integration with metainformant.dna module

**Usage:**
```bash
# Basic sequence analysis
python3 scripts/dna/run_dna_analysis.py --input sequences.fasta --output output/dna/basic

# Full analysis with all modules
python3 scripts/dna/run_dna_analysis.py --input sequences.fasta --analyze-composition --analyze-population --analyze-phylogeny

# Quality control and composition only
python3 scripts/dna/run_dna_analysis.py --input sequences.fasta --analyze-composition --min-length 100
```

**Options:**
- `--input`: Input DNA sequences (FASTA format)
- `--output`: Output directory (defaults to output/dna/)
- `--analyze-composition`: Perform sequence composition analysis
- `--analyze-population`: Perform population genetics analysis
- `--min-length`: Minimum sequence length filter
- `--max-length`: Maximum sequence length filter
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/dna/
├── sequences_summary.json        # Sequence statistics
├── composition_analysis.json     # GC content, nucleotide frequencies
├── population_analysis.json      # Population genetics metrics
├── filtered_sequences.fasta      # Quality-filtered sequences
└── analysis_report.json          # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.dna**: Core DNA analysis functionality
- **Core utilities**: I/O, logging, path management
- **Quality control**: Sequence filtering and validation
- **Visualization**: Plot generation for analysis results

## Dependencies

- **metainformant.dna**: DNA sequence analysis module
- **BioPython**: Sequence parsing and manipulation
- **NumPy/SciPy**: Statistical analysis
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [DNA Analysis Documentation](../../docs/dna/README.md)
- [Core Utilities](../../docs/core/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)

