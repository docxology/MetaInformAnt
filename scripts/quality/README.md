# Quality Control Scripts

Data quality assessment and validation workflow orchestrators.

## Directory Structure

```
scripts/quality/
├── run_quality_control.py         # Quality control workflow orchestrator
└── README.md                      # This file
```

## Quality Control Workflow (`run_quality_control.py`)

Comprehensive quality control workflow orchestrator for biological data validation and assessment.

**Features:**
- FASTQ quality assessment
- Sequence quality metrics
- Contamination detection
- Data integrity validation
- Quality reporting and visualization

**Usage:**
```bash
# FASTQ quality assessment
python3 scripts/quality/run_quality_control.py --fastq reads.fastq.gz --output output/quality/fastq

# Full quality control suite
python3 scripts/quality/run_quality_control.py --fastq reads.fastq.gz --check-contamination --generate-report
```

**Options:**
- `--fastq`: Input FASTQ file for quality assessment
- `--output`: Output directory (defaults to output/quality/)
- `--check-contamination`: Perform contamination detection
- `--generate-report`: Create comprehensive quality report
- `--threads`: Number of threads to use
- `--verbose`: Enable verbose logging

**Output Structure:**
```
output/quality/
├── fastq_assessment/              # FASTQ quality results
│   ├── quality_scores.json
│   ├── base_composition.json
│   └── quality_distributions.json
├── contamination_detection/       # Contamination analysis results
│   ├── contamination_report.json
│   ├── contaminant_sequences.fasta
│   └── contamination_summary.json
├── quality_reports/               # Generated quality reports
│   ├── quality_summary.json
│   ├── quality_plots/
│   └── quality_report.html
└── analysis_report.json           # Comprehensive analysis report
```

## Integration

Integrates with:
- **metainformant.quality**: Core quality control functionality
- **FASTQC/cutadapt**: Quality assessment tools
- **Core utilities**: I/O, logging, path management

## Dependencies

- **metainformant.quality**: Quality control module
- **Biopython**: Sequence handling
- **matplotlib/seaborn**: Visualization support

## Related Documentation

- [Quality Control Documentation](../../docs/quality/README.md)
- [METAINFORMANT CLI](../../docs/cli.md)
