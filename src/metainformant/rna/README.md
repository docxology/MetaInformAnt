# RNA Analysis Module

The `rna` module provides comprehensive tools for transcriptomic analysis, focusing on RNA sequencing data processing, quantification, and workflow orchestration. This module integrates with external tools like amalgkit while providing a consistent Python interface.

## Overview

This module handles RNA sequencing workflows from raw data to analyzed results:
- **Amalgkit Integration**: Modular wrapper around the amalgkit CLI toolkit
- **Workflow Management**: Complete pipeline planning and execution
- **Data Processing**: RNA-seq processing steps including quality control, alignment, and quantification
- **Metadata Handling**: Transcriptomic metadata retrieval and curation

## Key Components

### Amalgkit Integration (`amalgkit.py`)
Thin, modular wrapper around the external amalgkit CLI toolkit for transcriptomic meta-analysis.

**Key Features:**
- CLI availability checking and validation
- Command construction and execution
- Parameter validation and error handling
- Modular step execution (metadata, integration, configuration, etc.)

**Usage:**
```python
from metainformant.rna import check_cli_available, metadata, quant

# Check if amalgkit is available
available, help_text = check_cli_available()
if not available:
    print(f"Amalgkit not available: {help_text}")

# Run metadata retrieval
result = metadata({"threads": 4, "species": "Homo sapiens"})
print(f"Exit code: {result.returncode}")

# Run quantification
quant_result = quant({"input": "aligned.bam", "output": "counts.txt"})
```

### Workflow Management (`workflow.py`)
Complete workflow planning and execution for complex RNA-seq pipelines.

**Key Features:**
- Workflow configuration and validation
- Step planning and dependency resolution
- Execution orchestration with error handling
- Progress tracking and logging

**Usage:**
```python
from metainformant.rna import AmalgkitWorkflowConfig, plan_workflow, execute_workflow

# Configure workflow
config = AmalgkitWorkflowConfig(
    work_dir="/path/to/work",
    threads=8,
    species_list=["Homo sapiens", "Mus musculus"]
)

# Plan workflow steps
steps = plan_workflow(config)

# Execute complete workflow
results = execute_workflow(config)
for step, result in results.items():
    print(f"{step}: {result.returncode}")
```

## Supported Workflow Steps

The module supports all major amalgkit workflow steps:

### Metadata (`metadata`)
Retrieve and organize transcriptomic metadata from public databases.

**Parameters:**
- `threads`: Number of parallel threads
- `species`: Target species for metadata retrieval
- `output_dir`: Directory for metadata files

### Integration (`integrate`)
Integrate multiple data sources and annotations.

**Parameters:**
- `input_dirs`: Directories containing input data
- `annotation`: Reference annotation file
- `threads`: Processing threads

### Configuration (`config`)
Generate configuration files for downstream analysis.

**Parameters:**
- `species`: Target organism
- `assay`: Assay type (RNA-seq, etc.)
- `stranded`: Library strandedness

### Selection (`select`)
Select and filter relevant features for analysis.

**Parameters:**
- `input`: Input feature file
- `criteria`: Selection criteria
- `output`: Output file

### Data Retrieval (`getfastq`)
Download FASTQ files from public repositories.

**Parameters:**
- `accessions`: SRA/ENA accession numbers
- `output_dir`: Download directory
- `threads`: Download threads

### Quantification (`quant`)
Quantify gene/transcript expression levels.

**Parameters:**
- `input`: Aligned reads (BAM/SAM)
- `annotation`: Transcript annotation
- `method`: Quantification method (featureCounts, etc.)

### Merging (`merge`)
Merge results from parallel processing.

**Parameters:**
- `input_files`: Files to merge
- `output`: Merged output file
- `method`: Merge strategy

### Statistical Testing (`cstmm`, `csca`)
Statistical testing for differential expression.

**Parameters:**
- `counts`: Count matrix
- `design`: Experimental design matrix
- `method`: Statistical test method

### Curation (`curate`)
Data curation and quality assessment.

**Parameters:**
- `input`: Raw data directory
- `criteria`: Curation criteria
- `output`: Curated results

### Sanity Checking (`sanity`)
Validate data integrity and completeness.

**Parameters:**
- `input_dirs`: Directories to validate
- `checks`: Validation checks to perform

## Command Line Interface

The module provides a comprehensive CLI for RNA analysis:

```bash
# Check amalgkit availability
python -m metainformant rna check

# Run complete workflow
python -m metainformant rna run --work-dir /path/to/work --threads 8 --species "Homo sapiens"

# Run individual steps
python -m metainformant rna metadata --threads 4 --species "Homo sapiens"
python -m metainformant rna quant --input aligned.bam --output counts.txt

# Plan workflow without execution
python -m metainformant rna plan --work-dir /path/to/work --species "Homo sapiens"
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences
from metainformant.rna import workflow

# Use DNA sequences for RNA annotation
gene_sequences = sequences.read_fasta("genes.fasta")
config = AmalgkitWorkflowConfig(annotation_fasta=gene_sequences)
```

### With Visualization Module
```python
from metainformant.rna import workflow
from metainformant.visualization import heatmap, lineplot

# Visualize expression data
results = workflow.execute_workflow(config)
heatmap(results["expression_matrix"])
```

### With Statistical Analysis
```python
from metainformant.rna import cstmm

# Statistical analysis of expression data
# cstmm returns results with p-values directly
test_results = cstmm({"counts": counts_matrix, "design": design_matrix})
# Extract p-values from test_results as needed
```

## Configuration Management

The module supports comprehensive configuration management:

```python
from metainformant.core import config as core_config
from metainformant.rna import AmalgkitWorkflowConfig

# Load configuration from file
cfg = core_config.load_config("rna_config.yaml")

# Create workflow configuration
workflow_cfg = AmalgkitWorkflowConfig(
    work_dir=cfg.work_dir,
    threads=cfg.compute.threads,
    species_list=cfg.species,
    **cfg.amalgkit_params
)
```

## Error Handling and Validation

- **Defensive Imports**: Optional dependencies are imported defensively
- **CLI Validation**: Amalgkit availability is checked before execution
- **Parameter Validation**: All parameters are validated before execution
- **Error Propagation**: Clear error messages with context

## Performance Features

- **Parallel Execution**: Multi-threaded processing for large datasets
- **Streaming Processing**: Memory-efficient handling of large files
- **Progress Tracking**: Real-time progress reporting for long-running operations
- **Resource Management**: Proper cleanup and resource management

## Testing

The module includes comprehensive tests:
- CLI availability and functionality tests
- Workflow execution validation
- Parameter validation tests
- Integration tests with real data
- Performance benchmarks

## Dependencies

- **Required**: Python 3.11+ for full workflow functionality
- **Optional**: amalgkit CLI (external dependency)
- **Core**: metainformant.core for configuration and I/O

## Usage Examples

### Basic CLI Usage
```python
from metainformant.rna import check_cli_available, run_amalgkit

# Check availability
if check_cli_available()[0]:
    # Run metadata step
    result = run_amalgkit("metadata", {"threads": 4})
    print(f"Metadata completed with code: {result.returncode}")
```

### Complete Workflow
```python
from metainformant.rna import AmalgkitWorkflowConfig, execute_workflow

# Setup and run complete workflow
config = AmalgkitWorkflowConfig(
    work_dir="/tmp/rna_analysis",
    threads=8,
    species_list=["Homo sapiens"],
    steps=["metadata", "getfastq", "quant", "cstmm"]
)

results = execute_workflow(config)
print("Workflow completed successfully!")
```

This module provides a complete solution for transcriptomic analysis, from individual processing steps to complete workflow orchestration.
