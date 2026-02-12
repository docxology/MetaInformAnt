# RNA Workflow Sample Validation Guide

Complete guide to understanding and using the RNA-seq workflow sample validation system.

## Overview

The validation system tracks samples through the complete pipeline stages:
- **Download**: SRA files downloaded
- **Extraction**: FASTQ files extracted
- **Quantification**: Abundance files generated
- **Merge**: Sample included in merged abundance matrix

## Quick Start

### Automatic Validation

Validation runs automatically after key workflow steps:

```bash
# Validation runs automatically after getfastq and quant steps
python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml
```

Validation reports are saved to `work_dir/validation/`:
- `getfastq_validation.json` - Download and extraction validation
- `quant_validation.json` - Quantification validation

### Standalone Validation

Run validation without executing the workflow:

```bash
# Validate all stages
python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --validate

# Validate specific stage
python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --validate --validate-stage quantification
```

### Programmatic Validation

```python
from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.validation import validate_all_samples, get_sample_pipeline_status

# Load config
config = load_workflow_config("config/amalgkit/my_species.yaml")

# Validate all samples
result = validate_all_samples(config)
print(f"Validated: {result['validated']}/{result['total_samples']}")

# Check specific sample
status = get_sample_pipeline_status("SRR123456", config.work_dir)
print(f"Stage: {status['stage']}")
```

## Validation Report Schema

### Report Structure

```json
{
  "total_samples": 100,
  "validated": 95,
  "failed": 5,
  "validation_stage": "all",
  "missing_stages": {
    "download": 0,
    "extraction": 2,
    "quantification": 3,
    "merge": 5
  },
  "per_sample": {
    "SRR123456": {
      "valid": true,
      "sample_id": "SRR123456",
      "stages": {
        "download": true,
        "extraction": true,
        "quantification": true,
        "merge": false
      },
      "missing_stages": ["merge"],
      "current_stage": "quantification",
      "diagnostics": {
        "sra_file": "path/to/SRR123456.sra",
        "sra_size": 12345678,
        "fastq_files": ["path/to/SRR123456_1.fastq.gz", "path/to/SRR123456_2.fastq.gz"],
        "fastq_count": 2,
        "fastq_total_size": 98765432,
        "abundance_file": "path/to/abundance.tsv",
        "abundance_size": 12345
      }
    }
  },
  "summary": {
    "download": {
      "total": 100,
      "complete": 100,
      "missing": 0
    },
    "extraction": {
      "total": 100,
      "complete": 98,
      "missing": 2
    },
    "quantification": {
      "total": 100,
      "complete": 97,
      "missing": 3
    },
    "merge": {
      "total": 100,
      "complete": 95,
      "missing": 5
    }
  }
}
```

### Field Descriptions

- **total_samples**: Total number of samples in metadata
- **validated**: Number of samples that passed validation (extraction + quantification)
- **failed**: Number of samples that failed validation
- **validation_stage**: Stage being validated ("all", "download", "extraction", "quantification", "merge")
- **missing_stages**: Count of samples missing each stage
- **per_sample**: Detailed status for each sample
  - **valid**: Whether sample passed validation
  - **stages**: Boolean status for each pipeline stage
  - **missing_stages**: List of incomplete stages
  - **current_stage**: Current pipeline stage
  - **diagnostics**: File paths, sizes, and detailed information
- **summary**: Stage-specific summary statistics

## Troubleshooting

### Common Issues

#### Issue: "Metadata file not found"

**Symptoms**:
- Validation returns `total_samples: 0`
- Error message: "Metadata file not found"

**Solutions**:
1. Run `metadata` step to generate metadata file:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --steps metadata
   ```
2. Run `select` step to create selected metadata:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --steps select
   ```
3. Check that metadata files exist:
   ```bash
   ls output/amalgkit/my_species/work/metadata/
   ```

#### Issue: "No samples found in metadata"

**Symptoms**:
- Validation returns `total_samples: 0`
- Metadata file exists but has no sample IDs

**Solutions**:
1. Check metadata file format:
   ```bash
   head output/amalgkit/my_species/work/metadata/metadata.tsv
   ```
2. Ensure metadata has a column named `run`, `Run`, `SRA_Run`, `sample_id`, or `accession`
3. Verify samples were selected in `select` step

#### Issue: Samples Missing FASTQ Files

**Symptoms**:
- Validation shows `extraction: false` for samples
- `missing_stages` includes "extraction"

**Solutions**:
1. Check if `getfastq` step completed:
   ```bash
   ls output/amalgkit/my_species/fastq/getfastq/
   ```
2. Re-run `getfastq` step for failed samples:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --steps getfastq
   ```
3. Check for download errors in logs:
   ```bash
   cat output/amalgkit/my_species/work/logs/getfastq.log
   ```

#### Issue: Samples Missing Quantification Files

**Symptoms**:
- Validation shows `quantification: false` for samples
- `missing_stages` includes "quantification"

**Solutions**:
1. Check if `quant` step completed:
   ```bash
   ls output/amalgkit/my_species/quant/
   ```
2. Verify FASTQ files exist for samples:
   ```bash
   ls output/amalgkit/my_species/fastq/getfastq/SRR123456/
   ```
3. Re-run `quant` step:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --steps quant
   ```

#### Issue: Samples Not in Merged Matrix

**Symptoms**:
- Validation shows `merge: false` for samples
- Samples quantified but not in merged abundance file

**Solutions**:
1. Run `merge` step:
   ```bash
   python3 scripts/rna/run_workflow.py config/amalgkit/my_species.yaml --steps merge
   ```
2. Check merged abundance file:
   ```bash
   head output/amalgkit/my_species/merged/merged_abundance.tsv
   ```

### Interpreting Validation Results

#### All Samples Valid

```json
{
  "total_samples": 100,
  "validated": 100,
  "failed": 0
}
```

**Meaning**: All samples completed extraction and quantification stages.

#### Some Samples Failed

```json
{
  "total_samples": 100,
  "validated": 95,
  "failed": 5,
  "missing_stages": {
    "extraction": 2,
    "quantification": 3
  }
}
```

**Meaning**: 
- 95 samples are valid
- 2 samples missing FASTQ files (extraction stage)
- 3 samples missing abundance files (quantification stage)

**Action**: Check `per_sample` details for specific sample IDs and remediation steps.

#### Stage-Specific Validation

When validating a specific stage:

```bash
--validate --validate-stage extraction
```

Returns validation for that stage only:
- `validated`: Samples that passed that stage
- `failed`: Samples that failed that stage
- `missing_stages`: Only includes the specified stage

## Performance

### Expected Validation Time

- **Small datasets** (< 100 samples): < 1 second
- **Medium datasets** (100-1000 samples): 1-5 seconds
- **Large datasets** (> 1000 samples): 5-30 seconds

### Optimization

Validation is optimized for performance:
- **Streaming file reading**: Merge validation reads files in chunks
- **Early exit**: Stops searching when sample found
- **Progress indication**: Logs progress for large sample sets (> 10 samples)

### Impact on Workflow

- **Automatic validation**: Adds < 1% overhead to workflow execution
- **Standalone validation**: Runs independently, no workflow impact
- **Large files**: Merge validation limits search to first 1000 rows for very large files

## Examples

### Reading Validation Reports

```python
from metainformant.core.io import load_json
from pathlib import Path

# Load validation report
report_path = Path("output/amalgkit/my_species/work/validation/quant_validation.json")
report = load_json(report_path)

# Check overall status
print(f"Validated: {report['validated']}/{report['total_samples']}")

# Find failed samples
failed_samples = [
    sample_id for sample_id, status in report['per_sample'].items()
    if not status['valid']
]
print(f"Failed samples: {failed_samples}")

# Get diagnostics for a specific sample
sample_status = report['per_sample']['SRR123456']
print(f"Missing stages: {sample_status['missing_stages']}")
if 'quantification_issue' in sample_status['diagnostics']:
    print(sample_status['diagnostics']['quantification_issue'])
```

### Custom Validation Logic

```python
from metainformant.rna.validation import get_sample_pipeline_status

def check_sample_ready_for_analysis(sample_id: str, work_dir: Path) -> bool:
    """Check if sample is ready for downstream analysis."""
    status = get_sample_pipeline_status(sample_id, work_dir)
    
    # Require at least quantification
    return status['quantification'] is True

# Use in workflow
samples = ["SRR123456", "SRR789012"]
ready_samples = [
    s for s in samples 
    if check_sample_ready_for_analysis(s, work_dir)
]
```

### Batch Validation Script

```python
#!/usr/bin/env python3
"""Validate multiple species workflows."""

from pathlib import Path
from metainformant.rna.workflow import load_workflow_config
from metainformant.rna.validation import validate_all_samples

species_configs = [
    "config/amalgkit/amalgkit_pbarbatus.yaml",
    "config/amalgkit/amalgkit_camponotus_floridanus.yaml",
]

for config_path in species_configs:
    print(f"\nValidating {config_path}...")
    config = load_workflow_config(config_path)
    result = validate_all_samples(config)
    
    print(f"  Total: {result['total_samples']}")
    print(f"  Validated: {result['validated']}")
    print(f"  Failed: {result['failed']}")
    
    if result['failed'] > 0:
        print(f"  Missing stages: {result['missing_stages']}")
```

## Best Practices

1. **Run validation after each major step**: Use automatic validation or run `--validate` after `getfastq` and `quant`
2. **Check validation reports**: Review JSON reports for detailed diagnostics
3. **Address failures promptly**: Fix failed samples before proceeding to next steps
4. **Use stage-specific validation**: Validate specific stages during troubleshooting
5. **Monitor progress**: Watch for validation progress logs in large workflows

## Related Documentation

- **[getfastq Step Guide](amalgkit/steps/04_getfastq.md)**: Download and extraction details
- **[Workflow Guide](workflow.md)**: Complete workflow execution
- **[Configuration Guide](CONFIGURATION.md)**: Workflow configuration
- **[README](../src/metainformant/rna/README.md)**: Module overview


