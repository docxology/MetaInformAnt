# AI Agents in Amalgkit Workflow Script Development

This document outlines AI assistance in developing and maintaining METAINFORMANT's amalgkit workflow verification and monitoring scripts.

## AI Contributions

### Script Development
**Code Assistant Agent** implemented:
- **`verify_workflow.sh`**: Comprehensive workflow validation and status checking
  - Multi-species workflow execution verification
  - Step completion and output file validation
  - Error detection and diagnostic reporting
  - Status summary generation across all species
- **Workflow Monitoring**: Real-time progress tracking and logging analysis
- **Error Handling**: Robust error detection and recovery patterns
- **Documentation Integration**: Cross-referencing with workflow documentation

### Automation Features
**Code Assistant Agent** designed:
- **Batch Processing**: Multi-species workflow orchestration
- **Progress Tracking**: Real-time workflow execution monitoring
- **File Validation**: Comprehensive output file existence and integrity checks
- **Status Reporting**: Clear, actionable workflow status summaries
- **Diagnostic Tools**: Error identification and troubleshooting assistance

### Integration Patterns
**Documentation Agent** contributed to:
- Script usage documentation and examples
- Integration with METAINFORMANT's Python workflow API
- Cross-reference with amalgkit step documentation
- Troubleshooting guides for common workflow issues

## Script Architecture

### Workflow Verification (`verify_workflow.sh`)

**Purpose**: Validate amalgkit workflow execution and output completeness

**Key Features**:
- ✅ **Multi-Species Support**: Verifies workflows for multiple species simultaneously
- ✅ **Step Validation**: Checks completion of all 11 amalgkit workflow steps
- ✅ **Output Verification**: Validates existence of expected output files
- ✅ **Error Detection**: Identifies failed steps and missing outputs
- ✅ **Status Reporting**: Generates comprehensive workflow status summaries
- ✅ **Diagnostic Mode**: Provides detailed debugging information

**Validation Checks**:
```bash
# Per-species validation includes:
- Genome download completion
- Metadata retrieval and filtering
- Configuration file generation
- Sample selection and prioritization
- FASTQ file generation
- Transcript quantification
- Expression matrix merging
- Cross-species normalization (if applicable)
- Quality control and curation
- Sanity check execution
```

**AI Contributions**:
- Systematic validation logic for all workflow steps
- Error pattern recognition from production runs
- Diagnostic message generation for common issues
- Multi-species execution coordination

## Development Approach

### Modular Design
AI helped establish:
- **Function Decomposition**: Separate validation logic for each workflow step
- **Reusable Components**: Common validation patterns across steps
- **Error Handling**: Consistent error detection and reporting
- **Extensibility**: Easy addition of new validation checks

### Production Validation
AI analyzed:
- **Real Workflow Executions**: Patterns from actual multi-species runs
- **Common Failures**: Typical error scenarios and their indicators
- **File Patterns**: Expected output file structures and naming
- **Performance Characteristics**: Normal execution times and resource usage

### User Experience
AI optimized:
- **Clear Output**: Human-readable status messages and summaries
- **Actionable Diagnostics**: Specific troubleshooting recommendations
- **Progress Indicators**: Visual feedback during validation
- **Documentation Links**: References to relevant troubleshooting docs

## Workflow Validation Logic

### Step-by-Step Verification

#### 1. **Genome Validation**
```bash
# Checks for:
- genome/ncbi_dataset_api.zip existence
- genome/download.heartbeat timestamp
- Extracted genome files and annotations
```

#### 2. **Metadata Validation**
```bash
# Verifies:
- work/metadata/metadata_original.tsv
- Filtered metadata files by criteria
- Metadata filtering logs
```

#### 3. **Configuration Validation**
```bash
# Confirms:
- work/config_base/*.config files
- Configuration file completeness
- Parameter consistency
```

#### 4. **Selection Validation**
```bash
# Validates:
- work/pivot_qualified.tsv
- Selected sample lists
- Selection criteria application
```

#### 5. **FASTQ Validation**
```bash
# Checks:
- fastq/*.fastq.gz files
- FASTQ download logs
- File integrity and completeness
```

#### 6. **Quantification Validation**
```bash
# Verifies:
- quant/*/abundance.tsv files
- Kallisto index existence
- Quantification logs and metrics
```

#### 7. **Merge Validation**
```bash
# Confirms:
- merged/merge/merge.tpm.tsv
- merged/merge/merge.count.tsv
- Expression matrix integrity
```

#### 8. **Curation Validation**
```bash
# Validates:
- curate/*.uncorrected.*.tsv
- curate/*.corrected.*.tsv
- Quality control plots and reports
```

#### 9. **Sanity Check Validation**
```bash
# Checks:
- work/sanity/ validation reports
- Integrity check results
- Error and warning summaries
```

### Cross-Species Validation

For multi-species workflows:
```bash
# Additional checks:
- Ortholog table existence and format
- Cross-species normalization outputs
- Comparative analysis results
- Cross-species correlation plots
```

## Usage Examples

### Basic Workflow Verification
```bash
# Verify single species workflow
cd /path/to/metainformant
bash scripts/rna/amalgkit/verify_workflow.sh amellifera

# Verify all species
bash scripts/rna/amalgkit/verify_workflow.sh all
```

### Diagnostic Mode
```bash
# Detailed diagnostic output
bash scripts/rna/amalgkit/verify_workflow.sh pbarbatus --verbose

# Generate validation report
bash scripts/rna/amalgkit/verify_workflow.sh --report
```

### Integration with Python API
```python
from metainformant.rna.workflow import validate_workflow_outputs
from pathlib import Path

# Validate workflow programmatically
results = validate_workflow_outputs(
    work_dir=Path("output/amalgkit/pbarbatus"),
    check_all_steps=True
)

print(f"Validation: {'PASSED' if results.all_passed else 'FAILED'}")
for step, status in results.step_status.items():
    print(f"  {step}: {status}")
```

## Best Practices

### Regular Validation
- **After Each Step**: Validate completion before proceeding
- **Post-Workflow**: Comprehensive validation after full workflow
- **Pre-Analysis**: Verify data integrity before downstream analysis
- **Troubleshooting**: Use diagnostic mode to identify issues

### Production Workflows
- **Automated Checks**: Integrate validation into workflow automation
- **Failure Recovery**: Use diagnostics to identify restart points
- **Documentation**: Log validation results for reproducibility
- **Version Tracking**: Validate against expected output formats

## Troubleshooting Integration

### Common Issues Detected
AI-identified common workflow issues:
- Missing FASTQ files → Check getfastq logs and SRA availability
- Incomplete quantification → Verify kallisto index and FASTQ quality
- Failed merging → Check quantification output consistency
- Sanity check failures → Review input data integrity

### Diagnostic Recommendations
For each detected issue, the script provides:
- **Specific Error**: Clear description of the validation failure
- **Likely Cause**: Common reasons for the specific failure
- **Resolution Steps**: Actionable troubleshooting recommendations
- **Documentation Links**: References to relevant troubleshooting guides

## Maintenance and Updates

### Script Evolution
- **Workflow Changes**: Updates with amalgkit version changes
- **New Validations**: Additional checks based on production experience
- **Performance Improvements**: Optimization of validation speed
- **Documentation Sync**: Keep aligned with workflow documentation

### Community Contributions
- **Issue Reports**: Incorporate validation checks for reported issues
- **Best Practices**: Add validation patterns from community feedback
- **Feature Requests**: Implement requested diagnostic capabilities

## Related Documentation

This script documentation integrates with:
- **[README.md](README.md)**: Amalgkit script overview and usage
- **[docs/rna/amalgkit/comprehensive_guide.md](../../../docs/rna/amalgkit/comprehensive_guide.md)**: Complete workflow guide
- **[docs/rna/amalgkit/steps/](../../../docs/rna/amalgkit/steps/)**: Individual step documentation
- **[docs/rna/workflow.md](../../../docs/rna/workflow.md)**: Workflow orchestration
- **[CENTRALIZATION_COMPLETE.md](CENTRALIZATION_COMPLETE.md)**: Workflow consolidation status

## Quality Assurance

### Validation Testing
All validation scripts are tested against:
- ✅ **Successful Workflows**: Verify correct validation passing
- ✅ **Failed Workflows**: Ensure accurate failure detection
- ✅ **Partial Workflows**: Handle incomplete executions gracefully
- ✅ **Multi-Species**: Validate complex comparative analyses

### Production Validation
Scripts validated through:
- Real multi-species workflow executions (amellifera, cfloridanus, mpharaonis, pbarbatus, sinvicta)
- Large-scale data processing (100+ samples per species)
- Various failure scenarios and recovery patterns
- Performance testing on different system configurations

---

*These workflow verification scripts demonstrate effective AI-assisted development of production-ready bioinformatics workflow tooling, combining systematic validation logic with practical troubleshooting guidance.*

**Scripts Updated**: October 29, 2025  
**METAINFORMANT Version**: 1.0  
**Amalgkit Compatibility**: 0.12.19  
**Status**: ✅ Production-tested, comprehensive validation


