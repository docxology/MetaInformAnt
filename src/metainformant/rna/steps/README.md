# RNA Workflow Steps

This directory contains individual step implementations for the METAINFORMANT RNA analysis workflow, providing modular components that can be composed into complete analysis pipelines.

## Step Modules

Each module implements a specific step in the RNA-seq analysis workflow:

- **`amalgkit_step_*.py`**: Individual amalgkit workflow steps
- **`metadata.py`**: Metadata retrieval and organization
- **`integration.py`**: Data source integration
- **`config.py`**: Configuration generation
- **`select.py`**: Feature selection and filtering
- **`getfastq.py`**: FASTQ file retrieval
- **`quant.py`**: Expression quantification
- **`merge.py`**: Result merging and aggregation
- **`cstmm.py`**: Statistical testing (CSTMM)
- **`csca.py`**: Statistical testing (CSCA)
- **`curate.py`**: Data curation and quality control
- **`sanity.py`**: Sanity checking and validation

## Architecture

### Modular Design
Each step is designed as an independent module that:
- Accepts standardized input parameters
- Performs specific analysis operations
- Returns structured results
- Can be chained with other steps

### Parameter Standardization
All steps follow consistent parameter patterns:
```python
def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute a single workflow step."""
    # Step implementation
    return StepResult(success=True, outputs=outputs, metadata=metadata)
```

## Integration

### With Main RNA Module
```python
from metainformant.rna import workflow
from metainformant.rna.steps import metadata, quant

# Use individual steps in custom workflows
meta_result = metadata.run_step({"species": "Homo sapiens", "threads": 4})
quant_result = quant.run_step({"input": meta_result.outputs["bam_file"]})
```

### With Workflow Orchestration
```python
from metainformant.rna.workflow import plan_workflow, execute_workflow

# Steps are automatically orchestrated in complete workflows
config = AmalgkitWorkflowConfig(steps=["metadata", "getfastq", "quant"])
planned_steps = plan_workflow(config)
results = execute_workflow(config)
```

## Testing

Each step module includes comprehensive tests:
- Input validation and parameter checking
- Output format verification
- Error handling and edge cases
- Integration testing with real data

## Contributing

### Adding New Steps
1. Follow the established parameter interface
2. Include comprehensive error handling
3. Add unit and integration tests
4. Update workflow planning logic

### Step Interface Requirements
- Accept dictionary of parameters
- Return standardized StepResult objects
- Handle errors gracefully with informative messages
- Support both CLI and programmatic usage

## Related Documentation

- See `src/metainformant/rna/README.md` for complete workflow documentation
- See `src/metainformant/rna/workflow.py` for orchestration details
- See `docs/rna/steps.md` for detailed step-by-step documentation
