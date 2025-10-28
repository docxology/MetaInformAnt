# Configuration Directory

This directory contains repository-level configuration files used by METAINFORMANT for various tools and workflows, with a focus on RNA-seq analysis configurations.

## RNA-seq Analysis Configurations

### Amalgkit Workflow Configurations

#### Core Pbarbatus Dataset Configurations
- **`amalgkit_pbarbatus.yaml`**: Standard configuration for amalgkit RNA-seq analysis with Pbarbatus dataset (Apis barbatus)
- **`amalgkit_pbarbatus_direct.yaml`**: Direct execution configuration for streamlined Pbarbatus analysis
- **`amalgkit_pbarbatus_limited.yaml`**: Resource-limited configuration for Pbarbatus analysis on constrained systems
- **`amalgkit_pbarbatus_manual.yaml`**: Manual step-by-step configuration for Pbarbatus workflows
- **`amalgkit_pbarbatus_multi8.yaml`**: Multi-threaded (8 threads) configuration for high-performance Pbarbatus analysis
- **`amalgkit_pbarbatus_simple.yaml`**: Simplified configuration for basic Pbarbatus analysis workflows
- **`amalgkit_pbarbatus_test.yaml`**: Testing configuration for validation and development of Pbarbatus workflows

#### Template Configurations
- **`amalgkit_template.yaml`**: Template configuration file for creating new amalgkit workflows with customizable parameters

## Usage

These configuration files are used by the RNA analysis workflows and can be customized for different datasets, species, and computational requirements.

### Loading Configurations
```python
from metainformant.core import config
from metainformant.rna import AmalgkitWorkflowConfig

# Load configuration from file
cfg = config.load_config("config/amalgkit_pbarbatus.yaml")

# Create workflow configuration
workflow_cfg = AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    threads=cfg.compute.threads,
    species_list=cfg.species,
    **cfg.amalgkit_params
)
```

### Command Line Usage
```bash
# Run workflow with specific configuration
uv run python -m metainformant rna run \
  --config config/amalgkit_pbarbatus_multi8.yaml \
  --work-dir output/amalgkit/pbarbatus
```

## Configuration Parameters

### Common Parameters
- **`species_list`**: List of species for analysis (e.g., `["Apis_mellifera", "Apis_cerana"]`)
- **`work_dir`**: Working directory for analysis outputs
- **`threads`**: Number of parallel threads for processing
- **`memory_gb`**: Memory allocation in GB
- **`steps`**: List of workflow steps to execute

### Species-Specific Configurations
- **`Apis_barbatus` (Pbarbatus)**: Specialized configurations for the African honey bee species
- **`Apis_mellifera`**: Western honey bee configurations
- **`Apis_cerana`**: Eastern honey bee configurations

### Performance Profiles
- **Direct**: Streamlined execution with minimal overhead
- **Limited**: Resource-constrained environments
- **Multi-threaded**: High-performance parallel processing
- **Simple**: Basic functionality for quick analysis
- **Test**: Validation and development workflows

## Integration

- Used by `src/metainformant/rna/` module for workflow configuration and execution
- Referenced in documentation examples and tutorials
- Can be overridden by environment variables or custom configurations
- Supports integration with SLURM, PBS, and other job schedulers

## Modifying Configurations

When creating or modifying these files:
1. Test changes with small datasets first
2. Validate against existing test suites
3. Update corresponding documentation
4. Consider backward compatibility
5. Document any new parameters added
6. Add appropriate test configurations

## Best Practices

### Configuration Design
- Use descriptive parameter names
- Include comments explaining complex parameters
- Provide sensible defaults for optional parameters
- Group related parameters logically
- Support environment variable overrides

### Testing Configurations
- Create test-specific configurations for validation
- Include sample datasets for configuration testing
- Document expected outputs for each configuration
- Test both success and failure scenarios

### Version Control
- Track configuration changes in version control
- Document breaking changes in release notes
- Maintain backward compatibility when possible
- Use semantic versioning for configuration formats

## Related Documentation

- See `src/metainformant/rna/README.md` for RNA workflow configuration details
- See `docs/rna/configs.md` for detailed configuration parameter explanations
- See `docs/rna/workflow.md` for workflow execution patterns
- See `tests/test_rna_workflow_config.py` for configuration testing examples
