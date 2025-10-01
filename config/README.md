# Configuration Directory

This directory contains repository-level configuration files used by METAINFORMANT for various tools and workflows.

## Files

- **`amalgkit_pbarbatus_multi8.yaml`**: Multi-threaded configuration for amalgkit RNA-seq analysis with Pbarbatus dataset
- **`amalgkit_pbarbatus.yaml`**: Standard configuration for amalgkit RNA-seq analysis with Pbarbatus dataset
- **`amalgkit_template.yaml`**: Template configuration file for amalgkit workflows

## Usage

These configuration files are used by the RNA analysis workflows and can be customized for different datasets and computational requirements.

## Integration

- Used by `src/metainformant/rna/` module for workflow configuration
- Referenced in documentation and example workflows
- Can be overridden by environment variables or custom configurations

## Modifying Configurations

When modifying these files:
1. Test changes with small datasets first
2. Update corresponding documentation
3. Consider backward compatibility
4. Document any new parameters added

## Related Documentation

- See `src/metainformant/rna/README.md` for RNA workflow configuration details
- See `docs/rna/configs.md` for detailed configuration parameter explanations
