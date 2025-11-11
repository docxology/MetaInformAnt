# METAINFORMANT Documentation Guide

This guide provides comprehensive navigation and best practices for using METAINFORMANT documentation.

## Documentation Structure

METAINFORMANT documentation is organized hierarchically by domain and functionality:

```
docs/
├── README.md                    # Documentation overview and navigation
├── index.md                    # Hierarchical navigation index
├── architecture.md            # System design and architecture
├── cli.md                      # Command-line interface reference
├── setup.md                    # Installation and environment setup
├── testing.md                  # Testing documentation
├── core/                       # Core utilities documentation
│   ├── README.md              # Core module overview
│   ├── config.md              # Configuration management
│   ├── io.md                  # Input/output utilities
│   └── ...                    # Other core utilities
├── dna/                        # DNA analysis documentation
├── rna/                        # RNA analysis documentation
├── gwas/                       # GWAS documentation
└── ...                         # Other domain modules
```

## Navigation Strategies

### By Domain

If you know which biological domain you're working with:

1. **Start with domain index**: `docs/<domain>/index.md`
   - Provides overview of domain capabilities
   - Links to all domain-specific documentation
   - Examples: `docs/dna/index.md`, `docs/rna/index.md`

2. **Check domain README**: `docs/<domain>/README.md`
   - Comprehensive module documentation
   - API reference and usage examples
   - Integration patterns

3. **Explore specific topics**: `docs/<domain>/<topic>.md`
   - Focused documentation on specific functionality
   - Examples: `docs/dna/alignment.md`, `docs/rna/workflow.md`

### By Task

If you have a specific task or question:

1. **Quick Start**: Start with `QUICKSTART.md` in repository root
2. **Setup**: See `docs/setup.md` for environment configuration
3. **CLI Usage**: See `docs/cli.md` for command-line interface
4. **Architecture**: See `docs/architecture.md` for system design
5. **Testing**: See `docs/testing.md` for testing guidelines

### By Module

If you're looking for module-specific information:

1. **Source Documentation**: `src/metainformant/<module>/README.md`
   - Detailed module implementation
   - Function-level documentation
   - Internal architecture

2. **User Documentation**: `docs/<module>/`
   - User-facing documentation
   - Workflow guides
   - Integration examples

## Documentation Types

### Overview Documents

- **`README.md`**: Comprehensive module documentation with examples
- **`index.md`**: Navigation and quick reference
- **`AGENTS.md`**: AI contribution documentation (if applicable)

### Topic-Specific Documents

- **`<topic>.md`**: Focused documentation on specific functionality
- Examples: `alignment.md`, `workflow.md`, `visualization.md`

### Reference Documents

- **`cli.md`**: Command-line interface reference
- **`architecture.md`**: System design and architecture
- **`testing.md`**: Testing guidelines and practices

## Finding Information

### Quick Reference

| What You Need | Where to Look |
|--------------|---------------|
| Installation | `QUICKSTART.md`, `docs/setup.md` |
| Module overview | `docs/<module>/index.md` or `docs/<module>/README.md` |
| API reference | `src/metainformant/<module>/README.md` |
| Workflow guide | `docs/<module>/workflow.md` or `docs/<module>/README.md` |
| Configuration | `docs/<module>/config.md` or `docs/core/config.md` |
| Examples | `docs/<module>/README.md`, `docs/<module>/examples/` |
| CLI commands | `docs/cli.md` |
| Testing | `docs/testing.md` |

### Search Strategies

1. **Use the main index**: `docs/index.md` provides hierarchical navigation
2. **Check module README**: Each module has comprehensive documentation
3. **Look for examples**: Most modules include usage examples in README files
4. **Check integration docs**: Cross-module integration patterns are documented

## Reading Documentation

### Recommended Reading Order

1. **Start with overview**: Read `docs/<module>/index.md` or `docs/<module>/README.md`
2. **Review examples**: Look at code examples in the README
3. **Check specific topics**: Read topic-specific docs as needed
4. **Explore integration**: See how modules work together

### Understanding Code Examples

All code examples in documentation:
- Write outputs to `output/` by default
- Use real implementations (no mocks)
- Include error handling
- Follow project conventions

### Following Cross-References

Documentation uses consistent cross-reference patterns:
- **Relative paths**: `./<file>.md` for same directory
- **Domain paths**: `../<domain>/<file>.md` for other domains
- **Source links**: `../../src/metainformant/<module>/README.md` for source docs

## Module Documentation Patterns

### Standard Structure

Each module follows a consistent documentation structure:

1. **Overview**: What the module does
2. **Key Components**: Main functionality areas
3. **Usage Examples**: Practical code examples
4. **Integration**: How to use with other modules
5. **API Reference**: Function and class documentation
6. **Configuration**: Configuration options and environment variables
7. **Testing**: Testing approach and coverage

### Documentation Files

- **`index.md`**: Quick navigation and overview
- **`README.md`**: Comprehensive module documentation
- **`AGENTS.md`**: AI contribution documentation (if applicable)
- **`<topic>.md`**: Topic-specific documentation

## Best Practices

### For Readers

1. **Start broad, then narrow**: Begin with overview, then dive into specifics
2. **Follow examples**: Code examples show real usage patterns
3. **Check integration docs**: Understand how modules work together
4. **Read error messages**: Documentation often explains common errors

### For Contributors

1. **Update existing docs**: Don't create new root-level documentation
2. **Follow structure**: Use established documentation patterns
3. **Include examples**: All functionality should have usage examples
4. **Cross-reference**: Link to related documentation
5. **Keep current**: Update docs when code changes

## Common Documentation Locations

### Core Infrastructure

- **Configuration**: `docs/core/config.md`
- **I/O Operations**: `docs/core/io.md`
- **Logging**: `docs/core/logging.md`
- **Path Handling**: `docs/core/paths.md`

### Domain Modules

- **DNA**: `docs/dna/`
- **RNA**: `docs/rna/`
- **GWAS**: `docs/gwas/`
- **Protein**: `docs/protein/`
- **Networks**: `docs/networks/`
- **Multi-Omics**: `docs/multiomics/`
- **Single-Cell**: `docs/singlecell/`
- **Information Theory**: `docs/information/`
- **Life Events**: `docs/life_events/`

### Specialized Domains

- **Math**: `docs/math/`
- **ML**: `docs/ml/`
- **Quality**: `docs/quality/`
- **Visualization**: `docs/visualization/`
- **Simulation**: `docs/simulation/`
- **Ontology**: `docs/ontology/`
- **Phenotype**: `docs/phenotype/`
- **Epigenome**: `docs/epigenome/`
- **Ecology**: `docs/ecology/`

## Getting Help

If you can't find what you're looking for:

1. **Check the index**: `docs/index.md` for hierarchical navigation
2. **Search module README**: Comprehensive module documentation
3. **Review examples**: Look at `docs/<module>/examples/` or code examples
4. **Check source**: `src/metainformant/<module>/README.md` for detailed API docs
5. **Review tests**: `tests/test_<module>_*.py` for usage examples

## Documentation Maintenance

Documentation is maintained alongside code:
- Updated with each major release
- Reviewed for accuracy and completeness
- Cross-references validated
- Examples tested for correctness

## Summary

METAINFORMANT documentation is comprehensive and well-organized:
- **Hierarchical structure** by domain and functionality
- **Consistent patterns** across all modules
- **Practical examples** for all functionality
- **Clear navigation** through indexes and cross-references

Start with `docs/index.md` or `docs/README.md` for an overview, then navigate to specific modules and topics as needed.


