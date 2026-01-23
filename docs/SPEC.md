# Specification: docs

## Scope

Technical documentation for METAINFORMANT, organized by domain module. Contains architecture guides, API documentation, tutorials, and reference materials for all 19 bioinformatics modules. Built with Sphinx for HTML generation.

## Architecture

- **Dependency Level**: Documentation
- **Component Type**: Technical Documentation
- **Build System**: Sphinx with markdown support

### Directory Structure
```
docs/
├── conf.py                  # Sphinx configuration
├── index.md                 # Documentation entry point
├── architecture.md          # System architecture overview
├── cli.md                   # CLI documentation
├── FAQ.md                   # Frequently asked questions
├── TUTORIALS.md             # Step-by-step tutorials
├── ERROR_HANDLING.md        # Error handling patterns
├── NO_MOCKING_POLICY.md     # Testing policy
├── UV_SETUP.md              # UV package manager guide
└── {module}/                # Domain-specific documentation
    ├── index.md             # Module entry point
    ├── API.md               # API reference
    ├── ARCHITECTURE.md      # Module architecture
    └── workflow.md          # Workflow documentation
```

## Data Structures

### Documentation Types
- **index.md**: Entry points for domains and subdirectories
- **API.md**: Function signatures, parameters, return types
- **ARCHITECTURE.md**: Module design and component relationships
- **workflow.md**: Step-by-step workflow documentation
- **TUTORIALS.md**: Hands-on tutorials with real examples
- **FAQ.md**: Common questions and troubleshooting

### Domain Subdirectories
Each module has a corresponding docs subdirectory:
- core/, dna/, rna/, gwas/, protein/, epigenome/
- networks/, multiomics/, singlecell/, visualization/
- quality/, ml/, math/, information/, ontology/
- phenotype/, ecology/, simulation/, life_events/

## Interface

### Building Documentation
```bash
# Build HTML documentation
bash scripts/package/uv_docs.sh

# Generated docs output to docs/_build/html/
```

### Adding Documentation
1. Create or update files in appropriate `docs/{module}/` subdirectory
2. Use Markdown format with code examples
3. Cross-reference related documentation
4. Include REAL, RUNNABLE code examples

### Documentation Standards
- Markdown format for all documentation
- Code examples must be real and executable
- Keep documentation synchronized with source code
- Never create new root-level files (update existing or add to subdirectories)
