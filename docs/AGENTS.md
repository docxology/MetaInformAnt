# Agent Directives: docs

## Role
Documentation agent context for METAINFORMANT's technical documentation.

## Directory Structure
Documentation is organized by module domain:
- `core/` - Core infrastructure documentation (I/O, config, paths, logging)
- `dna/` - DNA sequence analysis documentation
- `rna/` - RNA-seq and amalgkit workflow documentation
- `gwas/` - GWAS pipeline documentation
- `visualization/` - Plotting and visualization documentation
- `{module}/` - Domain-specific documentation for each module

## Key Root Files
- `architecture.md` - System architecture overview
- `cli.md` - Command-line interface documentation
- `testing.md` - Testing guidelines and patterns
- `FAQ.md` - Frequently asked questions
- `TUTORIALS.md` - Step-by-step tutorials
- `ERROR_HANDLING.md` - Error handling patterns
- `NO_MOCKING_POLICY.md` - Testing policy (no mocks allowed)
- `UV_SETUP.md` - UV package manager setup guide

## Rules and Constraints

### Documentation Standards
- Use Markdown format for all documentation
- Include code examples that are REAL and RUNNABLE
- Keep documentation synchronized with code
- Cross-reference related documentation

### File Organization
- `index.md` - Entry point for each domain subdirectory
- `workflow.md` - Workflow documentation when applicable
- Module-specific `.md` files for detailed API documentation

### Update Policy
- NEVER create new root-level documentation files
- Update existing docs or add to appropriate `docs/{domain}/` subdirectory
- Keep documentation DRY (Don't Repeat Yourself)

## Sphinx Integration
- `conf.py` - Sphinx configuration for documentation generation
- Run `bash scripts/package/uv_docs.sh` to build documentation
