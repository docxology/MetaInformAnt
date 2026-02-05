# Specification: cursorrules

## Scope

Module-specific Cursor AI rules for consistent development patterns across METAINFORMANT. Contains domain-specific coding guidelines that are automatically loaded by Cursor AI when working in corresponding module directories.

## Architecture

- **Dependency Level**: Development Tooling
- **Component Type**: AI Assistant Rules
- **Integration**: Cursor AI Editor

### File Structure
```
cursorrules/
├── README.md               # Cursorrules system documentation
├── core.cursorrules        # Core infrastructure patterns
├── dna.cursorrules         # DNA sequence analysis patterns
├── rna.cursorrules         # RNA-seq and amalgkit patterns
├── gwas.cursorrules        # GWAS pipeline patterns
└── {module}.cursorrules    # Domain-specific patterns
```

## Data Structures

### Cursorrules File Format
Each `.cursorrules` file is a plain text file containing:
- Module-specific coding patterns
- Import conventions for the domain
- Required utility usage patterns
- Domain-specific validation rules
- Common pitfalls to avoid

### Available Cursorrules
| File | Domain |
|------|--------|
| core.cursorrules | Core infrastructure (I/O, config, logging) |
| dna.cursorrules | DNA sequence analysis, alignment |
| rna.cursorrules | RNA-seq, amalgkit workflows |
| gwas.cursorrules | GWAS pipelines, association |
| protein.cursorrules | Protein analysis, structure |
| epigenome.cursorrules | Methylation, ChIP-seq |
| networks.cursorrules | Biological networks |
| multiomics.cursorrules | Multi-omic integration |
| singlecell.cursorrules | Single-cell RNA-seq |
| visualization.cursorrules | Plotting, figures |
| quality.cursorrules | Quality control |
| ml.cursorrules | Machine learning |
| math.cursorrules | Population genetics theory |
| information.cursorrules | Information theory |
| ontology.cursorrules | GO analysis |
| phenotype.cursorrules | Trait analysis |
| ecology.cursorrules | Community diversity |
| simulation.cursorrules | Synthetic data |
| life_events.cursorrules | Event sequences |
| longread.cursorrules | Long-read sequencing (PacBio/Nanopore) |
| metagenomics.cursorrules | Amplicon, shotgun metagenomics |
| structural_variants.cursorrules | CNV/SV detection and annotation |
| spatial.cursorrules | Spatial transcriptomics |
| pharmacogenomics.cursorrules | Clinical pharmacogenomics |
| menu.cursorrules | Interactive menu system |

## Interface

### Usage
Cursorrules are automatically loaded by Cursor AI based on the working directory. No manual invocation required.

### Key Patterns Enforced
All cursorrules enforce these project-wide standards:
```
- Use metainformant.core.io for file operations
- Use metainformant.core.logging for logging
- Write outputs to output/ directory
- Use absolute imports from metainformant
- Never use mocks in tests
- Use uv for package management
- Use real implementations only (NO MOCKING policy)
```

### Creating New Cursorrules
1. Create `{module}.cursorrules` in this directory
2. Document module-specific import patterns
3. Define required utility usage
4. Include common patterns and anti-patterns
5. Reference related cursorrules files
