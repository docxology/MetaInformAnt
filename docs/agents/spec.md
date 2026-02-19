# Specification: cursorrules

## Scope

Module-specific Cursor AI rules for consistent development patterns across METAINFORMANT. Contains domain-specific coding guidelines that are automatically loaded by Cursor AI when working in corresponding module directories.

## Architecture

- **Dependency Level**: Development Tooling
- **Component Type**: AI Assistant Rules
- **Integration**: Cursor AI Editor

### File Structure

```text
docs/agents/
├── README.md               # Agents system documentation
├── spec.md                 # Specification
└── rules/                  # Module-specific rules
    ├── core.md             # Core patterns
    ├── dna.md              # DNA patterns
    └── {module}.md         # Domain-specific patterns
```

## Data Structures

### Agent Rule File Format

Each `.md` file is a markdown file containing:

- Module-specific coding patterns
- Import conventions for the domain
- Required utility usage patterns
- Domain-specific validation rules
- Common pitfalls to avoid

### Available Agent Rules

| File | Domain |
|------|--------|
| core.md | Core infrastructure (I/O, config, logging) |
| dna.md | DNA sequence analysis, alignment |
| rna.md | RNA-seq, amalgkit workflows |
| gwas.md | GWAS pipelines, association |
| protein.md | Protein analysis, structure |
| epigenome.md | Methylation, ChIP-seq |
| networks.md | Biological networks |
| multiomics.md | Multi-omic integration |
| singlecell.md | Single-cell RNA-seq |
| visualization.md | Plotting, figures |
| quality.md | Quality control |
| ml.md | Machine learning |
| math.md | Population genetics theory |
| information.md | Information theory |
| ontology.md | GO analysis |
| phenotype.md | Trait analysis |
| ecology.md | Community diversity |
| simulation.md | Synthetic data |
| life_events.md | Event sequences |
| longread.md | Long-read sequencing (PacBio/Nanopore) |
| metagenomics.md | Amplicon, shotgun metagenomics |
| structural_variants.md | CNV/SV detection and annotation |
| spatial.md | Spatial transcriptomics |
| pharmacogenomics.md | Clinical pharmacogenomics |
| metabolomics.md | Metabolomics analysis |
| menu.md | Interactive menu system |

## Interface

### Usage

Cursorrules are automatically loaded by Cursor AI based on the working directory. No manual invocation required.

### Key Patterns Enforced

All cursorrules enforce these project-wide standards:

```text
- Use metainformant.core.io for file operations
- Use metainformant.core.utils.logging for logging
- Write outputs to output/ directory
- Use absolute imports from metainformant
- Never use mocks in tests
- Use uv for package management
- Use real implementations only (NO MOCKING policy)
```

### Creating New Agent Identify Rules

1. Create `{module}.md` in `docs/agents/rules/`
2. Document module-specific import patterns
3. Define required utility usage
4. Include common patterns and anti-patterns
5. Reference related rule files

## 🧪 Testing Policy

- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
