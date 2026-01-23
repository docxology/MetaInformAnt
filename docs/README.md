# Documentation

METAINFORMANT documentation organized by biological domain.

## Navigation

| Entry Point | Description |
|-------------|-------------|
| [index.md](index.md) | Complete documentation index |
| [DOCUMENTATION_GUIDE.md](DOCUMENTATION_GUIDE.md) | How to navigate docs |
| [architecture.md](architecture.md) | System architecture |
| [cli.md](cli.md) | CLI reference |
| [setup.md](setup.md) | Installation guide |
| [testing.md](testing.md) | Testing documentation |

## Domain Documentation

| Domain | Directory | Description |
|--------|-----------|-------------|
| Core | [core/](core/) | Shared utilities (I/O, config, logging) |
| DNA | [dna/](dna/) | Sequence analysis, alignment, phylogeny |
| RNA | [rna/](rna/) | RNA-seq workflows, amalgkit |
| Protein | [protein/](protein/) | Protein sequences, structures |
| GWAS | [gwas/](gwas/) | Association studies |
| Visualization | [visualization/](visualization/) | Plotting, figures |
| Networks | [networks/](networks/) | Biological networks |
| Single-Cell | [singlecell/](singlecell/) | scRNA-seq analysis |
| Multi-Omics | [multiomics/](multiomics/) | Data integration |
| ML | [ml/](ml/) | Machine learning |
| Math | [math/](math/) | Population genetics theory |
| Information | [information/](information/) | Information theory |
| Life Events | [life_events/](life_events/) | Event sequence analysis |
| Ontology | [ontology/](ontology/) | GO analysis |
| Phenotype | [phenotype/](phenotype/) | Trait analysis |
| Ecology | [ecology/](ecology/) | Community ecology |
| Epigenome | [epigenome/](epigenome/) | Methylation, ChIP-seq |
| Simulation | [simulation/](simulation/) | Synthetic data |
| Quality | [quality/](quality/) | QC metrics |

## Key Guides

- [UV_SETUP.md](UV_SETUP.md) - Package management with UV
- [NO_MOCKING_POLICY.md](NO_MOCKING_POLICY.md) - Testing philosophy
- [ERROR_HANDLING.md](ERROR_HANDLING.md) - Error handling patterns
- [FAQ.md](FAQ.md) - Frequently asked questions

## Contributing

When adding documentation:

1. Place in appropriate `docs/<domain>/` directory
2. Follow existing formatting patterns
3. Include code examples
4. Update cross-references
5. Never create docs in repository root or `output/`

## Related

- [README.md](../README.md) - Project overview
- [QUICKSTART.md](../QUICKSTART.md) - Quick start guide
- [CLAUDE.md](../CLAUDE.md) - AI assistant guide
- [SPEC.md](../SPEC.md) - Technical specification
