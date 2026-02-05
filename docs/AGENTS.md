# Agent Directives: docs

## Role
Documentation agent context for METAINFORMANT's technical documentation.

## Directory Structure
Documentation is organized by module domain:
- `core/` - Core infrastructure (I/O, config, paths, logging, parallel, caching)
- `dna/` - DNA sequence analysis (alignment, phylogeny, population genetics)
- `rna/` - RNA-seq and amalgkit workflow documentation
- `gwas/` - GWAS pipeline (association, QC, visualization)
- `protein/` - Protein analysis (sequences, structures, AlphaFold, UniProt, InterPro)
- `epigenome/` - Epigenomics (methylation, ChIP-seq, ATAC-seq)
- `networks/` - Biological networks (PPI, regulatory, community detection, pathways)
- `multiomics/` - Multi-omic integration methods
- `singlecell/` - Single-cell RNA-seq (preprocessing, clustering, trajectory)
- `visualization/` - Plotting and visualization (80+ plot types)
- `quality/` - QC metrics and assessment
- `ml/` - Machine learning pipelines
- `math/` - Population genetics theory (coalescent, FST, LD, selection)
- `information/` - Information theory (entropy, mutual information)
- `ontology/` - Gene Ontology analysis and semantic similarity
- `phenotype/` - Trait analysis
- `ecology/` - Community diversity and ecology
- `simulation/` - Synthetic data generation
- `life_events/` - Event sequence analysis
- `metagenomics/` - Microbiome (amplicon, shotgun, functional annotation)
- `structural_variants/` - CNV/SV detection, annotation, visualization
- `longread/` - Long-read sequencing (PacBio/Nanopore, modified bases, assembly)
- `spatial/` - Spatial transcriptomics (Visium, MERFISH, Xenium)
- `pharmacogenomics/` - Clinical variants (CPIC, PharmGKB, ACMG, star alleles)
- `menu/` - Interactive menu and discovery system

## Key Root Files
- `architecture.md` - System architecture overview
- `cli.md` - Command-line interface documentation
- `testing.md` - Testing guidelines and patterns
- `FAQ.md` - Frequently asked questions
- `TUTORIALS.md` - Step-by-step tutorials
- `ERROR_HANDLING.md` - Error handling patterns
- `NO_MOCKING_POLICY.md` - Testing policy (no mocks allowed)
- `UV_SETUP.md` - UV package manager setup guide
- `setup.md` - Environment setup instructions
- `DISK_SPACE_MANAGEMENT.md` - Disk space management for large datasets

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
