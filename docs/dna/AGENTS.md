# Agent Directives: docs/dna

## Role

Documentation for the DNA sequence analysis module covering sequence processing, alignment algorithms (Needleman-Wunsch, Smith-Waterman, ClustalW), phylogenetics (UPGMA, NJ, maximum likelihood), population genetics (Tajima's D, Fst, heterozygosity, PCA, LD), variant analysis (SNP calling, filtering, VEP/SnpEff annotation), and NCBI genome retrieval.

## Module Scope

End-to-end DNA workflows:
- **Sequence operations**: FASTA/FASTQ parsing, reverse complement, translation (6-frame), k-mer spectra, repeat detection
- **Pairwise alignment**: Global (Needleman-Wunsch), local (Smith-Waterman), affine gap penalties, banded alignment
- **Multiple sequence alignment**: ClustalW-style progressive alignment, MAFFT wrapper, sequence weighting
- **Phylogenetics**: Distance-based (UPGMA, NJ), maximum likelihood (IQ-TREE wrapper), Bayesian (MrBayes), bootstrap support
- **Population genetics**: Diversity metrics (π, θ), differentiation (Fst, Dxy), selection (Tajima's D, iHS), haplotype analysis
- **Variant calling**: BAM → VCF (FreeBayes, GATK HaplotypeCaller), filtering (QUAL, DP, MQ), annotation (VEP, SnpEff)
- **Genome retrieval**: ENA/SRA accession download, indexing (BWA, STAR, kallisto), reference caching

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/dna/sequence/` | Sequence processing (FASTA I/O, GC%, k-mers, reverse complement, translation, codon tables) |
| `src/metainformant/dna/alignment/` | Pairwise alignment (NW, SW), MSA (ClustalW), distance matrices, scoring matrices (BLOSUM62, PAM250) |
| `src/metainformant/dna/phylogeny/` | Tree construction (UPGMA, NJ, ML), bootstrap resampling, Newick I/O, tree visualization |
| `src/metainformant/dna/population/` | Population stats (Tajima's D, Fst, π), PCA, LD decay, relatedness, admixture |
| `src/metainformant/dna/variants/` | VCF parsing, variant calling (BAM → VCF), filtering, VEP/SnpEff annotation, effect prediction |
| `scripts/dna/` | CLI wrappers: `call_variants.py`, `build_tree.py`, `compute_fst.py` |
| `tests/dna/test_dna_*.py` | Test suites for each submodule |

## Cross-Module Dependencies

### Upstream (dependencies on core/math)
- **Core I/O**: `metainformant.core.io` — FASTA/VCF/BAM reading, serialization
- **Core config**: `metainformant.core.utils.config` — Pipeline parameters
- **Math population**: `metainformant.math.fst`, `metainformant.math.coalescent` — Theoretical formulas

### Downstream (modules that depend on DNA)
- **Visualization**: `metainformant.visualization.phylogeny` — Tree plotting, Manhattan plots, heatmaps
- **GWAS**: `metainformant.gwas.analysis` — Variant association (uses VCF output from variants module)
- **RNA**: `metainformant.rna.retrieval` — ENA genome download (shared with dna.ncbi module)
- **Quality**: `metainformant.quality.assembly` — Assembly QC uses DNA sequence metrics
- **Multiomics**: `metainformant.multiomics.integration` — Cross-domain variant harmonization

## Rule Constraints

1. **Type hints mandatory** — All function signatures must include full type annotations.
2. **No test mocks** — Use real genomic data; synthetic only for edge cases (e.g., invalid FASTA).
3. **Reproducible examples** — Every doc page must have runnable code with expected output.
4. **NCBI citation** — When using ENA/SRA retrieval, include dataset accession and citation template.
5. **Coordinate systems** — Clearly document 0-based vs 1-based; BED is 0-based, VCF is 1-based.

## Related Tasks & Guides
- **Task docs**: [../../tasks/analyze_dna.md](../tasks/analyze_dna.md) — Quick reference for common DNA operations
- **Workflow**: `docs/dna/example_workflows.md` — End-to-end population genomics pipeline
- **API spec**: `docs/dna/SPEC.md` — Detailed function signatures and module boundaries
- **Community**: See also `docs/gwas/` (association), `docs/rna/` (expression), `docs/visualization/` (plotting)

## Maintenance Notes
- **Version**: Align releases with core module (semver tracked in `pyproject.toml`).
- **Backwards compatibility**: Sequence API stable; alignment parameters may extend without breaking.
- **Test coverage**: Aim for 90%+ branch coverage; see `tests/coverage/`.
- **Benchmarks**: Run `pytest tests/benchmark/` before major performance changes.

## AI Assistant Guidance
When updating DNA module docs:
1. First check `src/metainformant/dna/` for implementation changes
2. Update corresponding `docs/dna/*.md` page (sequence, alignment, etc.)
3. Add runnable example + expected stdout to every code block
4. Cross-link to related modules (rna, gwas, visualization) at page bottom
5. If adding new submodule (e.g., `dna/epigenetics/`), create `AGENTS.md`, `SPEC.md`, update `index.md`
