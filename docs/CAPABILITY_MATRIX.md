# MetaInformAnt capability matrix

This is the source-linked capability contract for the 0.4.0 release. A
function being importable does not make it a completed scientific workflow.
The matrix describes implementation readiness; biological inference requires a
separate current manifest, provenance-qualified inputs, statistical evidence,
and an evidence bundle.

| Surface | Source of truth | Capability | Release boundary |
| --- | --- | --- | --- |
| Core utilities, I/O, logging, configuration | `src/metainformant/core/` | Stable | Public package utilities and canonical imports only |
| DNA sequence and genome helpers | `src/metainformant/dna/` | Stable | NCBI access requires explicit contact or anonymous opt-in |
| RNA acquisition and Amalgkit orchestration | `src/metainformant/rna/engine/`, `src/metainformant/rna/amalgkit/` | Experimental / operational | Requires pinned Amalgkit capability, campaign-local SRA environment, locks, receipts, and resumable state |
| RNA quantification provenance and readiness | `src/metainformant/rna/engine/provenance.py`, `scripts/rna/check_pipeline_status.py` | Stable contract | Readable files without current provenance are not promoted |
| RNA cross-species comparison | `src/metainformant/rna/analysis/cross_species.py` | Experimental | Explicit sample alignment, orthology validation, and descriptive statistics are required |
| GWAS and BeeWAS client/reporting surfaces | `src/metainformant/gwas/`, `projects/apis_gwas/` | Experimental / nested release | Nested repository PR and provenance validators govern publication |
| Population genetics statistics and neutrality tests | `src/metainformant/popgen/` | Implemented / tested | Summary statistics, neutrality tests, Fst, genotype structure, LD summaries, demographic comparisons; tests under `tests/popgen/` |
| Sequencing quality control | `src/metainformant/quality/` | Implemented / tested | FASTQ analysis, contamination detection, composite quality scoring, QC reports; tests under `tests/quality/` |
| Machine learning pipelines | `src/metainformant/ml/` | Implemented / tested | Classification, regression, AutoML, feature engineering, dimensionality reduction, interpretability, local LLM inference; tests under `tests/ml/` |
| Multi-omics integration | `src/metainformant/multiomics/` | Implemented / tested | `MultiOmicsData` layer alignment, `integrate_omics_data`, joint PCA/NMF, canonical correlation; tests under `tests/multiomics/` |
| Biological networks | `src/metainformant/networks/` | Implemented / tested | Network construction, graph algorithms, community detection; tests under `tests/networks/` |
| Ontology (Gene Ontology) | `src/metainformant/ontology/` | Implemented / tested | GO parsing, querying, enrichment analysis, visualization; tests under `tests/ontology/` |
| Phenotype analysis | `src/metainformant/phenotype/` | Implemented / tested | Morphological, behavioral, and chemical traits, life-course trajectories, AntWiki integration; tests under `tests/phenotype/` |
| Community ecology | `src/metainformant/ecology/` | Implemented / tested | Diversity indices, ordination, species abundance distributions, functional traits; tests under `tests/ecology/` |
| Life-course event sequences | `src/metainformant/life_events/` | Implemented / tested | Temporal event modeling, embedding learning, survival analysis, trajectory comparison; tests under `tests/life_events/` |
| Long-read sequencing | `src/metainformant/longread/` | Implemented / tested | PacBio/ONT signal I/O, quality assessment, assembly, methylation calling, haplotype phasing, SV detection; tests under `tests/longread/` |
| Metagenomics | `src/metainformant/metagenomics/` | Implemented / tested | 16S/ITS amplicon profiling, shotgun analysis, community diversity, functional annotation, differential abundance; tests under `tests/metagenomics/` |
| Pharmacogenomics | `src/metainformant/pharmacogenomics/` | Implemented / tested | Star allele calling, metabolizer phenotyping, CPIC guideline lookups, drug interaction prediction; tests under `tests/pharmacogenomics/` |
| Spatial transcriptomics | `src/metainformant/spatial/` | Implemented / tested | Visium/MERFISH/Xenium I/O, spatial statistics, cell-cell communication, deconvolution; tests under `tests/spatial/` |
| Structural variants | `src/metainformant/structural_variants/` | Implemented / tested | CBS-based CNV detection, split/discordant-read SV calling, annotation, population genotyping; tests under `tests/structural_variants/` |
| Metabolomics | `src/metainformant/metabolomics/` | Implemented / tested | MGF/CSV readers and writers, metabolite identification, metabolite set enrichment, plots; no mzML/mzXML readers; tests under `tests/metabolomics/` |
| Cloud deployment helpers | `src/metainformant/cloud/` | Implemented / tested | GCP VM lifecycle management for pipeline workloads; tests under `tests/cloud/`; external GCP access requires credentials |
| MCP server | `src/metainformant/mcp/` | Implemented / tested | stdio JSON-RPC 2.0 server (initialize, tools/list, tools/call, resources), schema-validated tool registry; tests under `tests/mcp/` |
| Menu system | `src/metainformant/menu/` | Implemented / tested | Interactive CLI menu, script discovery, navigation, execution; tests under `tests/menu/` |
| Protein analysis | `src/metainformant/protein/` | Implemented / tested | Sequences, 3D structure, database integration, functional annotation; tests under `tests/protein/` |
| Epigenomics | `src/metainformant/epigenome/` | Implemented / tested | Methylation, ChIP-seq, ATAC-seq, chromatin state discovery; tests under `tests/epigenome/` |
| Single-cell RNA-seq | `src/metainformant/singlecell/` | Implemented / tested | Preprocessing, dimensionality reduction, clustering, trajectory inference; tests under `tests/singlecell/` |
| Synthetic data generation | `src/metainformant/simulation/` | Implemented / tested | Sequence evolution, population genetics, RNA-seq counts, agent-based ecosystem models; tests under `tests/simulation/` |
| Mathematical biology | `src/metainformant/math/` | Implemented / tested | Population genetics, epidemiology, evolutionary dynamics, quantitative genetics; tests under `tests/math/` |
| Information theory | `src/metainformant/information/` | Implemented / tested | Entropy, mutual information, estimation, geometry; tests under `tests/information/` |
| Visualization | `src/metainformant/visualization/` | Implemented / tested | 70+ plot types across the visualization subpackages; tests under `tests/visualization/` |
| eQTL pipelines | `src/metainformant/eqtl/` | Implemented / tested | Input construction, HISAT2+bcftools wrappers, stats parsing, pipeline orchestration; tests under `tests/eqtl/`; external tools must be present |
| Optional integrations and external services | Module-specific `README.md`, `SPEC.md`, and dependency declarations | Experimental | Missing optional dependencies must follow the documented contract |
| Scaffolds and placeholders | Source TODOs, explicit `NotImplementedError`, and module validation reports | Scaffold | Must not be presented as production or inference-ready |
| Removed 0.4.0 surfaces | `docs/MIGRATION_0.4.md` | Removed | No runtime compatibility shims; migration documentation is the sole reference |
| Descriptive matrices and permutation scores | RNA analysis outputs and their provenance | Descriptive only | Must not be labeled conserved/divergent biological inference without the required evidence |
| Biological inference and release claims | Current manifests, receipts, finalized matrices, and evidence bundle | Gated | Withheld while the producer is active or the cohort is partial |

## Promotion rules

1. Stable means the source, tests, documentation, and import contract agree;
   it does not guarantee external service availability.
2. Experimental means the implementation is usable under the documented
   contract but remains subject to fixture, external-tool, or scientific
   validation limits.
3. Scaffold means the surface is intentionally incomplete and cannot be used
   as evidence of a completed analysis.
4. Implemented / tested means the surface ships with source and an in-repo
   test directory under `tests/<domain>/`; it is a presence statement, not a
   stability contract. A Stable claim additionally requires the source, tests,
   documentation, and import contract to agree (rule 1).
5. A release claim may advance from executable readiness to cohort readiness,
   descriptive analysis, and biological inference only when each preceding
   state has its own current receipt and evidence.
