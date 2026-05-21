# Methods Matrix: Comprehensive Comparison of All 28 METAINFORMANT Modules

## Overview

This matrix provides a side-by-side comparison of all 28 core modules in METAINFORMANT across key dimensions: data types, scale, computational characteristics, outputs, and best-use scenarios. Use this table to quickly identify which module(s) match your analysis needs.

## Quick Navigation

- **[Summary Table](#summary-table)** — 28-module overview at a glance
- **[By Data Type](#by-data-type)** — Find modules by input data
- **[By Analysis Goal](#by-analysis-goal)** — Find modules by research question
- **[By Computational Scale](#by-computational-scale)** — Find modules by resource needs
- **[Detailed Module Cards](#detailed-module-cards)** — Per-module deep dive
- **[Cross-Module Dependencies](#cross-module-dependencies)** — How modules interact

---

## Summary Table

| # | Module | Data Type | Scale | Primary Methods | Output Format | Computational Intensity | Best-Use Scenario |
|---|--------|-----------|-------|-----------------|---------------|------------------------|-------------------|
| 1 | **core** | Any | Any | I/O, config, logging, parallel, caching | Config files, logs | Minimal | Shared infrastructure — used by ALL modules |
| 2 | **dna** | FASTA, VCF, FASTQ | Single → population | Alignment (NW/SW), phylogeny (NJ/UPGMA), population genetics (π, Fst, Tajima's D) | Newick, distance matrices, annotated VCF | Low–Medium | Sequence analysis, evolutionary biology, variant annotation |
| 3 | **rna** | FASTQ, BAM | Bulk samples → populations | Amalgkit streaming workflow, quantification (Kallisto), TMM normalization, cross-species ortholog mapping | TPM/counts matrices, manifests, logs | Medium–High | Bulk RNA-seq meta-analysis, ENA/SRA large-scale integration |
| 4 | **gwas** | VCF, phenotypes | Population (10³–10⁶) | Association (linear/logistic/Mixed LM), QC (MAF, HWE), population structure (PCA, kinship), fine-mapping (SuSiE), colocalization | Association statistics (TSV), Manhattan/QQ plots | Medium–High | Genome-wide association, variant-trait mapping, heritability |
| 5 | **protein** | FASTA, PDB, UniProt IDs | Single → population | Sequence analysis, AlphaFold structure fetching, InterPro domains, PDB download, RMSD alignment | PDB/CIF, annotated FASTA, JSON | Low–Medium | Protein structures, domain annotation, structural comparison |
| 6 | **epigenome** | BAM, BedGraph, narrowPeak | Samples → cohorts | Methylation (bismark), ChIP-seq peak calling, ATAC-seq accessibility, chromatin states | BigWig, peak BED, methylation calls | Medium–High | Epigenetic modifications, chromatin accessibility |
| 7 | **singlecell** | h5ad, count matrices | Cells (10³–10⁶) | Preprocessing (filtering, normalization), clustering (Leiden), DE, trajectory (PAGA), RNA velocity, cell typing | AnnData, clusters, embeddings, marker lists | High | Single-cell transcriptomics, cell-type discovery, developmental trajectories |
| 8 | **spatial** | Visium/Xenium H5, images | Tissues (10–1000 spots) | Spatial clustering, domain detection, cell deconvolution (NNLS/NMF), spatial autocorrelation (Moran's I), ligand-receptor signaling | Spatial coordinates, tissue maps, DE genes per domain | Very High | Spatial transcriptomics, tissue architecture, niche analysis |
| 9 | **multiomics** | Mixed (DNA+RNA+Protein) | Coordinated samples (10–1000) | Multi-omics container (MultiOmicsData), integration (concatenation, joint PCA/CCA/NMF, autoencoders), cross-omic correlation | Joint components, integrated loadings, biomarkers | High | Multi-omic integration, cross-modal biomarkers, systems biology |
| 10 | **networks** | Edgelist, adjacency | Networks (10²–10⁵ nodes) | Community detection (Louvain, Leiden), centrality (degree, betweenness, eigenvector), PPI/regulatory network construction, pathway enrichment | Graph data, community assignments, centrality scores | Medium–High | Biological networks, pathway analysis, network-based prioritization |
| 11 | **ml** | Feature matrices | Samples (10²–10⁶) | Classification (RF, SVM, LR, NB), regression, feature selection (MI, RFE, LASSO, variance), AutoML, LLM integration (Ollama) | Trained models, predictions, feature rankings | Medium–High | Machine learning classification/regression, feature discovery |
| 12 | **math** | Equations, parameters | Theoretical | Population genetics theory (HWE, drift, Fst, coalescent), epidemiology (SIR/SEIR), evolutionary dynamics (game theory), Bayesian inference | Analytical solutions, simulated datasets | Low–Medium | Mathematical modeling, theory development, simulation parameters |
| 13 | **information** | Discrete/continuous variables | Features (10²–10⁴) | Shannon/Rényi/Tsallis entropy, mutual information, KL divergence, Jensen-Shannon, transfer entropy, semantic similarity | Entropy values, MI matrices, information profiles | Low–Medium | Feature dependence, redundancy, semantic similarity, info-theoretic selection |
| 14 | **ontology** | Gene lists, OBO files | Genes/pathways (10–10⁴) | GO graph parsing, semantic similarity (Resnik, Lin, Jiang-Conrath), functional enrichment, pathway annotation | Enrichment results, similarity scores, annotated gene lists | Low–Medium | Functional enrichment, gene set analysis, semantic similarity |
| 15 | **phenotype** | Measurements, images | Specimens (10–10⁴) | Morphometrics (allometry, shape indices), behavioral sequences (ethograms), chemical profiles (GC-MS), electronic tracking (GPS/RFID), acoustic analysis (FFT) | Phenotypic profiles, indices, trajectories | Medium | Phenotypic characterization, trait analysis, life-history strategies |
| 16 | **ecology** | Species tables, environmental | Communities | Diversity indices (Shannon, Simpson), species richness, phylogenetic diversity (PD), environmental correlation (CCA/RDA) | Diversity tables, ordination plots | Medium | Community ecology, biodiversity assessment, species-environment relationships |
| 17 | **simulation** | Parameters, seeds | Synthetic datasets | Sequence evolution (substitution models), population genetics sims (coalescent), agent-based models, benchmark generation | Simulated FASTQ/VCF/Counts, synthetic networks | Medium–High | Benchmarking, method validation, null models, synthetic data |
| 18 | **quality** | FASTQ, BAM, assemblies | Reads/contigs (10⁴–10⁸) | FASTQ quality metrics (Phred, GC, adapters), contamination detection, assembly validation (BUSCO, N50), batch effect diagnostics | QC reports (JSON/HTML), quality plots | Low–Medium | Data quality assessment, preprocessing validation, contamination screening |
| 19 | **visualization** | Any (plot-ready) | Any | 70+ plot types: basic (line, scatter, bar, pie), statistical (histogram, box, violin, QQ), genomics (Manhattan, volcano), trees, networks, dimred (PCA/UMAP), timeseries, animations | PNG/PDF/SVG/HTML, interactive Plotly | Low–High (depends on plot) | Publication-quality figures, interactive dashboards, animated sequences |
| 20 | **longread** | FAST5, POD5, BAM | Reads (10³–10⁶) | Basecalling (ont-guppy/dorado integration), quality metrics (N50, accuracy), modified base detection (5mC, 6mA), SV calling, haplotype phasing, assembly (overlap-layout-consensus) | Basecall FASTQ, methylation calls, SV VCF, assemblies | High | Long-read sequencing analysis (PacBio/ONT), structural variants, methylation |
| 21 | **metagenomics** | FASTQ (amplicon/shotgun) | Samples (10–10⁴) | 16S/ITS amplicon (OTU/ASV), taxonomic classification, assembly (metaSPAdes), binning (MetaBAT), functional annotation (ORF, KEGG), pathway reconstruction, community profiling | OTU tables, MAGs, annotated genes, KEGG modules | High | Microbiome analysis, metagenome assembly, functional profiling |
| 22 | **structural_variants** | BAM, VCF | Samples (10–10⁴) | CNV detection (read depth, segmentation), SV calling (split-read, discordant pair), breakpoint refinement, gene/regulatory overlap, functional impact, multi-caller merging (SURVIVOR) | SV VCF, CNV BED, annotated regions, Circos plots | High | Structural variant detection, fusion genes, TAD disruption, clinical SV |
| 23 | **pharmacogenomics** | VCF, diplotype data | Patients (10–10⁴) | Star allele calling (CYP enzymes), metabolizer phenotype prediction, CPIC guideline integration, PharmGKB annotations, ACMG pathogenicity scoring, drug-gene interaction checking | Clinical reports, dosing recommendations, evidence tables | Medium | Clinical pharmacogenomics, drug response prediction, variant interpretation |
| 24 | **metabolomics** | mzML, mzXML, CSV | Samples (10–10³) | Peak detection, metabolite identification (mz matching), quantification, normalization (PQN, log), KEGG/Reactome pathway mapping, metabolite set enrichment (MSEA), metabolite-gene integration | Feature tables, identified metabolites, enriched pathways | Medium–High | Mass spec metabolomics, pathway analysis, metabolite-gene integration |
| 25 | **menu** | Scripts (any) | Any | Script discovery, interactive CLI menus, breadcrumb navigation, argument prompting | Interactive terminal UI | Minimal | Interactive workflow discovery, script execution helper |
| 26 | **cloud** | GCP resources | Any | GCP VM lifecycle (create/delete), Docker container build/run, file transfer (gsutil/SCP), cost estimation, preemptible VM management | Cloud resources, logs | Very High | Cloud deployment, large-scale pipeline orchestration, cost-optimized compute |
| 27 | **eqtl** *(cross-cutting)* | VCF + expression | Samples (100–10⁴) | cis/trans eQTL scanning, colocalization ( coloc ), mediation analysis, transcriptome SNP calling (from RNA BAM) | eQTL summary stats, colocalization PP4, mediation results | High | Expression quantitative trait loci, GWAS→RNA mechanistic link |
| 28 | **life_events** | Event sequences, time series | Individuals (10–10⁴) | Event sequence embedding, temporal pattern discovery, survival analysis, outcome prediction, life course modeling | Event embeddings, survival curves, predicted outcomes | Medium | Longitudinal life event analysis, survival modeling, temporal pattern mining |

---

## By Data Type: Quick Reference

### FASTA / Sequence Data
**Primary module**: `dna` (direct)
- Use `dna.sequences` for parsing, translation, ORF finding
- Phylogeny: `dna.phylogeny` for trees
- Distance: `dna.distances` for evolutionary distances

### VCF / Variant Data
**Primary modules**: `dna`, `gwas`, `structural_variants`, `pharmacogenomics`
- `dna`: parse, annotate, effect prediction
- `gwas`: association testing (requires phenotype)
- `structural_variants`: CNV/SV detection & annotation
- `pharmacogenomics`: clinical interpretation

### FASTQ / Sequencing Reads
**Primary modules**: `dna` (basic QC), `rna` (bulk), `metagenomics` (community), `longread` (ONT/PacBio), `quality` (QC)
- `rna`: bulk RNA-seq quantification (amalgkit)
- `metagenomics`: 16S/shotgun metagenomics
- `longread`: long-read basecalling, SV, methylation
- `quality`: FastQC-like metrics

### Count Matrices / Expression
**Primary modules**: `rna` (bulk), `singlecell` (scRNA), `spatial` (spatial), `protein` (proteomics), `metabolomics` (metabolites)
- `rna`: TPM/counts from bulk RNA
- `singlecell`: cell-level expression + clustering
- `spatial`: spatially-resolved expression
- `protein`: protein abundance from MS
- `metabolomics`: metabolite intensities

### Phenotype / Trait Measurements
**Primary module**: `phenotype`
- Direct analysis: morphological, behavioral, chemical, acoustic phenotypes
- Integration: feed into `gwas` for association

### Multi-Omic Data (DNA+RNA+Protein+Metabolites)
**Primary module**: `multiomics`
- Container: `MultiOmicsData` for harmonization
- Integration: joint PCA, CCA, NMF, autoencoders

---

## By Analysis Goal: Decision Map

| Research Question | Recommended Module(s) | Workflow |
|-------------------|----------------------|----------|
| **"Which genetic variants affect my trait?"** | `gwas` | QC → PCA → Association → Correction → Viz |
| **"What is the population structure?"** | `dna` or `gwas` (PCA) | `dna.population` or `gwas.structure.compute_pca` |
| **"How diverse is my sample?"** | `dna` (nucleotide diversity), `ecology` (species diversity) | `dna.population.nucleotide_diversity()` or `ecology.analysis.shannon_index()` |
| **"How are these sequences related?"** | `dna` (phylogeny) | MSA → distance → tree (NJ/UPGMA) |
| **"What genes are differentially expressed?"** | `rna`, `singlecell`, or direct DE tools | Use `rna` outputs → DESeq2/edgeR |
| **"Which cell types are in my tissue?"** | `singlecell` (with optional `spatial` mapping) | QC → Normalize → Cluster → Annotate |
| **"Where are genes expressed spatially?"** | `spatial` | Preprocess → Cluster domains → Map genes |
| **"How do multiple omics layers interact?"** | `multiomics` + `eqtl` | Align samples → Joint factorization → Cross-omic associations |
| **"What biological pathways are enriched?"** | `ontology` | Gene list → GO enrichment → Semantic similarity |
| **"Which genes are in this network module?"** | `networks` | Build network → Community detection → Module genes |
| **"Can I predict phenotype from genotype?"** | `ml` (or `gwas` PRS) | Feature selection → Classifier/regressor → CV |
| **"Is my sample contaminated?"** | `quality` | FastQC metrics → adapter/contamination detection |
| **"How does this protein fold?"** | `protein` (AlphaFold) | `fetch_alphafold_model(uniprot_id)` |
| **"Which microbes are present?"** | `metagenomics` | QC → Taxonomic classification → Phylogeny |
| **"What structural variants exist?"** | `structural_variants` (from long-read/BAM) | SV calling → Filter → Annotate |
| **"What drug interacts with this variant?"** | `pharmacogenomics` | Diplotype → Metabolizer phenotype → CPIC guideline |
| **"What is the entropy of this DNA segment?"** | `information` | Shannon entropy, MI between positions |
| **"How to visualize my GWAS results?"** | `visualization.genomics` | `manhattan_plot()`, `qq_plot()` |
| **"How to create animated evolution demo?"** | `visualization.animations` | `animate_evolution()` |

---

## By Computational Scale: Resource Planning

### Low Resource (≤4 GB RAM, <1 hour)
| Module | Typical Task | Time | Memory |
|--------|--------------|------|--------|
| `core` | I/O, config | seconds | MB |
| `dna` (small FASTA) | Alignment of short sequences | <1 min | 1 GB |
| `math` | Theoretical calculations | <1 min | negligible |
| `information` | Entropy on 1000s features | <5 min | 2 GB |
| `ontology` | GO enrichment (100 genes) | <5 min | 1 GB |
| `quality` (single FASTQ) | QC metrics | <5 min | 4 GB |

### Medium Resource (4–32 GB RAM, 1–12 hours)
| Module | Typical Task | Time | Memory |
|--------|--------------|------|--------|
| `dna` (genome-scale) | Phylogeny of 1000s sequences | 1–4 h | 8 GB |
| `protein` | Structure download + alignment | 1–4 h | 8 GB |
| `ecology` | Diversity indices (100 samples) | 30 min | 4 GB |
| `networks` (medium) | Community detection (1000 nodes) | 2 h | 16 GB |
| `ml` (RF on 10K samples) | Classification with CV | 2–6 h | 16 GB |
| `gwas` (10K samples) | Association with BOLT-LMM | 4–12 h | 32 GB |

### High Resource (32–128 GB RAM, 12–72 hours)
| Module | Typical Task | Time | Memory |
|--------|--------------|------|--------|
| `rna` (bulk, 100 samples) | Amalgkit full workflow | 12–48 h | 64 GB |
| `singlecell` (50K cells) | Clustering + trajectory | 24–48 h | 64 GB |
| `metagenomics` (metagenome assembly) | MetaSPAdes assembly | 24–72 h | 128 GB |
| `multiomics` (100 samples × 3 omics) | Joint factorization | 12–24 h | 64 GB |
| `eqtl` (1000 samples) | cis/trans scanning | 24–48 h | 64 GB |

### Very High Resource (≥128 GB RAM, days–weeks)
| Module | Typical Task | Time | Memory |
|--------|--------------|------|--------|
| `spatial` (whole slide) | Domain detection + deconvolution | 3–7 days | 256 GB |
| `longread` (nanopore run) | Basecalling + assembly of 1M reads | 3–10 days | 256 GB |
| `gwas` (500K samples) | Large biobank GWAS | 1–2 weeks | 512 GB+ (cluster) |
| `cloud` orchestration | 10,000 genome processing | 1–4 weeks | N/A (distributed) |
| `metagenomics` (1000 metagenomes) | Aggregate analysis | 2–3 weeks | 512 GB+ |

---

## Detailed Module Cards

Below are per-module summaries with capabilities, inputs/outputs, and example use cases.

---

### Module 1: core

| Attribute | Details |
|-----------|---------|
| **Category** | Infrastructure |
| **Data Types** | Any (I/O, config, logging) |
| **Scale** | Universal (used by all other modules) |
| **Key Functions** | `load_json()`, `dump_json()`, `get_logger()`, `parallel.map()`, `cache()` |
| **Input Files** | Config YAML/JSON, any data files |
| **Output Files** | Logs, manifests, cached artifacts |
| **Submodules** | `io`, `config`, `parallel`, `validation`, `paths`, `utils` |
| **Dependencies** | None (base module) |
| **Used By** | ALL 27 domain modules |

**When to use**: Never called directly by end-users except as underlying infrastructure. All domain modules use `core.io` for file operations and `core.parallel` for multithreading.

---

### Module 2: dna

| Attribute | Details |
|-----------|---------|
| **Category** | Sequences & Population Genetics |
| **Data Types** | FASTA, FASTQ, VCF |
| **Scale** | 1 → 10,000+ sequences |
| **Key Methods** | Pairwise alignment (global/local), MSA, NJ/UPGMA trees, Hamming distance, nucleotide diversity (π), Tajima's D, Fst, mutation analysis, ORF finding, translation, restriction enzymes |
| **Input** | Sequences (`seqs.fasta`), VCF (`variants.vcf.gz`), FASTQ (quality) |
| **Output** | Newick trees (`.nwk`), distance matrices (`.npy`), annotated VCF, translation products, consensus |
| **Submodules** | `sequences`, `alignment`, `msa`, `phylogeny`, `population`, `distances`, `translation`, `transcription`, `mutations`, `variants`, `calling`, `annotation`, `ncbi`, `fastq` |
| **External Deps** | Optional: MAFFT, Clustal, PHYLIP |
| **Typical Runtime** | Seconds (small) → hours (10K seqs) |
| **Primary Use Cases** | Sequence alignment, evolutionary analysis, population genetics, variant annotation, cross-species translation |

**Key workflows**:
- Phylogeny: `FASTA → MSA → distance matrix → tree → Newick`
- Population: `VCF → allele frequencies → diversity statistics`
- Variants: `VCF → effect prediction → functional annotation`

---

### Module 3: rna

| Attribute | Details |
|-----------|---------|
| **Category** | Bulk Transcriptomics |
| **Data Types** | FASTQ, BAM |
| **Scale** | 1 → 10,000+ samples |
| **Key Methods** | Amalgkit workflow (plan/execute), Kallisto quantification, TMM normalization, CSTMM (cross-species TMM), ortholog mapping, metadata curation, sanity checks |
| **Input** | FASTQ files, ENA/SRA accessions, species profiles |
| **Output** | `merged_abundance.tsv` (TPM/counts), `workflow/` directory with per-step logs, `amalgkit.manifest.jsonl` |
| **Submodules** | `amalgkit` (CLI wrapper), `engine` (workflow steps), `retrieval` (ENA/SRA), `analysis` (post-quant) |
| **External Deps** | `amalgkit` CLI, `kallisto`, `fastp`, `samtools` |
| **Typical Runtime** | Hours (1 sample) → weeks (8,300+ samples, streaming) |
| **Primary Use Cases** | Bulk RNA-seq quantification, meta-analysis across studies/species, ENA-first large-scale processing |

**Key workflows**:
- Standard: `FASTQ → quant (kallisto) → merge (cstmm) → curate`
- ENA streaming: `metadata → integrate → getfastq → quant → merge` (orchestrated)
- Cross-species: `mouse genes → orthologs → human → harmonized expression`

**Notable**: Uses `workflow.plan_workflow()` and `workflow.execute_workflow()`; intermediate manifests track all artifacts.

---

### Module 4: gwas

| Attribute | Details |
|-----------|---------|
| **Category** | Association Genetics |
| **Data Types** | VCF (genotypes), phenotype TSV, covariates TSV |
| **Scale** | 1,000 → 1,000,000+ samples |
| **Key Methods** | Quality control (MAF, missingness, HWE), PCA/kinship for population structure, association testing (linear/logistic regression, BOLT-LMM, SAIGE), multiple testing correction (Bonferroni, FDR, genomic control), fine-mapping (SuSiE), colocalization, LD score regression, heritability estimation |
| **Input** | Genotype VCF/BGEN, phenotype file (`sample_id, trait`), covariates (`age,sex,PC1..`), genome reference |
| **Output** | `association.tsv` (chrom, pos, beta, se, p), Manhattan plot, QQ plot, regional plots, fine-mapped credible sets |
| **Submodules** | `analysis` (association), `structure` (PCA, kinship), `finemapping` (SuSiE, coloc), `visualization` (manhattan, qq, regional), `workflow` (pipeline orchestration) |
| **External Deps** | PLINK2, BOLT-LMM, SAIGE, flashpca, SUSiE |
| **Typical Runtime** | Linear model: 10K samples × 1M SNPs = 2–5 min (8 cores); BOLT-LMM: 30 min; Fine-mapping: 1–2 h |
| **Primary Use Cases** | GWAS of complex traits, fine-mapping causal variants, colocalization with eQTL, heritability estimation |

**Workflow**:
1. **Genome prep**: Download/build reference (if needed)
2. **Variant acquisition**: Use existing VCF or call variants
3. **QC**: Filter by MAF, missingness, HWE
4. **Population structure**: PCA + kinship matrix
5. **Association**: Linear (quantitative) or logistic (binary) regression; mixed model if related
6. **Correction**: Bonferroni (stringent) or FDR (balanced)
7. **Viz**: Manhattan, QQ, regional locus plots
8. **Export**: Results TSV, summary JSON

**Integration**: Output feeds `eqtl` for colocalization, `multiomics` for joint analysis, `ontology` for enrichment.

---

### Module 5: protein

| Attribute | Details |
|-----------|---------|
| **Category** | Protein Sequences & Structures |
| **Data Types** | FASTA (sequences), PDB/CIF (structures), UniProt IDs |
| **Scale** | 1 → 10,000+ proteins |
| **Key Methods** | FASTA parsing/validation, amino acid composition, k-mer frequencies, UniProt ID mapping & FASTA fetch, PDB structure retrieval, AlphaFold model fetching, InterPro domain annotation, RMSD (Kabsch alignment), structure analysis (CA coordinates) |
| **Input** | Protein FASTA, UniProt accessions, PDB IDs, InterPro queries |
| **Output** | Composition JSON, PDB/CIF files, RMSD values, domain annotations (InterPro entries) |
| **Submodules** | `sequences`, `uniprot`, `pdb`, `alphafold`, `interpro`, `structure`, `structure_io` |
| **External Deps** | Networked APIs: UniProt REST, PDB, AlphaFold DB (no local install needed) |
| **Typical Runtime** | Seconds (fetch) to minutes (alignment) |
| **Primary Use Cases** | Protein sequence analysis, structure retrieval (experimental or AlphaFold), domain annotation, structural alignment (RMSD) |

**Examples**:
- `fetch_alphafold_model("P69905", output_dir)` → download AlphaFold structure
- `compute_rmsd_kabsch(coords_A, coords_B)` → RMSD in Ångströms
- `map_ids_uniprot(["P12345"])` → cross-database ID mapping

---

### Module 6: epigenome

| Attribute | Details |
|-----------|---------|
| **Category** | Epigenomics |
| **Data Types** | BAM (aligned), BedGraph (methylation), narrowPeak (ChIP), ATAC BAM |
| **Scale** | Samples: 10 → 10,000+ |
| **Key Methods** | Bisulfite methylation calling (bismark), CpG/CHG/CHH analysis, ChIP-seq peak calling (MACS2), ATAC-seq peak calling, chromatin state segmentation, differential methylation analysis, accessibility quantification |
| **Input** | Aligned BAM (bisulfite or ChIP/ATAC), reference genome, sample metadata |
| **Output** | Methylation calls (CpG + %), peak BED files, BigWig tracks, differential peaks, chromatin state segmentations |
| **Submodules** | `assays` (assay-specific: methylation, chip, atac), `chromatin_state`, `peak_calling` |
| **External Deps** | `bismark`, `bwa` (bisulfite), `macs2` (ChIP), `bedtools` |
| **Typical Runtime** | Hours per sample; 100 samples = days |
| **Primary Use Cases** | DNA methylation profiling, histone modification mapping, chromatin accessibility, epigenome-wide association studies (EWAS) |

**Assays supported**:
- **WGBS/RRBS**: methylation % per CpG
- **ChIP-seq**: peak calling + signal tracks
- **ATAC-seq**: open chromatin peaks
- **ChIPmentation**: low-input ChIP

---

### Module 7: singlecell

| Attribute | Details |
|-----------|---------|
| **Category** | Single-Cell Transcriptomics |
| **Data Types** | h5ad (AnnData), count matrices (CSV/TSV), CellRanger outputs |
| **Scale** | Cells: 1,000 → 1,000,000+; Genes: 20K |
| **Key Methods** | Filtering (cells/genes by counts), normalization (log1p, SCTransform), PCA, neighborhood graph, Leiden clustering, differential expression (Wilcoxon, MAST), UMAP/t-SNE, PAGA trajectory, RNA velocity scVelo, cell-type annotation (marker-based or reference mapping) |
| **Input** | Raw counts matrix, h5ad file, CellRanger `filtered_feature_bc_matrix` |
| **Output** | Processed AnnData (`.h5ad`), cluster assignments (`.tsv`), marker genes, embeddings (UMAP coordinates), trajectory graphs |
| **Submodules** | `analysis` (clustering, DE, trajectory), `data` (I/O), `velocity` (RNA velocity), `celltyping` (annotation), `workflow` (pipeline orchestration) |
| **External Deps** | scanpy, scVI-tools, scvelo (optional) |
| **Typical Runtime** | 10K cells: 1–4 h; 100K cells: 4–12 h; 1M cells: 1–2 days |
| **Primary Use Cases** | Cell-type discovery, developmental trajectory inference, cell-state dynamics, spatial mapping (via integration with `spatial`) |

**Workflow**:
1. **QC & filtering**: Remove low-quality cells/genes
2. **Normalization**: Library size + log-transform or SCTransform
3. **Feature selection**: HVG selection
4. **PCA** → **Neighborhood graph** → **Leiden clustering**
5. **UMAP/t-SNE** visualization
6. **DE** to find cluster markers
7. **Trajectory** (optional): PAGA/velocity
8. **Annotation**: marker-based or reference mapping

---

### Module 8: spatial

| Attribute | Details |
|-----------|---------|
| **Category** | Spatial Transcriptomics |
| **Data Types** | 10x Visium (h5), MERFISH/Xenium (HDF5), TIFF images, AnnData with spatial coords |
| **Scale** | Spots: 1K → 100K; Tissues: 1 → 100 |
| **Key Methods** | Multi-platform I/O (Visium/Xenium/MERFISH), spatial clustering (SPARK, BayesSpace), spatial domain detection, cell-type deconvolution (NNLS, NMF from reference scRNA), spatial autocorrelation (Moran's I, Geary's C, Getis-Ord G), neighborhood enrichment (CellChat, Giotto), ligand-receptor signaling analysis, spatial visualization (tissue overlays, gene maps) |
| **Input** | Counts matrix + spatial coordinates (tissue positions), raw image (optional) |
| **Output** | Spatial coordinates + cluster labels, deconvolved cell-type proportions per spot, spatial autocorrelation statistics, ligand-receptor interaction scores, tissue overlay plots |
| **Submodules** | `io` (platform readers), `analysis` (clustering, autocorrelation, deconvolution), `integration` (scRNA mapping), `visualization` (spatial plots, overlays) |
| **External Deps** | scanpy (preprocessing), squidpy (spatial stats), stlearn, giotto (optional) |
| **Typical Runtime** | Visium (5K spots): 1–4 h; Xenium (100K spots, 500 genes): 8–24 h |
| **Primary Use Cases** | Spatial transcriptomics analysis, tissue domain detection, cell-type deconvolution, spatially-resolved ligand-receptor interactions |

**Platforms supported**:
- **10x Visium**: Whole transcriptome, spot-based (55 µm)
- **10x Xenium**: Subcellular resolution, targeted panels
- **MERFISH**: 100s–1000s genes, high-plex

**Key outputs**: tissue image + spot coordinates + expression + cluster labels

---

### Module 9: multiomics

| Attribute | Details |
|-----------|---------|
| **Category** | Multi-Omic Integration |
| **Data Types** | Any combination: VCF (DNA), counts (RNA), protein abundances, metabolite intensities |
| **Scale** | Samples: 10 → 1,000+ (coordinated); Features: 10⁴–10⁶ per omic |
| **Key Methods** | Multi-omics data container (`MultiOmicsData`) with automatic sample matching, integration strategies: early (concatenation), intermediate (joint dimensionality reduction), late (ensemble), joint analysis: PCA, CCA, NMF, autoencoders (AE/VAE), multi-view clustering, cross-omics biomarker discovery, pathway-level integration |
| **Input** | Multiple omic tables (samples × features) with sample IDs; optional metadata |
| **Output** | Joint embeddings (samples × components), cross-omic loadings (features × components), ranked biomarker pairs, cluster assignments per omic |
| **Submodules** | `integration` (core algorithms), `analysis` (statistical tests), `pathways` (pathway integration), `visualization` (multi-omic plots) |
| **External Deps** | None (pure Python); uses numpy, scipy, sklearn, torch (optional for AEs) |
| **Typical Runtime** | 50 samples × 3 omics: 5–30 min (PCA/CCA); Autoencoder: 1–4 h |
| **Primary Use Cases** | Multi-omic biomarker discovery, cross-modal association, systems biology, data harmonization across platforms |

**Integration strategies**:
- **Early**: concatenate all features → single matrix; simple but high-dimensional
- **Intermediate**: joint factorization (PCA/CCA/NMF/AE) → shared latent space
- **Late**: train separate models per omic → ensemble predictions

**Example**:
```python
multi = MultiOmicsData()
multi.add_omics("genome", variants_df)
multi.add_omics("transcriptome", expr_df)
multi.add_omics("proteome", prot_df)
# Automatic sample alignment
joint = joint_pca(multi, n_components=20)
```

---

### Module 10: networks

| Attribute | Details |
|-----------|---------|
| **Category** | Biological Networks |
| **Data Types** | Edge list (TSV), adjacency matrix, graphML, BioPAX |
| **Scale** | Nodes: 10² → 10⁵+; Edges: 10² → 10⁶+ |
| **Key Methods** | Graph construction (PPI, co-expression, regulatory, metabolic), community detection (Louvain, Leiden, Girvan-Newman), centrality measures (degree, betweenness, closeness, eigenvector), network statistics (diameter, clustering coefficient), pathway enrichment on network modules, graph neural networks (planned), shortest paths, subgraph extraction |
| **Input** | Edgelist (node1, node2, [weight]), adjacency matrix, network file (graphML), optional node attributes |
| **Output** | NetworkX graph objects, community assignments, centrality scores, enriched pathways per community |
| **Submodules** | `analysis` (graph algorithms), `interaction` (PPI construction), `regulatory` (TF-target networks), `pathway` (KEGG/Reactome import), `visualization` (network plots) |
| **External Deps** | networkx, igraph (optional), cdlib (community detection), leidenalg |
| **Typical Runtime** | Small network (1K nodes): secs–min; Large network (100K nodes): hours |
| **Primary Use Cases** | Protein-protein interaction analysis, gene regulatory networks, functional module detection, network-based gene prioritization |

**Plots available**: spring layout, circular layout, hierarchical layout, edge-weighted graphs.

---

### Module 11: ml

| Attribute | Details |
|-----------|---------|
| **Category** | Machine Learning |
| **Data Types** | Feature matrix X (n × p), labels y (classification/regression) |
| **Scale** | Samples: 10² → 10⁶; Features: 10² → 10⁵ |
| **Key Methods** | Classification (random forest, SVM, logistic regression, naive Bayes), regression (linear, ridge, lasso, random forest), feature selection (mutual information, recursive elimination RFE, lasso, variance threshold), model evaluation (cross-validation, ROC/AUC, confusion matrix), interpretability (SHAP, feature importance), AutoML (hyperparameter optimization), LLM integration (Ollama for text/sequence tasks) |
| **Input** | X (NumPy/Pandas DataFrame/CSV), y (labels), optional feature names |
| **Output** | Trained model object, predictions, feature importance/coefficients, CV scores, classification/regression metrics |
| **Submodules** | `models` (classifiers/regressors), `features` (feature extraction/selection), `evaluation` (metrics, CV), `interpret` (SHAP, LIME), `automl` (Optuna-based), `llm` (Ollama integration) |
| **External Deps** | scikit-learn, xgboost, lightgbm, optuna, shap, ollama |
| **Typical Runtime** | RF on 10K × 1K: 1–10 min; Deep learning: hours; AutoML: days |
| **Primary Use Cases** | Predict phenotypes from genotypes/omics, feature discovery, biomarker identification, classification of disease states |

**Example**:
```python
from metainformant.ml import classification, features

# Feature selection
X_selected, mask = features.select_features(X, y, method="mutual_info", k=100)

# Train classifier
clf = classification.train_classifier(X_selected, y, method="rf")

# Cross-validate
cv_results = classification.cross_validate_classifier(clf, X_selected, y, cv=5)
```

---

### Module 12: math

| Attribute | Details |
|-----------|---------|
| **Category** | Mathematical Biology |
| **Data Types** | Parameters, equations, simulation parameters |
| **Scale** | Theoretical (symbolic/numeric) |
| **Key Methods** | Population genetics: Hardy-Weinberg equilibrium, inbreeding coefficient, F-statistics (Fst, Fis, Fit), genetic drift simulations, coalescent theory (time to MRCA), selection coefficients, selection tests; Epidemiology: SIR/SEIR compartmental models, R₀ estimation, epidemic curves; Evolutionary dynamics: fitness landscapes, replicator dynamics, game theory (ESS); Bayesian inference: MCMC, HMC (via stan/pymc) |
| **Input** | Parameters (pop size, mutation rate, migration rate, etc.), equations, or empirical data (allele frequencies) |
| **Output** | Analytical solutions, simulated datasets (allele frequencies over time), parameter estimates, confidence intervals |
| **Submodules** | `population_genetics` (HWE, drift, Fst, coalescent), `epidemiology` (SIR/SEIR), `evolutionary_dynamics` (game theory, fitness), `bayesian` (MCMC, HMC) |
| **External Deps** | numpy, scipy, optional: pymc, stan |
| **Typical Runtime** | Simulations: seconds–hours depending on complexity |
| **Primary Use Cases** | Theoretical modeling, null hypothesis generation, parameter exploration, educational demos, method validation (synthetic data) |

**Note**: Provides theoretical foundation for `dna.population` methods.

---

### Module 13: information

| Attribute | Details |
|-----------|---------|
| **Category** | Information Theory |
| **Data Types** | Discrete or continuous variables (vectors, matrices) |
| **Scale** | Features: 10² → 10⁴ |
| **Key Methods** | Shannon entropy (marginal, conditional), Rényi entropy (order α), Tsallis entropy (parameter q), mutual information (MI), conditional mutual information (CMI), normalized MI, Kullback-Leibler divergence, Jensen-Shannon divergence, transfer entropy (time series), semantic information content, network information flow (betweenness-like) |
| **Input** | Discrete/continuous variable vectors, joint probability tables, sequences |
| **Output** | Entropy values (bits/nats), MI matrices, normalized similarity scores, information profiles |
| **Submodules** | `metrics` (core measures), `integration` (cross-omics information), `workflow` (batch pipelines) |
| **External Deps** | numpy, scipy |
| **Typical Runtime** | O(n) for entropy; O(n²) for MI matrix of n features (use kNN estimator for continuous) |
| **Primary Use Cases** | Feature selection (non-linear dependence), sequence complexity, functional site detection, redundancy analysis, semantic similarity between GO terms |

**Applications**:
- `information.syntactic`: MI for feature selection in `ml.features`
- `information.genomic`: Sequence entropy for conserved regions
- `information.integration`: Cross-omics redundancy (which omic is most informative?)

---

### Module 14: ontology

| Attribute | Details |
|-----------|---------|
| **Category** | Functional Annotation & Enrichment |
| **Data Types** | Gene lists, OBO files (Gene Ontology), pathway databases |
| **Scale** | Genes: 10 → 10⁴ per analysis |
| **Key Methods** | GO graph parsing (OBO format), semantic similarity measures (Resnik, Lin, Jiang-Conrath), functional enrichment testing (hypergeometric/Fisher's exact), multiple testing correction (FDR), pathway annotation (KEGG/Reactome), gene set enrichment (GSEA-like) |
| **Input** | Gene list (Entrez/Ensembl IDs), optional background gene universe, GO OBO file (downloaded once) |
| **Output** | Enriched GO terms/ pathways (TSV with p_adj, enrichment ratio), similarity matrix (gene-gene or term-term), annotated gene lists per term |
| **Submodules** | `core` (OBO parser, graph), `query` (enrichment, similarity), `pathway_enrichment` (KEGG, Reactome) |
| **External Deps** | None (OBO download via HTTP); optional gseapy |
| **Typical Runtime** | Enrichment for 100 genes vs 20K background: <1 min |
| **Primary Use Cases** | GO enrichment after GWAS/eQTL/clustering, semantic similarity for gene clustering, pathway-level interpretation of omics results |

**Common pipelines**:
1. `gwas` lead genes → `ontology.enrichment` → significant GO terms
2. Cluster markers (singlecell) → enrichment → functional interpretation
3. Gene pairs (co-expression) → semantic similarity → functional linkage

---

### Module 15: phenotype

| Attribute | Details |
|-----------|---------|
| **Category** | Phenotypic Trait Analysis |
| **Data Types** | Measurement tables (CSV/TSV), images (future), tracking CSV (GPS/RFID) |
| **Scale** | Specimens: 10 → 10,000 |
| **Key Methods** | Morphometrics (linear measurements, allometric regression, shape indices, geometric morphometrics), behavioral analysis (ethograms, sequence analysis, transition matrices, time budgets), chemical profiling (GC-MS peaks, compound identification, distance matrices, marker detection), acoustic analysis (FFT spectral, syllable detection, call rate), electronic tracking (RFID, video tracking, movement ecology metrics, trajectory analysis), workflow orchestration (multi-step phenotype pipelines), visualization (trait correlation heatmaps, PCA, time series) |
| **Input** | Measurement table (`specimen_id, measurement_type, value, unit`), metadata (species, colony, date), optional raw data files (audio, video, GPS) |
| **Output** | Measurement profiles, derived indices (CI, wing loading, etc.), ethograms, transition matrices, chemical distance matrices, acoustic features, tracking trajectories, pipeline results (JSON) |
| **Submodules** | `morphological`, `behavior`, `chemical`, `electronic` (RFID/video/GPS), `sonic` (acoustic), `workflow`, `visualization`, `integration` (genotype mapping), `gwas_integration` (PheWAS, PRS), `analysis` (life course), `data` (I/O) |
| **External Deps** | None (pure Python); image/audio libs optional (PIL, librosa) |
| **Typical Runtime** | Small dataset (100 specimens × 10 traits): <1 min; Large (10K specimens): minutes–hours |
| **Primary Use Cases** | Comprehensive phenotyping for genetic studies, morphological comparisons, behavioral quantification, chemical ecology (CHC, pheromones), movement ecology, phenotype-genotype association (feed into `gwas` or `multiomics`) |

**Integration**: `phenotype.gwas_integration` provides PheWAS wrappers; `phenotype.integration` links to `multiomics`.

---

### Module 16: ecology

| Attribute | Details |
|-----------|---------|
| **Category** | Community Ecology |
| **Data Types** | Species abundance tables (samples × species), environmental variables, phylogenetic trees |
| **Scale** | Sites: 10 → 1,000; Species: 10 → 10,000 |
| **Key Methods** | Alpha diversity (species richness, Shannon, Simpson, phylogenetic diversity PD), beta diversity (Bray-Curtis, Jaccard, UniFrac), ordination (PCA, NMDS, PCoA, CCA, RDA), species-environment correlation (envfit, permutation tests), rarefaction curves, diversity partitioning, community similarity clustering |
| **Input** | Abundance matrix (samples × species/species), environmental variables table (samples × env), optional phylo tree for PD |
| **Output** | Diversity indices (vector), ordination coordinates (samples × axes), environmental fit vectors, rarefaction curves, beta diversity distance matrices |
| **Submodules** | `analysis` (diversity, ordination, envfit), `phylogenetic` ( Faith's PD, UniFrac), `visualization` (ordination plots, diversity barplots) |
| **External Deps** | scikit-bio (optional), pandas, numpy |
| **Typical Runtime** | 100 samples × 500 species: <5 min |
| **Primary Use Cases** | Biodiversity assessment, community composition comparison, environmental drivers of community structure, phylogenetic diversity analysis |

**Example**:
```python
from metainformant.ecology.analysis import shannon_index, bray_curtis, nmds
alpha = shannon_index(abundance_matrix)
beta_dm = bray_curtis(abundance_matrix)
nmds_coords = nmds(beta_dm, n_components=2)
```

---

### Module 17: simulation

| Attribute | Details |
|-----------|---------|
| **Category** | Synthetic Data Generation |
| **Data Types** | Parameters (pop size, mutation rate, selection coefficient); output: sequences, genotypes, expression |
| **Scale** | Configurable (small validation → large benchmark) |
| **Key Methods** | Sequence evolution (substitution models: Jukes-Cantor, Kimura, GTR), coalescent simulations (ms, msprime wrapper), population genetics sims (selection, drift, migration), agent-based models (ecosystem/phenotype agents), benchmark dataset generation (ground truth), null model generation |
| **Input** | Simulation parameters (model, length, N, μ, etc.), random seed, model specifications |
| **Output** | Simulated FASTA, VCF, count matrices, agent trajectories, metadata |
| **Submodules** | `models` (evolutionary, population, agent-based), `workflow` (pipeline runners), `benchmark` (synthetic benchmarks for testing) |
| **External Deps** | msprime (coalescent), optional: SLiM (forward-time) |
| **Typical Runtime** | Coalescent (10K sequences, 10K bp): 1–10 min; Agent-based (10K agents, 1000 steps): minutes–hours |
| **Primary Use Cases** | Method validation (truth-known), power analysis, null hypothesis generation, benchmark dataset creation, educational demos |

**Example**: Simulate neutral dataset to test Tajima's D → should be ~0:
```python
from metainformant.simulation.models.coalescent import simulate_sequences
seqs = simulate_sequences(n=100, length=1000, theta=10)
from metainformant.dna.population import tajimas_d
td = tajimas_d(seqs)  # Expect ~0 if neutral
```

---

### Module 18: quality

| Attribute | Details |
|-----------|---------|
| **Category** | Data Quality Assessment |
| **Data Types** | FASTQ, BAM, assemblies (FASTA), count matrices |
| **Scale** | Reads: 10⁴ → 10⁸+ |
| **Key Methods** | Quality metrics: per-base Phred scores, GC content, adapter content, sequence length distribution, duplication levels, overrepresented sequences; contamination detection (kmer-based vs. reference), assembly validation (BUSCO completeness, N50, L50), batch effect diagnostics (PCA/UMAP colored by batch), coverage depth histograms |
| **Input** | FASTQ (raw reads), BAM (aligned), assembled FASTA, count matrix |
| **Input** | FASTQ, BAM, FASTA (assembly), count matrix |
| **Output** | QC report (JSON/HTML), quality plots (PNG/SVG), summary statistics table |
| **Submodules** | `analysis` (metrics computation), `io` (FASTQ/BAM parsers), `reporting` (JSON/HTML) |
| **External Deps** | fastqc, samtools, quast/busco (optional wrappers) |
| **Typical Runtime** | Single FASTQ (1G): 2–5 min; 100 FASTQs: parallelized 1–2 h |
| **Primary Use Cases** | Sequencing quality control, contamination screening, assembly validation, batch effect detection, preprocessing validation |

**Key outputs**: quality_scores.png, adapter_content.png, gc_distribution.png, qc_summary.json

---

### Module 19: visualization

| Attribute | Details |
|-----------|---------|
| **Category** | Plotting & Graphics |
| **Data Types** | Any plot-ready data (arrays, DataFrames, specialized objects) |
| **Scale** | Points: 10² → 10⁷+; Genes: 10² → 10⁵ |
| **Key Methods** | 70+ plot types organized into submodules (see below); static (Matplotlib/Seaborn) and interactive (Plotly) backends; publication styling (Set2 palette, Helvetica, 300+ DPI); multi-panel dashboard composition; animation engine |
| **Input** | Plot-specific: arrays for line/scatter; matrices for heatmaps; trees (Newick); graphs (NetworkX); embeddings (n × 2/3) |
| **Output** | PNG/TIFF/PDF/SVG (static), HTML (interactive), MP4/GIF (animation) |
| **Submodules** | `basic/` (12 plots: line, scatter, bar, pie, area, step, heatmap, etc.), `statistical/` (hist, box, violin, qq, roc, density), `genomics/` (manhattan, volcano, regional, ideogram, circos), `expression/` (expr heatmap, MA, enrichment), `dimred/` (PCA, UMAP, t-SNE, scree), `networks/` (network graphs, circular, hierarchical), `trees/` (rectangular, circular, radial trees), `timeseries/` (line+CI, forecast, ACF, seasonal), `multidim/` (pairplot, parallel coords, radar), `quality/` (FASTQ QC plots), `information/` (entropy, MI heatmap), `animations/` (6 animation types), `dashboards/` (multi-panel layout), `styling/` (themes, palettes) |
| **External Deps** | matplotlib, seaborn, plotly, networkx, Bio.Phylo |
| **Typical Runtime** | Simple plot (<10⁴ points): <1 s; Large heatmap (10⁶ cells, sparse): 5–30 s; Animation (100 frames): 1–5 min |
| **Primary Use Cases** | Publication figures, exploratory data analysis, interactive dashboards, animated demos, all-domain visualization needs |

**See detailed comparison**: [visualization_approaches.md](visualization_approaches.md)

---

### Module 20: longread

| Attribute | Details |
|-----------|---------|
| **Category** | Long-Read Sequencing |
| **Data Types** | FAST5 (ONT), POD5 (ONT), BAM (subreads), FASTA/Q |
| **Scale** | Reads: 1K → 1M+ |
| **Key Methods** | FAST5/POD5 I/O and signal extraction, basecalling integration (ont-guppy/dorado), read quality metrics (N50, accuracy, throughput), modified base detection (5mC, 6mA methylation calling), structural variant detection (split-read, read-pair, depth), haplotype phasing (WhatsHap-style), de novo assembly (overlap-layout-consensus, hifiasm-like), long-read-specific visualization (dotplots, methylation tracks) |
| **Input** | Raw FAST5/POD5, subread BAM, basecalled FASTQ |
| **Output** | Basecalled FASTQ, methylation calls (BED), SV VCF, phased haplotypes, assembled contigs (FASTA), assembly graphs |
| **Submodules** | `io` (format converters), `quality` (N50, accuracy), `analysis` (SVs, methylation, phasing), `assembly` (OLC, consensus), `visualization` (long-read plots), `workflow` (end-to-end) |
| **External Deps** | `guppy`/`dorado` (ONT basecaller), `minimap2`, `samtools`, `sniffles2`, `pbsv` |
| **Typical Runtime** | Basecall 1M reads (Guppy): 1–4 h; SV calling on 30X human: 4–12 h; Assembly (hi-fi): 1–3 days |
| **Primary Use Cases** | Structural variant detection (large indels, translocations), methylation profiling without bisulfite, haplotype-resolved assembly, full-length isoform sequencing (Iso-Seq) |

**Platforms**: Oxford Nanopore (ONT), PacBio HiFi/CLR; both raw signal and basecalled formats.

---

### Module 21: metagenomics

| Attribute | Details |
|-----------|---------|
| **Category** | Metagenomics & Microbiome |
| **Data Types** | FASTQ (amplicon 16S/ITS, shotgun), SRA runs |
| **Scale** | Samples: 10 → 10,000+ |
| **Key Methods** | Amplicon: denoising (DADA2, Deblur), OTU clustering (UCLUST/VSEARCH), ASV generation, taxonomic classification (BLAST, naive Bayes), tree building (FastTree); Shotgun: assembly (metaSPAdes, MEGAHIT), binning (MetaBAT2, CONCOCT, MaxBin), gene prediction (Prodigal), functional annotation (eggNOG-mapper, DIAMOND, HMMER), metabolic pathway reconstruction (KEGG modules), community profiling (relative abundance, alpha/beta diversity) |
| **Input** | 16S/ITS FASTQ (amplicon), shotgun metagenome FASTQ, reference databases (SILVA/Greengenes for 16S, eggNOG/KEGG for functions) |
| **Output** | ASV/OTU table (samples × features), taxonomic assignments, MAGs (Metagenome-Assembled Genomes) in FASTA, functional annotation tables (KO, COG), pathway completeness charts |
| **Submodules** | `amplicon` (16S/ITS pipeline), `shotgun` (assembly, binning, profiling), `functional` (annotation, pathways), `visualization` (Krona, ordination, rarefaction) |
| **External Deps** | QIIME2 (amplicon), metaSPAdes/MEGAHIT (assembly), MetaBAT2/CONCOCT (binning), DIAMOND, HMMER |
| **Typical Runtime** | 16S 100 samples (DADA2): 2–6 h; Shotgun assembly (10G sample): 6–24 h |
| **Primary Use Cases** | 16S/ITS microbial community profiling, metagenomic assembly and binning, functional metagenomics (KEGG pathways), microbiome-disease association, strain-level analysis |

**Quality control**: `quality` module for per-sample FastQC before proceeding.

---

### Module 22: structural_variants

| Attribute | Details |
|-----------|---------|
| **Category** | Structural Variation |
| **Data Types** | BAM (aligned reads), VCF (existing SV calls) |
| **Scale** | Samples: 10 → 10,000+ |
| **Key Methods** | CNV detection: read depth analysis, segmentation (CBS, HMM), CNVnator-style analysis; SV calling: split-read (Pindel, Manta), discordant pair (LUMPY), read-depth (CNV), de novo assembly-based (SPAdes); breakpoint refinement (microhomology, insert size); gene/regulatory overlap annotation (gene fusion, exon deletion, enhancer disruption); functional impact scoring (pathogenicity, dosage sensitivity scores); TAD disruption (CTCF site breakage, loop loss); quality filtering (blacklist regions, depth threshold, supporting read count); multi-caller merging (reciprocal overlap, SURVIVOR); visualization (Circos plots, coverage plots) |
| **Input** | Aligned BAM (or CRAM), reference genome, gene annotation GTF, optional existing SV VCFs to merge |
| **Output** | SV VCF (structural variants: DEL, DUP, INV, TRA, BND), CNV BED (copy number segments), annotated regions (overlap with genes/regulatory), quality-filtered VCF, Circos plot (SV karyogram) |
| **Submodules** | `detection` (CNV, SV calling, breakpoints), `annotation` (gene/reg overlap, functional impact), `filtering` (quality, blacklist, merge), `visualization` (Circos, coverage plots) |
| **External Deps** | `samtools`, `bedtools`, `survivor`, `manta`, `lumpy`, `cnvnator` (wrappers) |
| **Typical Runtime** | Single WGS BAM (30X): Manta 2–4 h; CNV detection: 1–2 h; Multi-caller merge: minutes |
| **Primary Use Cases** | Structural variant detection (indels >50bp, inversions, translocations, CNVs), clinical cytogenetics, cancer genome analysis, gene fusion detection, non-coding regulatory disruption |

**Note**: Often downstream of `longread` (SV from long reads) or used with `gwas` (SV-based association).

---

### Module 23: pharmacogenomics

| Attribute | Details |
|-----------|---------|
| **Category** | Clinical Pharmacogenomics |
| **Data Types** | VCF (genotypes), diplotype tables, clinical annotations |
| **Scale** | Patients: 10 → 10,000+ |
| **Key Methods** | Star allele determination (CPIC-defined allele definitions per gene), diplotype construction (phasing via pedigree or statistical), metabolizer phenotype prediction (NORMAL/INTERMEDIATE/POOR/ULTRA-RAPID for CYP enzymes), CPIC guideline integration (drug-gene pairs, dosing recommendations by phenotype), PharmGKB clinical annotations (evidence levels, PMC IDs), FDA drug label parsing (Table of Pharmacogenomic Biomarkers), ACMG/AMP variant classification (pathogenicity criteria), drug-gene interaction analysis (polypharmacy checks), clinical report generation (with disclaimers, confidence) |
| **Input** | Patient genotype VCF (CYP2D6, CYP2C19, VKORC1, etc.), optional clinical context (drugs prescribed, indications) |
| **Output** | Clinical report (PDF/HTML) with: star alleles, diplotype, metabolizer phenotype, drug recommendations (dose adjustments, alternatives), evidence citations (CPIC level A/B, PharmGKB), variant pathogenicity classification |
| **Submodules** | `alleles` (star allele, diplotype, phenotype), `annotations` (CPIC, PharmGKB, FDA), `clinical` (interactions, reporting, ACMG scoring), `visualization` (clinical plots) |
| **External Deps** | ClinVar API (variant pathogenicity), gnomAD (population freq), CPIC guidelines DB, PharmGKB API |
| **Typical Runtime** | Single patient VCF parse + annotation: seconds–minutes |
| **Primary Use Cases** | Clinical pharmacogenomic testing, drug dosing guidance, genotype-driven prescribing, drug-gene interaction checking, clinical decision support |

**Genes covered**: CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, NUDT15, etc.

**Guideline levels**: CPIC level A (strong) → actionable; level B (moderate) → may affect.

**Example**:
```python
from metainformant.pharmacogenomics import alleles, annotations
diplotype = alleles.determine_diplotype(vcf, gene="CYP2D6")
phenotype = alleles.predict_metabolizer(diplotype)
guidelines = annotations.get_cpic_guidelines(drug="codeine", phenotype=phenotype)
```

---

### Module 24: metabolomics

| Attribute | Details |
|-----------|---------|
| **Category** | Metabolomics |
| **Data Types** | mzML, mzXML, CSV (peak tables), metabolites (HMDB/KEGG IDs) |
| **Scale** | Samples: 10 → 1,000; Metabolites: 50 → 10,000 |
| **Key Methods** | MS data parsing (mzML/mzXML), peak detection (centroiding), metabolite identification (mz matching to HMDB/KEGG/MassBank), quantification (peak area/height), data normalization (PQN, log, probabilistic quotient), missing value imputation (kNN, half-min), batch effect correction (ComBat), batch effect detection (PCA), KEGG/Reactome pathway mapping, metabolite set enrichment analysis (MSEA, hypergeometric), metabolite-gene integration (via `multiomics`), visualization (volcano, PCA, heatmap) |
| **Input** | mzML files (raw MS), sample metadata, reference metabolite database (HMDB/KEGG IDs + m/z), optional quantification table (CSV) |
| **Output** | Peak table (samples × metabolites), identified metabolites (with HMDB/KEGG IDs), normalized abundance matrix, enriched pathways (KEGG modules), integrated scores (metabolite–gene association if multi-omics) |
| **Submodules** | `io` (mzML parser), `analysis` (peak detection, ID, quantification, normalization), `pathways` (KEGG/Reactome mapping, MSEA), `visualization` (metabolomics-specific plots), `integration` (with multiomics) |
| **External Deps** | msconvert (ProteoWizard), matchpy (mz matching), optional: MZmine/XCMS (external) |
| **Typical Runtime** | 100 samples × 5K metabolites: 30 min–2 h |
| **Primary Use Cases** | Mass spectrometry-based metabolomics profiling, biomarker discovery, metabolic pathway analysis, metabolite-gene integration |

**Normalization methods**: total ion count (TIC), probabilistic quotient (PQN), batch correction (ComBat).

**Pathway mapping**: map metabolite → KEGG compound → pathway (hsa00010 for glycolysis, etc.).

---

### Module 25: menu

| Attribute | Details |
|-----------|---------|
| **Category** | Interactive Discovery |
| **Data Types** | Scripts metadata (any) |
| **Scale** | Scripts: 10 → 1,000+ |
| **Key Methods** | Script discovery (scan `scripts/` directory), metadata extraction (docstrings, argparse help), categorization by domain/pipeline, interactive terminal menu with breadcrumb navigation, script execution with argument prompting/validation, dynamic menu data structures (MenuItem, Menu, MenuSystem, MenuHistory) |
| **Input** | Discovered scripts from `scripts/` with docstring metadata |
| **Output** | Interactive terminal UI; executed scripts with populated CLI arguments |
| **Submodules** | `core.discovery` (script scanner), `core.executor` (script runner), `ui.display` (terminal formatting), `ui.navigation` (menu state) |
| **External Deps** | None (terminal-based, uses `curses`/`prompt_toolkit`) |
| **Typical Runtime** | Instant response (menu-driven) |
| **Primary Use Cases** | Interactive exploration of available scripts, guided workflow execution, parameter prompting, menu-driven analysis launch |

**Example**:
```bash
# Launch interactive menu
uv run python -m metainformant menu

# Navigate: "RNA-seq" → "Amalgkit workflow" → prompts for species/Fastq dir
# Then launch script with populated arguments automatically.
```

---

### Module 26: cloud

| Attribute | Details |
|-----------|---------|
| **Category** | Cloud Infrastructure |
| **Data Types** | GCP resources (VMs, disks, buckets) |
| **Scale** | VMs: 1 → 10,000+ |
| **Key Methods** | VM lifecycle (create, delete, wait for SSH), Docker build & run, file transfer (SCP, `gsutil` rsync), cost estimation (per-VM per-hour), preemptible VM orchestration (70-80% savings), genome index caching (GCS `gs://metainformant-genomes/`), pipeline deployment (RNA amalgkit at scale), result collection (compress + download), status monitoring + streaming logs |
| **Input** | Cloud configuration YAML (project, zone, machine type, Docker image, pipeline), input data location (GCS/local) |
| **Output** | Running/completed VMs, output directories (local or GCS), cost logs, execution manifests |
| **Submodules** | `cloud_config` (32-field dataclass), `gcp_deployer` (VM orchestration class) |
| **External Deps** | `gcloud` CLI (thin wrapper), `docker`, `gsutil`, `ssh` |
| **Typical Runtime** | VM create: 1–2 min; Pipeline run: hours–days; Download results: 10 min–hours |
| **Primary Use Cases** | Large-scale genomic analysis at scale (e.g., 8,300 RNA-seq samples), cost-optimized compute (preemptible VMs), reproducible pipeline deployment, genome indexing at scale (gcp/prep_genomes.py) |

**Economics**: Preemptible VMs cost ~$0.01–0.03 per CPU-hour vs. normal $0.05–0.15.

**Key consumers**: `rna` module for ENA streaming workflows (orchestration), `gwas` for large cohorts.

---

### Module 27: eqtl

| Attribute | Details |
|-----------|---------|
| **Category** | eQTL Integration *(cross-cutting)* |
| **Data Types** | VCF (genotypes) + expression matrix (counts/TPM) |
| **Scale** | Samples: 100 → 10,000+ (expression cohort) |
| **Key Methods** | cis-eQTL scanning (window ±1 Mb around gene), trans-eQTL scanning (genome-wide, multiple testing extreme), colocalization ( coloc , eCAVIAR, enloc), mediation analysis, transcriptome SNP calling (from RNA BAM), eQTL mapping with Amalgkit data (cross-species TMM normalization), fine-mapping integration with GWAS, summary-statistic colocalization (fast) |
| **Input** | Genotype VCF (common variants), expression matrix (genes × samples, matched samples), optional GWAS summary stats (for colocalization) |
| **Output** | eQTL results (snp_id, gene_id, beta, p, qvalue), colocalization PP4 (posterior probability shared causal variant), mediated effect sizes, transcriptome VCF (from RNA BAM) |
| **Submodules** | Core logic in `gwas.finemapping.eqtl` and `multiomics.analysis.integration`; scripts in `scripts/eqtl/` |
| **External Deps** | `plink2` (geno QC), `matrixeqtl` (fast eQTL), `coloc` R package, `SUSIE` (fine-mapping) |
| **Typical Runtime** | cis-eQTL (10K samples, 20K genes, 1M SNPs): 1–3 hours; trans-eQTL (all pairs): days; coloc: minutes–hours |
| **Primary Use Cases** | Mapping genetic regulatory effects (cis and trans), colocalization between GWAS and eQTL, mechanistic interpretation of GWAS hits (which gene mediates?), cross-species eQTL (via Amalgkit normalization) |

**Note**: eQTL is cross-cutting, not a standalone module; scripts live in `scripts/eqtl/`, but core ETL and logic under `gwas/finemapping/`.

---

### Module 28: life_events

| Attribute | Details |
|-----------|---------|
| **Category** | Temporal Event Analysis |
| **Data Types** | Event sequences (ordered lists), time-to-event data, longitudinal measurements |
| **Scale** | Individuals: 100 → 10,000+; Events per person: 1 → 100+ |
| **Key Methods** | Event sequence embedding (categorical event → vector via word2vec/GloVe), temporal pattern discovery (clustering of event sequences), survival analysis (Cox PH, Kaplan-Meier), outcome prediction modeling (events → survival/risk), life course modeling (sequence analysis, optimal matching), time-varying covariates, recurrent events, competing risks, sequence visualization (calendar plots, state-transition diagrams) |
| **Input** | Event logs (subject_id, event_type, timestamp), survival data (time, event_observed, covariates), longitudinal measurements |
| **Output** | Event embeddings (subject × embedding_dim), survival curves, hazard ratios, predicted risk scores, clustered event sequence groups |
| **Submodules** | `core` (event data structures), `models` (survival, prediction), `workflow` (analysis pipelines) |
| **External Deps** | lifelines (survival), scikit-learn (clustering), optional: gensim (embeddings) |
| **Typical Runtime** | Cox PH (1K samples, 20 covariates): <1 min; embedding (10K sequences, 100 event types): 5–15 min |
| **Primary Use Cases** | Life history analysis (e.g., ant colony events), longitudinal health records, time-to-event outcomes, sequential pattern mining, event-based risk prediction |

**Example**: Modeling ant colony development from worker event logs; predicting colony fate from early events.

---

## Cross-Module Dependencies

| From → To | Dependency Type | Example Use |
|-----------|----------------|-------------|
| `rna` → `dna` | Translation (DNA→RNA), ortholog mapping | `dna.transcription.dna_to_rna()` |
| `eqtl` → `gwas` | Uses GWAS summary stats for coloc | `gwas.finemapping.coloc.coloc()` |
| `eqtl` → `multiomics` | Shared integration logic | `multiomics.analysis.integration` |
| `multiomics` → `gwas` | Variant data as DNA layer | `multiomics.add_omics("genome", gwas_variants)` |
| `multiomics` → `rna` | RNA expression as transcriptome layer | `multiomics.add_omics("transcriptome", tpm_df)` |
| `spatial` → `singlecell` | scRNA reference for deconvolution | `spatial.integration.map_sc_to_spatial()` |
| `visualization` → all | All domains use for plotting | `visualization.genomics.manhattan_plot()` |
| `quality` → all pre-processing | QC first | `quality.analysis.fastqc()` |

**Rules**:
- Domain modules (`dna`, `rna`, `gwas`, etc.) **do not import each other directly**; cross-domain coordination happens via `multiomics`, `eqtl`, or shared configs.
- All modules use `core` infrastructure (`core.io`, `core.parallel`, `core.utils.config`).

---

## Summary: How to Pick the Right Module

1. **Start with your data format**:
   - FASTA/VCF → `dna`
   - FASTQ (bulk RNA) → `rna`
   - h5ad/counts (single cells) → `singlecell`
   - Spatial coordinates + counts → `spatial`
   - Multiple omic tables → `multiomics`

2. **Consider your question**:
   - "Which variants associate?" → `gwas`
   - "What pathways enriched?" → `ontology`
   - "How diverse is community?" → `ecology`
   - "Predict from features?" → `ml`
   - "Visualize results?" → `visualization.{submodule}`

3. **Evaluate scale**:
   - Small (<1K samples): most modules OK
   - Medium (1K–100K): `rna`, `gwas`, `singlecell` OK with parallel
   - Large (>100K, population biobank): need `cloud`, `parallel`, or HPC

4. **Check integration needs**:
   - Single-omic: use domain module alone
   - Multi-omic: `multiomics` container
   - GWAS+RNA: `eqtl` pipeline
   - Multi-omic outcomes: `multiomics` + `ml`

---

## Related Documentation

- **[COMPARISON_GUIDES.md](../COMPARISON_GUIDES.md)** — Master index of all comparison guides
- **[dna vs RNA vs Transcriptome](dna_vs_rna_vs_transcriptome.md)** — Detailed domain comparison
- **[GWAS vs Phenotype vs Multi-omics](gwas_vs_phenotype_vs_multiomics.md)** — Study design & power comparison  
- **[Visualization Approaches](visualization_approaches.md)** — Plot selection guide
- **[docs/dna/](../dna/index.md)** — DNA module full documentation
- **[docs/rna/](../rna/index.md)** — RNA module full documentation
- **[Each module's `index.md` or `SPEC.md`]** — Module-specific specs

---

*Generated from module agents docs and SPEC files. For detailed per-module API, see individual module `API.md` or `index.md`.*
