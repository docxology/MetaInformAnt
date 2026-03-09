# Comprehensive RNA-seq Pipeline Guide

This guide details the METAINFORMANT RNA-seq pipeline, covering end-to-end methods from data retrieval to quantified expression matrices. This documentation is designed to be cited as a reproducible methods reference for publications.

## 1. Pipeline Overview

The METAINFORMANT RNA-seq pipeline automates large-scale transcriptomic analysis across 28 Hymenoptera species (21 ant and 7 bee/wasp species) using the [amalgkit](https://github.com/kfuku52/amalgkit) framework. The pipeline is implemented in Python 3.13+ and is managed via the `metainformant` package (`src/metainformant/rna/`).

### Key Features

* **Automated Data Retrieval**: Two-tier download from ENA (primary) and NCBI SRA (fallback)
* **Pseudo-alignment Quantification**: Kallisto-based transcript abundance estimation
* **Tissue Normalization**: Automatic harmonization of tissue/organ labels across studies
* **Index Complexity Management**: Automatic filtering of problematic transcripts for robust indexing
* **Streaming Architecture**: Download → quantify → delete per-sample for disk efficiency
* **SQLite Progress Tracking**: Resumable, concurrent-safe state management via WAL-mode database
* **Cross-species Analysis**: TMM normalization (CSTMM) and correlation analysis (CSCA) via orthogroup mappings

### Software Dependencies

| Tool | Version | Purpose |
|---|---|---|
| [amalgkit](https://github.com/kfuku52/amalgkit) | 0.16.0 | RNA-seq workflow framework (metadata, download, quant, merge, curate) |
| [Kallisto](https://pachterlab.github.io/kallisto/) | ≥0.50.0 | Pseudo-alignment and transcript quantification |
| [SRA Toolkit](https://github.com/ncbi/sra-tools) | ≥3.0 | NCBI SRA data access (fasterq-dump fallback) |
| [curl](https://curl.se/) | ≥7.0 | ENA direct HTTP/FTP downloads |
| Python | ≥3.13 | Pipeline orchestration |
| SQLite | ≥3.35 | Progress database (WAL mode) |

## 2. Data Retrieval

### Sample Discovery

For each species, publicly available RNA-seq datasets are discovered by querying the NCBI Sequence Read Archive (SRA) with species-specific Entrez search strings:

```
"<Species name>"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]
```

This retrieves all Illumina-based RNA-seq runs for the target species. Metadata including tissue type, BioProject, sample accession, and estimated file size are extracted via the `amalgkit metadata` step.

### Sample Selection and Filtering

Samples are filtered by:
- **Platform**: Illumina sequencing only
- **Strategy**: RNA-Seq experiments only
- **File size**: Configurable maximum (default: 20 GB per sample) to skip excessively large files
- **Tissue requirement**: Optional tissue annotation filtering (configurable per species)

### Download Strategy: ENA-First with NCBI Fallback

All species use a two-tier download strategy managed by `StreamingPipelineOrchestrator`:

1. **ENA primary** — The European Nucleotide Archive Portal API is queried for direct HTTPS links to compressed FASTQ files (`.fastq.gz`). Downloads use `curl` with retry logic and MD5 verification. This bypasses the slow `prefetch` + `fasterq-dump` SRA extraction pipeline.

2. **NCBI fallback** — If ENA download fails (no ENA record, network error, or corrupt file), the pipeline falls back to `fasterq-dump` from the NCBI SRA Toolkit for direct FASTQ extraction.

**Concurrency**: Up to N parallel download workers (configurable, typically 12–24) operate simultaneously, with a single shared SQLite progress database tracking sample states (`pending` → `downloading` → `quantifying` → `quantified` or `failed`).

## 3. Reference Genome Preparation

### Genome Download

Reference genomes are downloaded from NCBI RefSeq using the assembly accession specified in each species configuration file. The `GenomePreparator` class (`src/metainformant/rna/amalgkit/genome_prep.py`) handles:

1. Downloading the `_rna_from_genomic.fna.gz` transcriptome file from NCBI FTP
2. Validating file integrity
3. Building the Kallisto transcriptome index

### Index Complexity Management

For genomes with high repetitive content (notably *Harpegnathos saltator*), standard Kallisto indexing can produce indices with excessive equivalence class sizes (Max EC > 3,000), causing `kallisto quant` to hang.

The `IndexComplexityManager` (`src/metainformant/rna/amalgkit/index_prep.py`) automatically:
1. Filters out `XR_` and `NR_` (non-coding RNA) transcript accessions
2. Removes transcripts shorter than 200 bp
3. Deduplicates identical sequences
4. Rebuilds the Kallisto index with reduced complexity

This is applied automatically for all species.

## 4. Quantification

Transcript-level abundance is estimated using **Kallisto** pseudo-alignment (Bray et al., 2016). For each sample:

1. Compressed FASTQ files are streamed to Kallisto with automatic single-end/paired-end detection
2. The `amalgkit quant --batch` mode executes Kallisto with the species-specific transcriptome index
3. Output: `abundance.tsv` per sample containing estimated counts, TPM, and effective lengths per transcript

### Per-Sample Processing Loop

The streaming orchestrator implements an immediate cleanup pattern:
```
For each sample:
  1. Download FASTQ from ENA (or NCBI fallback)
  2. Quantify against Kallisto index → abundance.tsv
  3. Delete FASTQ files immediately
  4. Update progress database
```

This ensures only 1–N samples' FASTQ files exist at any time, enabling processing of datasets much larger than available disk space.

## 5. Tissue Normalization

Raw tissue annotations from SRA metadata are heterogeneous (e.g., "whole body", "whole_body", "Whole Body", "whole animal"). The `TissueNormalizer` (`src/metainformant/rna/amalgkit/tissue_normalizer.py`) harmonizes these using:

1. **Mapping file** (`config/amalgkit/tissue_mapping.yaml`): A manually curated dictionary of >200 tissue label mappings organized hierarchically
2. **Patch file** (`config/amalgkit/tissue_patches.yaml`): Sample-level, BioProject-level, and BioSample-level overrides for cases not captured by the general mapping

Normalization is applied destructively to the `tissue` column in the metadata TSV before caching, ensuring downstream `cstmm` and `csca` analyses operate on consistent tissue labels.

## 6. Expression Matrix Merging

After quantification, per-sample `abundance.tsv` files are merged into a single expression matrix using `amalgkit merge`. The merged matrix (`merged_abundance.tsv`) contains:
- Rows: transcript IDs
- Columns: sample accessions (SRR IDs)
- Values: estimated counts (or TPM, configurable)

## 7. Quality Curation

The `amalgkit curate` step applies quality control and filtering:
- Outlier detection via inter-sample correlation
- Library size normalization assessment
- Bias correction for batch effects
- Output: curated expression matrix with QC metrics

## 8. Cross-Species Analysis

### CSTMM (Cross-Species TMM Normalization)

For comparative analyses, expression matrices from multiple species are normalized together using TMM (Robinson & Oshlack, 2010) via orthogroup mappings. Requires:
- `Orthogroups.tsv` mapping file (from OrthoFinder or equivalent)
- Per-species merged expression matrices

### CSCA (Cross-Species Correlation Analysis)

Generates cross-species PCA and correlation analyses. Operates on CSTMM-normalized matrices, grouping by species and tissue with orthogroup-level feature alignment.

## 9. Orchestration Architecture

### Entry Points

| Script | Purpose | Use Case |
|---|---|---|
| `scripts/rna/run_all_species.py` | Multi-species sequential orchestrator | Production runs (all 28 species) |
| `scripts/rna/run_workflow.py` | Single-species full workflow | Individual species / debugging |
| `scripts/rna/check_pipeline_status.py` | Status dashboard | Monitoring active runs |

### StreamingPipelineOrchestrator

The core orchestrator (`src/metainformant/rna/engine/streaming_orchestrator.py`) implements:

- **Species iteration**: Processes species from a prioritized list (smallest first by default)
- **Per-species workflow**: Metadata → tissue normalization → genome index verification → concurrent download+quant
- **Concurrent workers**: ThreadPoolExecutor with configurable worker count
- **Progress database**: SQLite with WAL mode for concurrent read/write safety
- **Automatic reconciliation**: On startup, scans filesystem for already-quantified samples and updates the database

### Progress Database Schema

```sql
CREATE TABLE samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    species TEXT NOT NULL,
    srr_id TEXT NOT NULL UNIQUE,
    state TEXT NOT NULL DEFAULT 'pending',
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
-- States: pending → downloading → quantifying → quantified | failed
```

### Cloud Deployment

For massive datasets, the pipeline runs inside a Docker container on Google Cloud Platform:
- **VM**: n2-standard-32 (32 vCPUs, 128 GB RAM, 1 TB SSD)
- **Management**: `scripts/cloud/deploy_gcp.py` CLI for VM lifecycle
- **Container**: `metainformant-pipeline` Docker image with all bioinformatics tools pre-installed
- **Volume mount**: `/opt/MetaInformAnt` → container `/app`

## 10. Species Coverage

The pipeline is configured for 28 Hymenoptera species across two major clades:

### Formicidae (Ants) — 21 species

| Species | Config | Assembly |
|---|---|---|
| *Acromyrmex echinatior* | `amalgkit_acromyrmex_echinatior.yaml` | GCF_000204515.1 |
| *Anoplolepis gracilipes* | `amalgkit_anoplolepis_gracilipes.yaml` | GCF_047496725.1 |
| *Atta cephalotes* | `amalgkit_atta_cephalotes.yaml` | GCF_000143395.1 |
| *Camponotus floridanus* | `amalgkit_camponotus_floridanus.yaml` | GCF_003227725.1 |
| *Cardiocondyla obscurior* | `amalgkit_cardiocondyla_obscurior.yaml` | GCF_019399895.1 |
| *Dinoponera quadriceps* | `amalgkit_dinoponera_quadriceps.yaml` | GCF_001313825.1 |
| *Formica exsecta* | `amalgkit_formica_exsecta.yaml` | GCF_003651465.1 |
| *Harpegnathos saltator* | `amalgkit_harpegnathos_saltator.yaml` | GCF_003227715.2 |
| *Linepithema humile* | `amalgkit_linepithema_humile.yaml` | GCF_040581485.1 |
| *Monomorium pharaonis* | `amalgkit_monomorium_pharaonis.yaml` | GCF_013373865.1 |
| *Nylanderia fulva* | `amalgkit_nylanderia_fulva.yaml` | GCF_005281655.2 |
| *Odontomachus brunneus* | `amalgkit_odontomachus_brunneus.yaml` | GCF_010583005.1 |
| *Ooceraea biroi* | `amalgkit_ooceraea_biroi.yaml` | GCF_003672135.1 |
| *Pogonomyrmex barbatus* | `amalgkit_pogonomyrmex_barbatus.yaml` | GCF_000187915.1 |
| *Solenopsis invicta* | `amalgkit_solenopsis_invicta.yaml` | GCF_016802725.1 |
| *Temnothorax americanus* | `amalgkit_temnothorax_americanus.yaml` | GCF_048541705.1 |
| *Temnothorax curvispinosus* | `amalgkit_temnothorax_curvispinosus.yaml` | GCF_003070985.1 |
| *Temnothorax longispinosus* | `amalgkit_temnothorax_longispinosus.yaml` | GCF_030848805.1 |
| *Temnothorax nylanderi* | `amalgkit_temnothorax_nylanderi.yaml` | GCF_030848795.1 |
| *Vollenhovia emeryi* | `amalgkit_vollenhovia_emeryi.yaml` | GCF_000949405.1 |
| *Wasmannia auropunctata* | `amalgkit_wasmannia_auropunctata.yaml` | GCF_000956235.1 |

### Other Hymenoptera (Bees/Wasps) — 7 species

| Species | Config | Assembly |
|---|---|---|
| *Apis mellifera* | `amalgkit_apis_mellifera.yaml` | GCF_003254395.2 (Amel_HAv3.1) |
| *Athalia rosae* | `amalgkit_athalia_rosae.yaml` | GCF_000515845.1 |
| *Bombus terrestris* | `amalgkit_bombus_terrestris.yaml` | GCF_910591885.1 |
| *Megachile rotundata* | `amalgkit_megachile_rotundata.yaml` | GCF_000226025.1 |
| *Nasonia vitripennis* | `amalgkit_nasonia_vitripennis.yaml` | GCF_009193385.2 |
| *Polistes canadensis* | `amalgkit_polistes_canadensis.yaml` | GCF_022905055.1 |
| *Polistes fuscatus* | `amalgkit_polistes_fuscatus.yaml` | GCF_000220695.1 |

## 11. Pipeline Progress (March 8, 2026)

| Species | Quantified | Failed | Pending | Total | Status |
|---|---:|---:|---:|---:|---|
| *A. mellifera* | 5,000 | 459 | 82 active | 5,541 | 🏃 Active |
| *H. saltator* | 368 | 0 | 0 | 368 | ✅ Complete |
| *S. invicta* | 325 | 10 | 0 | 335 | ✅ Complete |
| *T. americanus* | 30 | 0 | 300 | 330 | ⏳ Pending |
| *C. floridanus* | 292 | 0 | 0 | 292 | ✅ Complete |
| *A. cephalotes* | 219 | 1 | 0 | 220 | ✅ Complete |
| *O. biroi* | 217 | 0 | 0 | 217 | ✅ Complete |
| *C. obscurior* | 151 | 2 | 5 | 158 | ✅ Nearly |
| *T. longispinosus* | 148 | 0 | 0 | 148 | ✅ Complete |
| *L. humile* | 113 | 0 | 0 | 113 | ✅ Complete |
| *T. nylanderi* | 103 | 12 | 0 | 115 | ✅ Complete |
| *M. pharaonis* | 98 | 0 | 0 | 98 | ✅ Complete |
| *P. barbatus* | 70 | 0 | 62 | 132 | ⏳ Pending |
| *A. echinatior* | 44 | 0 | 0 | 44 | ✅ Complete |
| *N. fulva* | 40 | 0 | 0 | 40 | ✅ Complete |
| *W. auropunctata* | 33 | 0 | 0 | 33 | ✅ Complete |
| *F. exsecta* | 23 | 0 | 0 | 23 | ✅ Complete |
| *O. brunneus* | 19 | 0 | 0 | 19 | ✅ Complete |
| *V. emeryi* | 15 | 0 | 0 | 15 | ✅ Complete |
| *T. curvispinosus* | 14 | 0 | 29 | 43 | ⏳ Pending |
| *D. quadriceps* | 13 | 0 | 0 | 13 | ✅ Complete |
| *A. gracilipes* | 7 | 0 | 0 | 7 | ✅ Complete |
| **TOTAL** | **7,342** | **484** | **478** | **8,304** | **88% done** |

## 12. Output Directory Structure

```
output/amalgkit/<species>/
├── work/
│   ├── metadata/
│   │   └── metadata.tsv          # SRA metadata with tissue normalization
│   ├── index/
│   │   └── <species>_transcripts.idx  # Kallisto transcriptome index
│   ├── quant/
│   │   └── <SRR_ID>/
│   │       └── abundance.tsv     # Per-sample Kallisto output
│   ├── merge/
│   │   └── merged_abundance.tsv  # Combined expression matrix
│   ├── curate/
│   │   └── <species>/tables/     # QC-filtered expression data
│   └── sanity/                   # Validation outputs
├── shared/genome/<Species>/      # Reference genome and index files
└── pipeline_progress.db          # SQLite progress tracking database
```

## 13. Reproducibility

### Configuration Version Control

All species configurations are version-controlled in `config/amalgkit/` as YAML files. Each config specifies:
- NCBI assembly accession and annotation release
- Entrez search query for sample discovery
- Per-step parameters (threads, output directories, redo flags)
- Genome FTP URLs for deterministic reference downloads

### Resumability

The SQLite progress database enables:
- **Crash recovery**: Pipeline resumes from the last known state
- **Incremental processing**: New samples added to SRA are automatically discovered and processed
- **Multi-machine coordination**: Database can be synced between local and cloud environments

### Tissue Normalization Reproducibility

All tissue mappings and patches are declaratively specified in:
- `config/amalgkit/tissue_mapping.yaml` — General tissue label harmonization
- `config/amalgkit/tissue_patches.yaml` — Sample/BioProject-level overrides

---

## Citation

If using this pipeline, please cite:

- **Amalgkit**: Fukushima K. amalgkit: a toolkit for amalgamating RNA-seq data. GitHub: https://github.com/kfuku52/amalgkit
- **Kallisto**: Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology. 2016;34(5):525-527.
- **METAINFORMANT**: [Repository URL and publication reference]

---
*Generated March 7, 2026 for METAINFORMANT RNA-seq V2*
