# RNA Module Function Index

Quick reference table for all Python functions in the METAINFORMANT RNA module.

## Amalgkit Step Functions

High-level wrappers for amalgkit CLI subcommands.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `metadata` | `metainformant.rna.amalgkit` | Retrieve RNA-seq metadata from NCBI SRA/ENA | [01_metadata.md](steps/01_metadata.md) |
| `config` | `metainformant.rna.amalgkit` | Generate configuration files | [02_config.md](steps/02_config.md) |
| `select` | `metainformant.rna.amalgkit` | Filter SRA entries by quality | [03_select.md](steps/03_select.md) |
| `getfastq` | `metainformant.rna.amalgkit` | Download and convert SRA to FASTQ | [04_getfastq.md](steps/04_getfastq.md) |
| `integrate` | `metainformant.rna.amalgkit` | Integrate FASTQ paths into metadata | [05_integrate.md](steps/05_integrate.md) |
| `quant` | `metainformant.rna.amalgkit` | Quantify transcript abundances | [06_quant.md](steps/06_quant.md) |
| `merge` | `metainformant.rna.amalgkit` | Merge quantification results | [07_merge.md](steps/07_merge.md) |
| `cstmm` | `metainformant.rna.amalgkit` | Cross-species TMM normalization | [08_cstmm.md](steps/08_cstmm.md) |
| `curate` | `metainformant.rna.amalgkit` | Quality control and batch correction | [09_curate.md](steps/09_curate.md) |
| `csca` | `metainformant.rna.amalgkit` | Cross-species correlation analysis | [10_csca.md](steps/10_csca.md) |
| `sanity` | `metainformant.rna.amalgkit` | Validate workflow outputs | [11_sanity.md](steps/11_sanity.md) |

## Step Runner Functions

Amalgkit steps are invoked via CLI wrappers in `metainformant.rna.amalgkit`. Step execution logic lives in `metainformant.rna.engine.workflow_steps`, FASTQ retrieval in `metainformant.rna.engine.sra_extraction`, and pipeline orchestration in `metainformant.rna.engine.pipeline`.

## Workflow Functions

Workflow planning and execution.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `load_workflow_config` | `metainformant.rna.workflow` | Load configuration from YAML | [API.md](../API.md#load_workflow_config) |
| `plan_workflow` | `metainformant.rna.workflow` | Plan workflow steps (dry-run) | [API.md](../API.md#plan_workflow) |
| `plan_workflow_with_params` | `metainformant.rna.workflow` | Plan workflow with overrides | [API.md](../API.md#plan_workflow_with_params) |
| `execute_workflow` | `metainformant.rna.workflow` | Execute complete workflow | [API.md](../API.md#execute_workflow) |
| `AmalgkitWorkflowConfig` | `metainformant.rna.workflow` | Configuration dataclass | [API.md](../API.md#amalgkitworkflowconfig) |

## Genome Preparation Functions

Genome download, transcriptome preparation, and indexing.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `prepare_genome_for_quantification` | `metainformant.rna.genome_prep` | Complete genome setup pipeline | [genome_preparation.md](genome_preparation.md) |
| `prepare_transcriptome_for_kallisto` | `metainformant.rna.genome_prep` | Extract RNA FASTA from genome | [genome_preparation.md](genome_preparation.md) |
| `build_kallisto_index` | `metainformant.rna.genome_prep` | Build kallisto index | [genome_preparation.md](genome_preparation.md) |
| `find_rna_fasta_in_genome_dir` | `metainformant.rna.genome_prep` | Locate RNA FASTA in genome dir | [genome_preparation.md](genome_preparation.md) |
| `download_rna_fasta_from_ftp` | `metainformant.rna.genome_prep` | Download RNA FASTA from FTP | [genome_preparation.md](genome_preparation.md) |
| `download_cds_fasta_from_ftp` | `metainformant.rna.genome_prep` | Download CDS FASTA from FTP | [genome_preparation.md](genome_preparation.md) |
| `extract_transcripts_from_gff` | `metainformant.rna.genome_prep` | Extract transcripts from GFF | [genome_preparation.md](genome_preparation.md) |
| `get_expected_index_path` | `metainformant.rna.genome_prep` | Get expected index path | [genome_preparation.md](genome_preparation.md) |
| `verify_genome_status` | `metainformant.rna.genome_prep` | Check genome/index status | [genome_preparation.md](genome_preparation.md) |
| `orchestrate_genome_setup` | `metainformant.rna.genome_prep` | Run genome setup for all species | [genome_setup_guide.md](genome_setup_guide.md) |

## Orchestration Functions

Multi-species workflow management and monitoring.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `discover_species_configs` | `metainformant.rna.orchestration` | Find all species configs | [API.md](../API.md#discover_species_configs) |
| `run_workflow_for_species` | `metainformant.rna.orchestration` | Execute workflow for one species | [API.md](../API.md#run_workflow_for_species) |
| `check_workflow_status` | `metainformant.rna.orchestration` | Check workflow completion status | [API.md](../API.md#check_workflow_status) |
| `cleanup_unquantified_samples` | `metainformant.rna.orchestration` | Clean up quantified FASTQs | [API.md](../API.md#cleanup_unquantified_samples) |
| `monitor_workflows` | `metainformant.rna.orchestration` | Monitor multiple workflows | [API.md](../API.md#monitor_workflows) |

## Utility Functions

CLI interaction and parameter handling.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `check_cli_available` | `metainformant.rna.amalgkit` | Check if amalgkit is on PATH | [API.md](../API.md#check_cli_available) |
| `ensure_cli_available` | `metainformant.rna.amalgkit` | Ensure amalgkit available (auto-install) | [API.md](../API.md#ensure_cli_available) |
| `build_cli_args` | `metainformant.rna.amalgkit` | Convert params to CLI args | [API.md](../API.md#build_cli_args) |
| `build_amalgkit_command` | `metainformant.rna.amalgkit` | Build complete command | [API.md](../API.md#build_amalgkit_command) |
| `run_amalgkit` | `metainformant.rna.amalgkit` | Execute amalgkit subcommand | [API.md](../API.md#run_amalgkit) |

## Processing Functions

Sample-level processing pipelines.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `quantify_sample` | `metainformant.rna.engine.workflow_steps` | Quantify single sample | [API.md](../API.md#quantify_sample) |
| `convert_sra_to_fastq` | `metainformant.rna.engine.sra_extraction` | Convert SRA to FASTQ | [API.md](../API.md#convert_sra_to_fastq) |
| `delete_sample_fastqs` | `metainformant.rna.engine.sra_extraction` | Delete sample FASTQs | [API.md](../API.md#delete_sample_fastqs) |
| `run_download_quant_workflow` | `metainformant.rna.engine.pipeline` | Unified download-quantify workflow (sequential/parallel) | [API.md](../API.md#run_download_quant_workflow) |

## Monitoring Functions

Workflow progress and sample status tracking.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `count_quantified_samples` | `metainformant.rna.monitoring` | Count quantified and total samples | [API.md](../API.md#count_quantified_samples) |
| `get_sample_status` | `metainformant.rna.monitoring` | Get detailed status for a single sample | [API.md](../API.md#get_sample_status) |
| `analyze_species_status` | `metainformant.rna.monitoring` | Comprehensive analysis of species workflow status | [API.md](../API.md#analyze_species_status) |
| `find_unquantified_samples` | `metainformant.rna.monitoring` | Find all unquantified samples | [API.md](../API.md#find_unquantified_samples) |
| `check_active_downloads` | `metainformant.rna.monitoring` | Check for samples currently being downloaded | [API.md](../API.md#check_active_downloads) |
| `check_workflow_progress` | `metainformant.rna.monitoring` | Get workflow progress summary | [API.md](../API.md#check_workflow_progress) |
| `assess_all_species_progress` | `metainformant.rna.monitoring` | Assess progress for all species | [API.md](../API.md#assess_all_species_progress) |
| `initialize_progress_tracking` | `metainformant.rna.monitoring` | Initialize progress tracking | [API.md](../API.md#initialize_progress_tracking) |

## Environment Functions

Tool availability and environment validation.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `check_amalgkit` | `metainformant.rna.environment` | Check if amalgkit is available | [API.md](../API.md#check_amalgkit) |
| `check_sra_toolkit` | `metainformant.rna.environment` | Check if SRA Toolkit is installed | [API.md](../API.md#check_sra_toolkit) |
| `check_kallisto` | `metainformant.rna.environment` | Check if kallisto is installed | [API.md](../API.md#check_kallisto) |
| `check_metainformant` | `metainformant.rna.environment` | Check if metainformant package is installed | [API.md](../API.md#check_metainformant) |
| `check_virtual_env` | `metainformant.rna.environment` | Check if running inside a virtual environment | [API.md](../API.md#check_virtual_env) |
| `check_rscript` | `metainformant.rna.environment` | Check if Rscript is available | [API.md](../API.md#check_rscript) |
| `check_dependencies` | `metainformant.rna.environment` | Check all required dependencies | [API.md](../API.md#check_dependencies) |
| `validate_environment` | `metainformant.rna.environment` | Comprehensive environment validation | [API.md](../API.md#validate_environment) |

## Cleanup Functions

Partial download cleanup and file naming fixes.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `cleanup_partial_downloads` | `metainformant.rna.cleanup` | Clean up partial downloads | [API.md](../API.md#cleanup_partial_downloads) |
| `fix_abundance_naming` | `metainformant.rna.cleanup` | Create symlink for abundance file naming | [API.md](../API.md#fix_abundance_naming) |
| `fix_abundance_naming_for_species` | `metainformant.rna.cleanup` | Fix abundance naming for all samples | [API.md](../API.md#fix_abundance_naming_for_species) |

## Discovery Functions

Species discovery and configuration generation.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `search_species_with_rnaseq` | `metainformant.rna.discovery` | Search NCBI SRA for species with RNA-seq data | [API.md](../API.md#search_species_with_rnaseq) |
| `get_genome_info` | `metainformant.rna.discovery` | Get genome assembly information | [API.md](../API.md#get_genome_info) |
| `generate_config_yaml` | `metainformant.rna.discovery` | Generate amalgkit YAML configuration | [API.md](../API.md#generate_config_yaml) |

## Quick Links

- [Complete API Reference](../API.md) - Detailed function documentation with signatures
- [Step Documentation](steps/README.md) - Comprehensive step guides
- [Workflow Guide](../workflow.md) - Workflow planning and execution
- [Configuration Guide](../CONFIGURATION.md) - Configuration management

