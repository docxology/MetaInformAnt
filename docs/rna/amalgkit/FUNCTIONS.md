# RNA Module Function Index

Quick reference table for all Python functions in the METAINFORMANT RNA module.

## Amalgkit Step Functions

High-level wrappers for amalgkit CLI subcommands.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `metadata` | `metainformant.rna.amalgkit` | Retrieve RNA-seq metadata from NCBI SRA/ENA | [metadata.md](steps/metadata.md) |
| `integrate` | `metainformant.rna.amalgkit` | Integrate FASTQ paths into metadata | [integrate.md](steps/integrate.md) |
| `config` | `metainformant.rna.amalgkit` | Generate configuration files | [config.md](steps/config.md) |
| `select` | `metainformant.rna.amalgkit` | Filter SRA entries by quality | [select.md](steps/select.md) |
| `getfastq` | `metainformant.rna.amalgkit` | Download and convert SRA to FASTQ | [getfastq.md](steps/getfastq.md) |
| `quant` | `metainformant.rna.amalgkit` | Quantify transcript abundances | [quant.md](steps/quant.md) |
| `merge` | `metainformant.rna.amalgkit` | Merge quantification results | [merge.md](steps/merge.md) |
| `cstmm` | `metainformant.rna.amalgkit` | Cross-species TMM normalization | [cstmm.md](steps/cstmm.md) |
| `curate` | `metainformant.rna.amalgkit` | Quality control and batch correction | [curate.md](steps/curate.md) |
| `csca` | `metainformant.rna.amalgkit` | Cross-species correlation analysis | [csca.md](steps/csca.md) |
| `sanity` | `metainformant.rna.amalgkit` | Validate workflow outputs | [sanity.md](steps/sanity.md) |

## Step Runner Functions

Low-level step execution with explicit directory control.

| Function | Module | Description | Documentation |
|----------|--------|-------------|---------------|
| `run_metadata` | `metainformant.rna.steps.metadata` | Execute metadata step | [metadata.md](steps/metadata.md) |
| `run_integrate` | `metainformant.rna.steps.integrate` | Execute integrate step | [integrate.md](steps/integrate.md) |
| `run_config` | `metainformant.rna.steps.config` | Execute config step | [config.md](steps/config.md) |
| `run_select` | `metainformant.rna.steps.select` | Execute select step | [select.md](steps/select.md) |
| `run_getfastq` | `metainformant.rna.steps.getfastq` | Execute getfastq step | [getfastq.md](steps/getfastq.md) |
| `run_quant` | `metainformant.rna.steps.quant` | Execute quant step | [quant.md](steps/quant.md) |
| `run_merge` | `metainformant.rna.steps.merge` | Execute merge step | [merge.md](steps/merge.md) |
| `run_cstmm` | `metainformant.rna.steps.cstmm` | Execute cstmm step | [cstmm.md](steps/cstmm.md) |
| `run_curate` | `metainformant.rna.steps.curate` | Execute curate step | [curate.md](steps/curate.md) |
| `run_csca` | `metainformant.rna.steps.csca` | Execute csca step | [csca.md](steps/csca.md) |
| `run_sanity` | `metainformant.rna.steps.sanity` | Execute sanity step | [sanity.md](steps/sanity.md) |

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
| `quantify_sample` | `metainformant.rna.steps.quant` | Quantify single sample | [API.md](../API.md#quantify_sample) |
| `convert_sra_to_fastq` | `metainformant.rna.steps.getfastq` | Convert SRA to FASTQ | [API.md](../API.md#convert_sra_to_fastq) |
| `delete_sample_fastqs` | `metainformant.rna.steps.getfastq` | Delete sample FASTQs | [API.md](../API.md#delete_sample_fastqs) |
| `run_download_quant_workflow` | `metainformant.rna.steps.process_samples` | Unified download-quantify workflow (sequential/parallel) | [API.md](../API.md#run_download_quant_workflow) |

## Quick Links

- [Complete API Reference](../API.md) - Detailed function documentation with signatures
- [Step Documentation](steps/README.md) - Comprehensive step guides
- [Workflow Guide](../workflow.md) - Workflow planning and execution
- [Configuration Guide](../CONFIGURATION.md) - Configuration management

