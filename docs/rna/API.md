# RNA Module API Reference

Complete function and method reference for the METAINFORMANT RNA analysis module.

## Quick Navigation

- [Amalgkit Step Functions](#amalgkit-step-functions) - High-level step wrappers
- [Step Runner Functions](#step-runner-functions) - Low-level step execution
- [Workflow Functions](#workflow-functions) - Workflow planning and execution
- [Genome Preparation Functions](#genome-preparation-functions) - Genome setup and indexing
- [Orchestration Functions](#orchestration-functions) - Multi-species workflow management
- [Utility Functions](#utility-functions) - CLI helpers and checks
- [Processing Functions](#processing-functions) - Sample processing pipelines
- [Monitoring Functions](#monitoring-functions) - Workflow progress and sample status tracking
- [Environment Functions](#environment-functions) - Tool availability and environment validation
- [Cleanup Functions](#cleanup-functions) - Partial download cleanup and file naming fixes
- [Discovery Functions](#discovery-functions) - Species discovery and configuration generation

---

## Amalgkit Step Functions

High-level wrapper functions for each amalgkit subcommand. These functions provide a clean Python interface to the amalgkit CLI.

### `metadata`

```python
def metadata(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Retrieve RNA-seq metadata from NCBI SRA/ENA databases.

**Parameters**:
- `params`: Dictionary of amalgkit parameters (e.g., `{"out_dir": "work", "search_string": "..."}`)
- `**kwargs`: Additional arguments passed to `run_amalgkit` (e.g., `work_dir`, `log_dir`, `check`)

**Returns**: `subprocess.CompletedProcess[str]` with execution results

**See Also**: [Step Documentation: metadata](amalgkit/steps/metadata.md)

---

### `integrate`

```python
def integrate(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Integrate FASTQ file paths into metadata tables.

**See Also**: [Step Documentation: integrate](amalgkit/steps/integrate.md)

---

### `config`

```python
def config(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Generate configuration files for metadata selection.

**See Also**: [Step Documentation: config](amalgkit/steps/config.md)

---

### `select`

```python
def select(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Select and filter SRA entries based on quality criteria.

**See Also**: [Step Documentation: select](amalgkit/steps/select.md)

---

### `getfastq`

```python
def getfastq(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Download and convert SRA files to FASTQ format.

**See Also**: [Step Documentation: getfastq](amalgkit/steps/getfastq.md)

---

### `quant`

```python
def quant(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Quantify transcript abundances using kallisto or salmon.

**See Also**: [Step Documentation: quant](amalgkit/steps/quant.md)

---

### `merge`

```python
def merge(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Merge per-sample quantification results into expression matrices.

**See Also**: [Step Documentation: merge](amalgkit/steps/merge.md)

---

### `cstmm`

```python
def cstmm(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Cross-species TMM (Trimmed Mean of M-values) normalization.

**See Also**: [Step Documentation: cstmm](amalgkit/steps/cstmm.md)

---

### `curate`

```python
def curate(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Quality control, outlier detection, and batch effect correction.

**See Also**: [Step Documentation: curate](amalgkit/steps/curate.md)

---

### `csca`

```python
def csca(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Cross-species correlation analysis and visualization.

**See Also**: [Step Documentation: csca](amalgkit/steps/csca.md)

---

### `sanity`

```python
def sanity(
    params: AmalgkitParams | None = None,
    **kwargs: Any
) -> subprocess.CompletedProcess[str]
```

**Purpose**: Validate workflow outputs and check data integrity.

**See Also**: [Step Documentation: sanity](amalgkit/steps/sanity.md)

---

## Step Runner Functions

Low-level step execution functions that provide more control over step execution. These are used internally by workflow orchestration but can also be called directly.

### `run_metadata`

```python
def run_metadata(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.metadata`

**Purpose**: Execute metadata retrieval step with explicit directory control.

---

### `run_integrate`

```python
def run_integrate(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.integrate`

---

### `run_config`

```python
def run_config(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.config`

---

### `run_select`

```python
def run_select(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.select`

---

### `run_getfastq`

```python
def run_getfastq(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.getfastq`

**Note**: Includes robust retry logic and fallback mechanisms for failed downloads.

---

### `run_quant`

```python
def run_quant(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.quant`

---

### `run_merge`

```python
def run_merge(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.merge`

---

### `run_cstmm`

```python
def run_cstmm(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.cstmm`

---

### `run_curate`

```python
def run_curate(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.curate`

---

### `run_csca`

```python
def run_csca(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.csca`

---

### `run_sanity`

```python
def run_sanity(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.steps.sanity`

---

## Workflow Functions

Functions for planning and executing complete RNA-seq workflows.

### `load_workflow_config`

```python
def load_workflow_config(config_file: str | Path) -> AmalgkitWorkflowConfig
```

**Module**: `metainformant.rna.workflow`

**Purpose**: Load workflow configuration from YAML file.

**Parameters**:
- `config_file`: Path to YAML configuration file

**Returns**: `AmalgkitWorkflowConfig` dataclass instance

**See Also**: [Configuration Guide](CONFIGURATION.md)

---

### `plan_workflow`

```python
def plan_workflow(config: AmalgkitWorkflowConfig) -> list[tuple[str, AmalgkitParams]]
```

**Module**: `metainformant.rna.workflow`

**Purpose**: Generate ordered list of workflow steps with parameters (dry-run planning).

**Parameters**:
- `config`: Workflow configuration

**Returns**: List of `(step_name, params)` tuples in execution order

**See Also**: [Workflow Guide](workflow.md)

---

### `plan_workflow_with_params`

```python
def plan_workflow_with_params(
    config: AmalgkitWorkflowConfig,
    step_params: dict[str, AmalgkitParams],
) -> list[tuple[str, AmalgkitParams]]
```

**Module**: `metainformant.rna.workflow`

**Purpose**: Plan workflow with explicit per-step parameter overrides.

---

### `execute_workflow`

```python
def execute_workflow(
    config: AmalgkitWorkflowConfig,
    *,
    check: bool = False
) -> list[int]
```

**Module**: `metainformant.rna.workflow`

**Purpose**: Execute complete workflow from configuration.

**Parameters**:
- `config`: Workflow configuration
- `check`: If True, raise exception on step failure

**Returns**: List of return codes (one per step)

**See Also**: [Workflow Guide](workflow.md)

---

### `AmalgkitWorkflowConfig`

```python
@dataclass
class AmalgkitWorkflowConfig:
    work_dir: Path
    threads: int = 6
    species_list: list[str] = field(default_factory=list)
    log_dir: Path | None = None
    manifest_path: Path | None = None
    per_step: dict[str, AmalgkitParams] = field(default_factory=dict)
    auto_install_amalgkit: bool = False
    genome: dict[str, Any] | None = None
    filters: dict[str, Any] = field(default_factory=dict)
```

**Module**: `metainformant.rna.workflow`

**Purpose**: Configuration dataclass for workflow execution.

**See Also**: [Configuration Guide](CONFIGURATION.md)

---

## Genome Preparation Functions

Functions for downloading genomes, preparing transcriptomes, and building kallisto indexes.

### `prepare_genome_for_quantification`

```python
def prepare_genome_for_quantification(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
    build_index: bool = True,
    kmer_size: int = 31,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Complete genome preparation pipeline (download → extract → index).

**Returns**: Dictionary with `success`, `fasta_path`, `index_path`, `error` keys

**See Also**: [Genome Preparation](amalgkit/genome_preparation.md)

---

### `prepare_transcriptome_for_kallisto`

```python
def prepare_transcriptome_for_kallisto(
    genome_dir: Path,
    species_name: str,
    work_dir: Path,
    *,
    accession: str | None = None,
) -> Path | None
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Extract and prepare RNA FASTA file from genome package.

**Returns**: Path to prepared FASTA file or None if failed

---

### `build_kallisto_index`

```python
def build_kallisto_index(
    fasta_path: Path,
    index_path: Path,
    *,
    kmer_size: int = 31,
    check_existing: bool = True,
) -> bool
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Build kallisto index from transcriptome FASTA.

**Parameters**:
- `fasta_path`: Path to transcriptome FASTA file
- `index_path`: Output path for kallisto index
- `kmer_size`: K-mer size (31 for standard reads, 23 for short reads)
- `check_existing`: Skip if index already exists

**Returns**: True if index was built successfully or already exists

---

### `find_rna_fasta_in_genome_dir`

```python
def find_rna_fasta_in_genome_dir(
    genome_dir: Path,
    accession: str
) -> Path | None
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Locate RNA FASTA file in extracted genome directory.

---

### `download_rna_fasta_from_ftp`

```python
def download_rna_fasta_from_ftp(
    ftp_url: str,
    genome_dir: Path,
    accession: str,
    assembly_name: str | None = None,
    config: dict[str, Any] | None = None,
) -> Path | None
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Download RNA FASTA directly from NCBI FTP.

---

### `download_cds_fasta_from_ftp`

```python
def download_cds_fasta_from_ftp(
    ftp_url: str,
    genome_dir: Path,
    accession: str,
    assembly_name: str | None = None,
) -> Path | None
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Download CDS FASTA directly from NCBI FTP.

---

### `extract_transcripts_from_gff`

```python
def extract_transcripts_from_gff(
    gff_path: Path,
    genome_fasta: Path,
    output_fasta: Path,
) -> bool
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Extract transcript sequences from GFF annotation using gffread.

---

### `get_expected_index_path`

```python
def get_expected_index_path(work_dir: Path, species_name: str) -> Path
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Get expected kallisto index path for a species.

---

### `verify_genome_status`

```python
def verify_genome_status(
    genome_dir: Path,
    work_dir: Path,
    species_name: str,
    accession: str | None = None,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Check genome download, transcriptome preparation, and index status.

**Returns**: Dictionary with status flags and paths

---

### `orchestrate_genome_setup`

```python
def orchestrate_genome_setup(
    config_dir: Path = Path("config/amalgkit"),
    *,
    species: str | None = None,
    skip_download: bool = False,
    skip_prepare: bool = False,
    skip_build: bool = False,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.genome_prep`

**Purpose**: Run complete genome setup pipeline for all or specific species.

**See Also**: [Genome Setup Guide](amalgkit/genome_setup_guide.md)

---

## Orchestration Functions

Functions for managing multi-species workflows and monitoring progress.

### `discover_species_configs`

```python
def discover_species_configs(
    config_dir: Path = Path("config/amalgkit")
) -> dict[str, dict[str, Any]]
```

**Module**: `metainformant.rna.orchestration`

**Purpose**: Discover all species configuration files in config directory.

**Returns**: Dictionary mapping species names to config dictionaries

---

### `run_workflow_for_species`

```python
def run_workflow_for_species(
    config_path: Path,
    steps: Sequence[str] | None = None,
    *,
    check: bool = False,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.orchestration`

**Purpose**: Execute workflow steps for a single species.

**Parameters**:
- `config_path`: Path to species workflow config file
- `steps`: List of steps to run (default: all steps)
- `check`: If True, stop on first failure

**Returns**: Dictionary with `success`, `completed`, `failed`, `return_codes` keys

---

### `check_workflow_status`

```python
def check_workflow_status(
    config_path: Path,
    *,
    detailed: bool = False,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.orchestration`

**Purpose**: Check workflow status for a species. Unified interface that delegates to monitoring functions.

**Parameters**:
- `config_path`: Path to species workflow config file
- `detailed`: If True, return detailed analysis via `analyze_species_status()`; if False, return progress summary via `check_workflow_progress()`

**Returns**: Dictionary with status information (format depends on `detailed` parameter)

**Note**: This is a convenience wrapper. For direct access:
- Use `check_workflow_progress()` for progress summary
- Use `analyze_species_status()` for detailed analysis

---

### `cleanup_unquantified_samples`

```python
def cleanup_unquantified_samples(
    config_path: Path,
    *,
    log_dir: Path | None = None,
) -> tuple[int, int]
```

**Module**: `metainformant.rna.orchestration`

**Purpose**: Quantify downloaded samples and cleanup FASTQs. Finds all samples with FASTQ files but no quantification results, quantifies each sample, and deletes FASTQ files after successful quantification.

**Parameters**:
- `config_path`: Path to species workflow config file
- `log_dir`: Optional log directory

**Returns**: Tuple of `(quantified_count, failed_count)`

---

### `monitor_workflows`

```python
def monitor_workflows(
    species_configs: dict[str, Path],
    watch_interval: int = 60,
) -> None
```

**Module**: `metainformant.rna.orchestration`

**Purpose**: Real-time monitoring of multiple species workflows.

**Parameters**:
- `species_configs`: Dictionary mapping species_id -> config_path
- `watch_interval`: Update interval in seconds (default: 60)

**Note**: This function runs continuously until interrupted (Ctrl+C). It displays a real-time dashboard showing progress for all monitored species.

---

## Utility Functions

Core utilities for CLI interaction and parameter handling.

### `check_cli_available`

```python
def check_cli_available() -> tuple[bool, str]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Check if `amalgkit` CLI is available on PATH.

**Returns**: Tuple of `(is_available: bool, help_text_or_error: str)`

---

### `ensure_cli_available`

```python
def ensure_cli_available(
    *,
    auto_install: bool = False
) -> tuple[bool, str, dict | None]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Ensure amalgkit CLI is available, optionally attempting auto-install.

**Returns**: Tuple of `(ok: bool, message: str, install_record: dict | None)`

---

### `build_cli_args`

```python
def build_cli_args(
    params: AmalgkitParams | None,
    *,
    for_cli: bool = False
) -> list[str]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Convert parameter dictionary to CLI argument list.

**Parameters**:
- `params`: Parameter mapping
- `for_cli`: If True, use snake_case flags (for actual CLI); if False, use kebab-case (for display)

**Returns**: List of CLI argument strings

---

### `build_amalgkit_command`

```python
def build_amalgkit_command(
    subcommand: str,
    params: AmalgkitParams | None = None
) -> list[str]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Build complete amalgkit command as token list.

**Example**: `build_amalgkit_command("metadata", {"threads": 8})` → `["amalgkit", "metadata", "--threads", "8"]`

---

### `run_amalgkit`

```python
def run_amalgkit(
    subcommand: str,
    params: AmalgkitParams | None = None,
    *,
    work_dir: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = False,
    capture_output: bool = True,
    log_dir: str | Path | None = None,
    step_name: str | None = None,
) -> subprocess.CompletedProcess[str]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Execute amalgkit subcommand with optional logging and directory management.

**Parameters**:
- `subcommand`: Amalgkit subcommand name (e.g., "metadata", "quant")
- `params`: Parameter dictionary
- `work_dir`: Working directory (created if missing)
- `env`: Additional environment variables
- `check`: Raise exception on non-zero exit
- `capture_output`: Capture stdout/stderr
- `log_dir`: Directory for timestamped log files
- `step_name`: Optional label for log filenames

**Returns**: `subprocess.CompletedProcess[str]` with execution results

---

## Processing Functions

Functions for sample-level processing pipelines.

### `quantify_sample`

```python
def quantify_sample(
    sample_id: str,
    metadata_rows: list[dict[str, Any]],
    quant_params: Mapping[str, Any],
    *,
    log_dir: Path | None = None,
    step_name: str | None = None,
) -> tuple[bool, str, Path | None]
```

**Module**: `metainformant.rna.steps.quant`

**Purpose**: Quantify a single sample using amalgkit quant.

**Returns**: Tuple of `(success: bool, message: str, abundance_file: Path | None)`

---

### `convert_sra_to_fastq`

```python
def convert_sra_to_fastq(
    sample_id: str,
    sra_file: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    log_dir: Path | None = None,
) -> tuple[bool, str, list[Path]]
```

**Module**: `metainformant.rna.steps.getfastq`

**Purpose**: Convert a local SRA file to FASTQ format. Prefers `parallel-fastq-dump` (works better with local files) and falls back to `fasterq-dump` if needed. Automatically compresses output FASTQ files.

**Args**:
- `sample_id`: SRA accession ID (e.g., "SRR1234567")
- `sra_file`: Path to the SRA file
- `output_dir`: Directory where FASTQ files should be written
- `threads`: Number of threads for conversion (default: 4)
- `log_dir`: Optional directory for log files

**Returns**: Tuple of `(success: bool, message: str, fastq_files: list[Path])`. `fastq_files` contains paths to created FASTQ files (may be empty if failed).

**Notes**:
- Automatically detects and uses the real `fasterq-dump` binary (not wrapper scripts)
- Passes `--size-check off` to `fasterq-dump` to prevent "disk-limit exceeded" errors
- Automatically compresses output using `pigz` or `gzip`
- Checks for existing FASTQ files before conversion

---

### `delete_sample_fastqs`

```python
def delete_sample_fastqs(
    sample_id: str,
    fastq_dir: Path
) -> None
```

**Module**: `metainformant.rna.steps.getfastq`

**Purpose**: Delete FASTQ files for a specific sample.

---

### `run_download_quant_workflow`

```python
def run_download_quant_workflow(
    metadata_path: str | Path,
    getfastq_params: Mapping[str, Any] | None = None,
    quant_params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    num_workers: int = 1,
    max_samples: int | None = None,
    skip_completed: bool = True,
    progress_monitor: DownloadProgressMonitor | None = None,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.steps.process_samples`

**Purpose**: Unified function for download-quantify-delete workflows. Supports both sequential and parallel processing modes.

**Parameters**:
- `metadata_path`: Path to metadata TSV with sample list
- `getfastq_params`: Parameters for amalgkit getfastq step
- `quant_params`: Parameters for amalgkit quant step
- `work_dir`: Working directory for amalgkit commands
- `log_dir`: Directory for step logs
- `num_workers`: Number of parallel download workers (default: 1 for sequential mode)
  - `num_workers=1`: Sequential mode (one sample at a time, maximum disk efficiency)
  - `num_workers>1`: Parallel mode (N downloads in parallel, sequential quantification)
- `max_samples`: Optional limit on number of samples to process
- `skip_completed`: If True, skip samples that are already quantified (sequential mode only)
- `progress_monitor`: Optional progress monitor for tracking downloads

**Returns**: Dictionary with processing statistics:
- `total_samples`: Total number of samples
- `processed`: Number of samples successfully processed
- `skipped`: Number of samples skipped (already done)
- `failed`: Number of samples that failed
- `failed_runs`: List of run IDs that failed

**Processing Modes**:
- **Sequential** (`num_workers=1`): Process one sample at a time. Maximum disk efficiency: only one sample's FASTQs exist at a time.
- **Parallel** (`num_workers>1`): Multiple download workers fetch FASTQ files in parallel, then a single quantification worker processes them sequentially. Maximizes throughput while preventing disk exhaustion.

**See Also**: [Architecture Documentation](ARCHITECTURE.md#processing-workflows)

---

## Type Definitions

### `AmalgkitParams`

```python
AmalgkitParams = Mapping[str, Any]
```

**Module**: `metainformant.rna.amalgkit`

**Purpose**: Type alias for amalgkit parameter dictionaries.

---

## Monitoring Functions

Functions for tracking workflow progress and sample status.

### `count_quantified_samples`

```python
def count_quantified_samples(config_path: Path) -> tuple[int, int]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Count quantified and total samples for a species.

**Returns**: Tuple of `(quantified_count, total_count)`

---

### `get_sample_status`

```python
def get_sample_status(config_path: Path, sample_id: str) -> dict[str, Any]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Get detailed status for a single sample.

**Returns**: Dictionary with status information:
- `quantified`: bool
- `has_fastq`: bool
- `has_sra`: bool
- `is_downloading`: bool
- `status`: str ("quantified", "downloading", "has_fastq", "has_sra", "undownloaded")

---

### `analyze_species_status`

```python
def analyze_species_status(config_path: Path) -> dict[str, Any]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Comprehensive analysis of species workflow status.

**Returns**: Dictionary with comprehensive status information:
- `total_in_metadata`: int
- `quantified`: int
- `quantified_and_deleted`: int
- `quantified_not_deleted`: int
- `downloading`: int
- `failed_download`: int
- `undownloaded`: int
- `categories`: dict mapping category -> list of sample_ids

---

### `find_unquantified_samples`

```python
def find_unquantified_samples(config_path: Path) -> list[str]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Find all unquantified samples.

**Returns**: List of sample IDs that are not quantified

---

### `check_active_downloads`

```python
def check_active_downloads() -> set[str]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Check for samples currently being downloaded.

**Returns**: Set of sample IDs that are actively downloading

---

### `check_workflow_progress`

```python
def check_workflow_progress(config_path: Path) -> dict[str, Any]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Get workflow progress summary for a species.

**Returns**: Dictionary with progress information:
- `quantified`: int
- `total`: int
- `percentage`: float
- `remaining`: int
- `downloading`: int (number of samples currently downloading)
- `has_files`: int (number of samples with downloaded files but not quantified)

**Note**: This is called by `check_workflow_status()` when `detailed=False`. Use `check_workflow_status()` for the unified interface.

---

### `assess_all_species_progress`

```python
def assess_all_species_progress(
    config_dir: Path,
    *,
    repo_root: Path | None = None,
) -> dict[str, dict[str, Any]]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Assess progress for all species in config directory.

**Returns**: Dictionary mapping species_id -> progress information

---

### `initialize_progress_tracking`

```python
def initialize_progress_tracking(
    config_path: Path,
    *,
    tracker=None,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.monitoring`

**Purpose**: Initialize progress tracking for a species.

**Returns**: Dictionary with initialization results

---

## Environment Functions

Functions for checking tool availability and environment validation.

### `check_amalgkit`

```python
def check_amalgkit() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if amalgkit is available and get version.

**Returns**: Tuple of `(is_available: bool, message: str)`

---

### `check_sra_toolkit`

```python
def check_sra_toolkit() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if SRA Toolkit is installed.

**Returns**: Tuple of `(is_available: bool, message: str)`

---

### `check_kallisto`

```python
def check_kallisto() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if kallisto is installed.

**Returns**: Tuple of `(is_available: bool, message: str)`

---

### `check_metainformant`

```python
def check_metainformant() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if metainformant package is installed.

**Returns**: Tuple of `(is_available: bool, message: str)`

---

### `check_virtual_env`

```python
def check_virtual_env() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if running inside a virtual environment.

**Returns**: Tuple of `(is_in_venv: bool, message: str)`

---

### `check_rscript`

```python
def check_rscript() -> tuple[bool, str]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check if Rscript is available.

**Returns**: Tuple of `(is_available: bool, message: str)`

---

### `check_dependencies`

```python
def check_dependencies() -> dict[str, tuple[bool, str]]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Check all required dependencies for RNA-seq workflows.

**Returns**: Dictionary mapping dependency name -> `(is_available: bool, message: str)`

---

### `validate_environment`

```python
def validate_environment() -> dict[str, Any]
```

**Module**: `metainformant.rna.environment`

**Purpose**: Comprehensive environment validation.

**Returns**: Dictionary with validation results:
- `all_passed`: bool
- `dependencies`: dict mapping name -> `(is_available, message)`
- `recommendations`: list of strings with recommendations

---

## Cleanup Functions

Functions for cleaning up partial downloads and fixing file naming issues.

### `cleanup_partial_downloads`

```python
def cleanup_partial_downloads(
    config_path: Path,
    *,
    dry_run: bool = False,
) -> dict[str, Any]
```

**Module**: `metainformant.rna.cleanup`

**Purpose**: Clean up partial downloads for a species.

**Parameters**:
- `config_path`: Path to species workflow config file
- `dry_run`: If True, only report what would be deleted

**Returns**: Dictionary with `deleted`, `freed_mb`, `errors` keys

---

### `fix_abundance_naming`

```python
def fix_abundance_naming(quant_dir: Path, sample_id: str) -> bool
```

**Module**: `metainformant.rna.cleanup`

**Purpose**: Create symlink from `abundance.tsv` to `{SRR}_abundance.tsv` for amalgkit merge.

**Parameters**:
- `quant_dir`: Directory containing quantification results
- `sample_id`: Sample ID (e.g., 'SRR1234567')

**Returns**: True if symlink was created or already exists

---

### `fix_abundance_naming_for_species`

```python
def fix_abundance_naming_for_species(
    config_path: Path,
) -> tuple[int, int]
```

**Module**: `metainformant.rna.cleanup`

**Purpose**: Fix abundance naming for all samples in a species.

**Returns**: Tuple of `(created_count, already_exists_count)`

---

## Discovery Functions

Functions for discovering species with RNA-seq data and generating configurations.

### `search_species_with_rnaseq`

```python
def search_species_with_rnaseq(
    search_query: str,
    *,
    max_records: int = 10000,
) -> dict[str, dict[str, Any]]
```

**Module**: `metainformant.rna.discovery`

**Purpose**: Search NCBI SRA for species with RNA-seq data.

**Parameters**:
- `search_query`: NCBI Entrez search query
- `max_records`: Maximum number of records to retrieve

**Returns**: Dictionary mapping species names to metadata

**Raises**: `ImportError` if Biopython is not available

---

### `get_genome_info`

```python
def get_genome_info(taxonomy_id: str, species_name: str) -> dict[str, Any] | None
```

**Module**: `metainformant.rna.discovery`

**Purpose**: Get genome assembly information for a species.

**Parameters**:
- `taxonomy_id`: NCBI taxonomy ID
- `species_name`: Scientific name

**Returns**: Genome information dictionary or None

---

### `generate_config_yaml`

```python
def generate_config_yaml(
    species_name: str,
    species_data: dict[str, Any],
    genome_info: dict[str, Any] | None = None,
    *,
    repo_root: Path | None = None,
) -> str
```

**Module**: `metainformant.rna.discovery`

**Purpose**: Generate amalgkit YAML configuration for a species.

**Parameters**:
- `species_name`: Scientific name
- `species_data`: RNA-seq metadata
- `genome_info`: Genome assembly metadata (optional)
- `repo_root`: Repository root directory for paths (optional)

**Returns**: YAML configuration string

---

## See Also

- [Function Index](amalgkit/FUNCTIONS.md) - Quick reference table
- [Step Documentation](amalgkit/steps/README.md) - Detailed step guides
- [Workflow Guide](workflow.md) - Workflow planning and execution
- [Configuration Guide](CONFIGURATION.md) - Configuration management

