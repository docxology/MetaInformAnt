# AI Agents in Core Infrastructure Development

This document outlines AI assistance in developing METAINFORMANT's core infrastructure and shared utilities.

## AI Contributions

### Core Architecture Design
**Code Assistant Agent** (grok-code-fast-1) designed:
- Modular core infrastructure organized by functionality
- Consistent API patterns across all utility modules
- Comprehensive error handling and validation frameworks
- Performance optimization patterns for large-scale data processing

### Infrastructure Components
**Code Assistant Agent** implemented:
- Configuration management with YAML/TOML parsing and environment overrides
- Comprehensive I/O utilities (JSON, JSONL, CSV, TSV, Parquet, downloads)
- Structured logging framework with context support and multiple outputs
- Parallel processing utilities with thread management and progress tracking
- Path handling and security validation with containment checks
- JSON-based caching with TTL support and thread safety
- Config-driven processing workflows with validation and error reporting
- Database integration helpers for PostgreSQL connections
- Text processing utilities for normalization and encoding
- Content hashing utilities for file integrity verification
- Symbolic mapping and context discovery utilities for repo-wide navigation
- Symbol indexing and cross-referencing for functions and classes
- Enhanced configuration discovery with schema extraction
- Output pattern discovery and directory structure mapping

### Symbolic Mapping and Discovery
**Code Assistant Agent** implemented:
- `discovery.py`: Comprehensive symbolic mapping and context discovery module
  - Function discovery with AST-based signature extraction
  - Config file discovery with domain filtering and metadata
  - Output pattern identification per module
  - Call graph construction for entry points
  - Symbol usage tracking across repository
  - Module dependency analysis with import extraction
  - Workflow discovery for all domain modules
- `symbols.py`: Symbol indexing and cross-referencing module
  - Function and class indexing across entire repository
  - Symbol definition lookup with fuzzy matching
  - Reference finding with context extraction
  - Signature extraction with type hints
  - Metadata retrieval (docstrings, decorators, parameters)
  - Caching for performance optimization
- Enhanced `config.py`: Configuration discovery extensions
  - Config file discovery with domain filtering
  - Config schema extraction and structure analysis
  - Module-to-config mapping
  - Template listing and discovery
- Enhanced `paths.py`: Output pattern discovery extensions
  - Output pattern discovery per module
  - Output location finding with pattern matching
  - Module output base path resolution
  - Complete output directory structure mapping

### Quality Assurance
**Documentation Agent** assisted with:
- Core utility documentation
- API reference generation
- Usage examples and best practices
- Integration guides and patterns
- Discovery and symbols module documentation
- Context discovery usage patterns and examples

## Development Approach

- **Modular Design**: AI helped design flexible core modules
- **Consistent Patterns**: Established reusable patterns across utilities
- **Error Prevention**: Intelligent validation and type checking
- **Performance Focus**: Efficient algorithms for large-scale data

## Quality Assurance

- Human oversight ensures infrastructure reliability and security
- AI assistance accelerates development while maintaining standards
- Comprehensive testing validates core functionality

This core infrastructure provides a solid foundation for METAINFORMANT's diverse domain modules.

## Core Function Signatures

### Configuration Management (`config.py`)
- `load_config_file(config_path: str | Path) -> dict[str, Any]`
- `load_mapping_from_file(config_path: str | Path) -> dict[str, Any]`
- `apply_env_overrides(config: Mapping[str, Any], *, prefix: str = "AK") -> dict[str, Any]`
- `merge_configs(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]`
- `coerce_config_types(config: dict[str, Any], type_map: dict[str, type]) -> dict[str, Any]`
- `discover_config_files(repo_root: str | Path, domain: str | None = None) -> list[dict[str, Any]]`
- `get_config_schema(config_path: str | Path) -> dict[str, Any]`
- `find_configs_for_module(module_name: str, repo_root: str | Path | None = None) -> list[dict[str, Any]]`
- `list_config_templates(repo_root: str | Path | None = None) -> list[dict[str, Any]]`

### I/O Operations (`io.py`)
- `ensure_directory(path: str | Path) -> Path`
- `open_text_auto(path: str | Path, mode: str = "rt", encoding: str = "utf-8") -> io.TextIOBase`
- `load_json(path: str | Path) -> Any`
- `dump_json(obj: Any, path: str | Path, *, indent: int | None = None, atomic: bool = True) -> None`
- `dump_json_gz(obj: Any, path: str | Path, *, indent: int | None = None) -> None`
- `load_json_gz(path: str | Path) -> Any`
- `read_parquet(path: str | Path, **kwargs) -> Any`
- `write_parquet(df: Any, path: str | Path, **kwargs) -> None`
- `read_jsonl(path: str | Path) -> Iterator[dict[str, Any]]`
- `write_jsonl(rows: Iterable[Mapping[str, Any]], path: str | Path, *, atomic: bool = True) -> None`
- `read_delimited(path: str | Path, *, delimiter: str = ",") -> Iterator[dict[str, str]]`
- `write_delimited(rows: Iterable[Mapping[str, Any]], path: str | Path, *, delimiter: str = ",", atomic: bool = True) -> None`
- `read_csv(path: str | Path, **kwargs) -> Any`
- `write_csv(data: Any, path: str | Path, **kwargs) -> None`
- `read_tsv(path: str | Path) -> list[list[str]]`
- `write_tsv(data, path: str | Path) -> None`
- `download_file(url: str, dest_path: str | Path, *, chunk_size: int = 8192, timeout: int = 30) -> bool`
- `download_json(url: str, *, timeout: int = 30) -> Any`
- `download_text(url: str, *, timeout: int = 30) -> str | None`
- `download_csv(url: str, *, timeout: int = 30, **kwargs) -> Any`
- `batch_download(urls: list[str], dest_dir: str | Path, *, timeout: int = 30) -> dict[str, bool]`

### Path Management (`paths.py`)
- `expand_and_resolve(path: str | Path) -> Path`
- `is_within(path: str | Path, parent: str | Path) -> bool`
- `ensure_directory(path: Path) -> None`
- `prepare_file_path(file_path: Path) -> None`
- `is_safe_path(path: str) -> bool`
- `get_file_extension(filename: str) -> str`
- `change_extension(path: str, new_extension: str) -> Path`
- `find_files_by_extension(directory: str | Path, extension: str) -> list[Path]`
- `get_file_size(path: str | Path) -> int`
- `get_directory_size(path: str | Path) -> int`
- `sanitize_filename(filename: str) -> str`
- `create_temp_file(suffix: str = "", prefix: str = "tmp", directory: str | Path | None = None) -> Path`
- `discover_output_patterns(module_name: str) -> dict[str, Any]`
- `find_output_locations(repo_root: str | Path, pattern: str | None = None) -> list[Path]`
- `get_module_output_base(module_name: str) -> str`
- `list_output_structure(repo_root: str | Path) -> dict[str, Any]`

### Logging Framework (`logging.py`)
- `get_logger(name: str) -> logging.Logger`
- `setup_logger(name: str, log_file: str | None = None, level: str = "INFO") -> logging.Logger`
- `get_logger_with_level(name: str, level: str | int | None = None) -> logging.Logger`
- `configure_logging_from_env(default_level: str = "INFO") -> None`
- `log_with_metadata(logger: logging.Logger, message: str, metadata: dict, *, level: str = "INFO", structured: bool = False) -> None`

### Parallel Processing (`parallel.py`)
- `thread_map(func: Callable, items: Iterable, max_workers: int | None = None, **kwargs: Any) -> list`
- `thread_map_unordered(func: Callable, items: Iterable, max_workers: int | None = None, **kwargs: Any) -> Iterator`
- `parallel_batch(func: Callable, items: Iterable, batch_size: int = 10, max_workers: int | None = None, **kwargs: Any) -> list`
- `cpu_count() -> int`

### Caching System (`cache.py`) - **NOT IMPLEMENTED**
*Planned: JSON-based caching with TTL support*
- `JsonCache(cache_dir: str | Path, ttl_seconds: int = 3600)`
- `JsonCache.get(key: str) -> Any`
- `JsonCache.set(key: str, value: Any) -> None`
- `JsonCache.clear() -> None`
- `JsonCache.cleanup_expired() -> None`

### Progress Tracking (`progress.py`)
- `progress_bar(iterable: Iterator[Any] | None = None, total: int | None = None, desc: str | None = None, **kwargs: Any) -> tqdm`
- `task_context(task_name: str, total_steps: int | None = None) -> Iterator[Any]`
- `log_progress(current: int, total: int | None, message: str = "") -> None`

### Hashing Utilities (`hash.py`)
- `sha256_bytes(data: bytes) -> str`
- `sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str`
- `deterministic_seed(data: str) -> int`
- `sha256_string(s: str) -> str`
- `file_hash_comparison(file1: str | Path, file2: str | Path) -> bool`
- `hash_directory(path: str | Path, pattern: str = "**/*") -> dict[str, str]`
- `verify_file_integrity(file_path: str | Path, expected_hash: str) -> bool`

### Text Processing (`text.py`)
- `normalize_whitespace(s: str) -> str`
- `slugify(s: str) -> str`
- `safe_filename(name: str) -> str`
- `clean_whitespace(text: str) -> str`
- `remove_control_chars(text: str) -> str`
- `standardize_gene_name(gene_name: str) -> str`
- `format_species_name(species_name: str) -> str`
- `clean_sequence_id(sequence_id: str) -> str`
- `extract_numbers(text: str) -> list[float]`
- `truncate_text(text: str, max_length: int, suffix: str = "...") -> str`
- `count_words(text: str) -> int`
- `extract_email_addresses(text: str) -> list[str]`

### Workflow Management (`workflow.py`) - **NOT IMPLEMENTED**
*Planned: Config-driven processing workflows*
- `download_and_process_data(url: str, processor: Callable, output_dir: str | Path) -> Any`
- `validate_config_file(config_path: str | Path) -> tuple[bool, list[str]]`
- `create_sample_config(output_path: str | Path, sample_type: str = "basic") -> None`
- `run_config_based_workflow(config_path: str | Path, **kwargs) -> dict[str, Any]`
- `BaseWorkflowOrchestrator.__init__(config: dict[str, Any], working_dir: str | Path | None = None)`
- `BaseWorkflowOrchestrator.execute_step(step_name: str, **kwargs) -> dict[str, Any]`
- `BaseWorkflowOrchestrator.run_workflow() -> dict[str, Any]`
- `BaseWorkflowOrchestrator.validate_workflow_config() -> tuple[bool, list[str]]`

### Symbolic Mapping and Discovery (`discovery.py`)
- `discover_functions(module_path: str | Path, pattern: str | None = None) -> list[FunctionInfo]`
- `discover_configs(repo_root: str | Path, domain: str | None = None) -> list[ConfigInfo]`
- `discover_output_patterns(module_name: str) -> OutputPattern`
- `build_call_graph(entry_point: str | Path, repo_root: str | Path | None = None) -> dict[str, list[str]]`
- `find_symbol_usage(symbol_name: str, repo_root: str | Path) -> list[SymbolUsage]`
- `get_module_dependencies(module_path: str | Path) -> ModuleDependency`
- `discover_workflows(repo_root: str | Path | None = None) -> list[dict[str, Any]]`

### Symbol Indexing (`symbols.py`)
- `index_functions(repo_root: str | Path, use_cache: bool = True) -> dict[str, list[SymbolDefinition]]`
- `index_classes(repo_root: str | Path, use_cache: bool = True) -> dict[str, list[SymbolDefinition]]`
- `find_symbol(symbol_name: str, symbol_type: str = "function", repo_root: str | Path | None = None) -> list[SymbolDefinition]`
- `get_symbol_signature(symbol_path: str | Path, symbol_name: str) -> str | None`
- `find_symbol_references(symbol_name: str, repo_root: str | Path) -> list[SymbolReference]`
- `get_symbol_metadata(symbol_path: str | Path, symbol_name: str) -> dict[str, Any]`
- `fuzzy_find_symbol(symbol_name: str, symbol_type: str = "function", repo_root: str | Path | None = None, threshold: float = 0.6) -> list[tuple[str, float]]`

### Error Handling (`errors.py`)
- `retry_with_backoff(func: Callable[..., T], *args: Any, max_attempts: int = 3, backoff_factor: float = 2.0, **kwargs: Any) -> T`
- `error_context(context_msg: str, reraise: bool = True) -> Iterator[None]`
- `safe_execute(func: Callable[..., T], *args: Any, default: T | None = None, **kwargs: Any) -> T | None`
- `validate_not_none(value: Any, name: str = "value") -> None`
- `validate_type(value: Any, expected_type: type | tuple[type, ...], name: str = "value") -> None`
- `validate_range(value: float, min_val: float | None = None, max_val: float | None = None, name: str = "value") -> None`
- `validate_path_exists(path: str | Path, name: str = "path") -> Path`
- `validate_path_is_file(path: str | Path, name: str = "path") -> Path`
- `validate_path_is_dir(path: str | Path, name: str = "path") -> Path`
- `validate_path_within(parent: str | Path, path: str | Path, name: str = "path") -> Path`
- `validate_not_empty(value: str | list | dict, name: str = "value") -> None`
- `validate_schema(data: dict[str, Any], schema: dict[str, Any], name: str = "data") -> None`
- `validate_json_schema(data: dict[str, Any], schema_path: str | Path) -> None`
- `validator(func: Callable[[Any], bool]) -> Callable[[Any], None]`

### Database Integration (`db.py`) - **NOT IMPLEMENTED**
*Planned: Optional database client helpers (PostgreSQL)*
- `PostgresConnection(host: str, port: int, database: str, user: str, password: str)`
- `PostgresConnection.connect() -> psycopg2.extensions.connection`
- `PostgresConnection.execute_query(query: str, params: tuple = ()) -> list[dict]`
- `PostgresConnection.bulk_insert(table: str, columns: list[str], data: list[tuple]) -> None`

### Disk Operations (`disk.py`) - **NOT IMPLEMENTED**
*Planned: Disk space and file system utilities*
- `get_disk_usage(path: str | Path) -> dict[str, Any]`
- `get_free_space(path: str | Path) -> int`
- `ensure_disk_space(path: str | Path, required_bytes: int) -> bool`
- `cleanup_temp_files(temp_dir: str | Path, max_age_hours: int = 24) -> int`
