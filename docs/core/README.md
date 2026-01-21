# Core Utilities Documentation

This directory contains documentation for METAINFORMANT's core infrastructure and shared utilities.

## Overview

The core module provides essential utilities used across all METAINFORMANT domains, including configuration management, I/O operations, logging, and parallel processing.

## Core Utilities

- [config](config.md): environment, typed env loaders, Postgres config
- [io](io.md): JSON/JSONL, CSV/TSV, gzip-aware open, mkdir
- [logging](logging.md): consistent logger factory
- [text](text.md): normalization, slugify, safe filenames
- [parallel](parallel.md): thread-based map
- [hash](hash.md): SHA-256 hashing helpers
- [paths](paths.md): path expansion, resolution, and containment checks
- [cache](cache.md): simple JSON cache with TTL
- [db](db.md): PostgreSQL client initialization

```mermaid
graph TD
    subgraph "Configuration Management"
        CONFIGconfig.pyYaml/tomlParsingEnvOverridesValidation[config.py_YAML/TOML parsing_Env overrides_Validation]
        SYMBOLSsymbols.pyFunctionDiscoveryAstParsingSignatureExtraction[symbols.py_Function discovery_AST parsing_Signature extraction]
        DISCOVERYdiscovery.pyConfigDiscoveryWorkflowMappingModuleIntegration[discovery.py_Config discovery_Workflow mapping_Module integration]
    end

    subgraph "I/O & Data"
        IOio.pyJson/csv/parquetGzipHandlingAtomicWrites[io.py_JSON/CSV/Parquet_Gzip handling_Atomic writes]
        CACHEcache.pyTtlCachingThreadSafetyJsonStorage[cache.py_TTL caching_Thread safety_JSON storage]
        DBdb.pyPostgresqlConnectionPoolingQueryExecution[db.py_PostgreSQL_Connection pooling_Query execution]
    end

    subgraph "System Utilities"
        PATHSpaths.pyPathResolutionSecurityValidationContainmentChecks[paths.py_Path resolution_Security validation_Containment checks]
        TEXTtext.pyTextProcessingNormalizationEncoding[text.py_Text processing_Normalization_Encoding]
        HASHhash.pySha-256FileIntegrityContentHashing[hash.py_SHA-256_File integrity_Content hashing]
    end

    subgraph "Execution & Control"
        LOGGINGlogging.pyStructuredLoggingContextSupportMultipleOutputs[logging.py_Structured logging_Context support_Multiple outputs]
        PARALLELparallel.pyThreadPoolsMapOperationsResourceManagement[parallel.py_Thread pools_Map operations_Resource management]
        WORKFLOWworkflow.pyConfig-basedWorkflowsStepOrchestrationErrorHandling[workflow.py_Config-based workflows_Step orchestration_Error handling]
    end

    CONFIG --> IO
    CONFIG --> SYMBOLS
    SYMBOLS --> DISCOVERY
    IO --> CACHE
    IO --> DB
    PATHS --> IO
    TEXT --> IO
    HASH --> CACHE
    LOGGING --> WORKFLOW
    PARALLEL --> WORKFLOW

    CONFIG -->|Configuration| DNAdnaModule[DNA Module]
    IO -->|I/O Operations| RNArnaModule[RNA Module]
    LOGGING -->|Logging| GWASgwasModule[GWAS Module]
    PARALLEL -->|Parallel Processing| MATHmathModule[Math Module]
    PATHS -->|Path Management| VIZ[Visualization]
    CACHE -->|Caching| MLmlModule[ML Module]


    class CONFIG,SYMBOLS,DISCOVERY,IO,CACHE,DB,PATHS,TEXT,HASH,LOGGING,PARALLEL,WORKFLOW core
    class DNA,RNA,GWAS,MATH,VIZ,ML domain
```

## Conventions

- **Outputs**: write to `output/` by default; accept an override path parameter.
- **Configs**: read from `config/` using `core.config` helpers; allow env overrides.
- **Data**: read inputs from `data/` unless explicit paths are provided.

## Documentation Files

### Infrastructure Utilities
- **`cache.md`**: JSON-based caching with TTL support
- **`config.md`**: Configuration management and environment variables
- **`db.md`**: Database integration and connection management
- **`hash.md`**: Content hashing and file integrity
- **`io.md`**: Input/output utilities and file format support
- **`logging.md`**: Consistent logging configuration and best practices
- **`parallel.md`**: Parallel processing and thread management
- **`paths.md`**: Path handling and validation utilities
- **`text.md`**: Text processing and normalization

## Related Source Code

- See `src/metainformant/core/` for implementation details
- See `tests/test_core_*.py` for test coverage
- See `src/metainformant/core/README.md` for module-specific documentation

## Usage Examples

Core utilities are used throughout METAINFORMANT:

```python
from metainformant.core import io, parallel, hash

io.ensure_directory("output")
checksums = parallel.thread_map(hash.sha256_file, ["README.md"])  # list[str]
```

Additional examples:

```python
from metainformant.core import config, io, logging, download_and_process_data

# Load configuration with environment overrides
cfg = config.load_config("config.yaml")

# Setup consistent logging
logger = logging.setup_logger("my_module")

# Robust I/O operations
data = io.read_json("data.json")
io.write_json(data, "output.json")

# High-level config-driven processing
results = download_and_process_data({
    "downloads": {"data": {"url": "https://example.com/data.csv"}},
    "processing": {"analyze": {"type": "analysis"}}
})
```

## Integration

Core utilities are integrated across all domains:
- **Configuration** used by all modules for parameter management
- **I/O operations** provide consistent file handling
- **Logging** ensures uniform error reporting and debugging
- **Path utilities** enforce repository structure conventions
- **Caching** improves performance for expensive operations

## Testing

Comprehensive tests ensure infrastructure reliability:
- Configuration parsing and validation
- I/O format compatibility and error handling
- Logging functionality and performance
- Path handling security and correctness
- Cache operation and TTL functionality

## Contributing

When adding new core functionality:
1. Update relevant utility documentation
2. Add comprehensive tests in `tests/test_core_*.py`
3. Ensure backward compatibility
4. Follow established patterns for error handling and validation

This documentation provides complete coverage of METAINFORMANT's core infrastructure utilities.
