# Core Examples

This directory contains foundational examples demonstrating METAINFORMANT's core utilities that are used across all biological analysis domains.

## Overview

These examples teach the fundamental patterns for configuration management, file I/O, logging, and path handling that form the foundation of METAINFORMANT applications.

## Examples

### Configuration Management (`example_config.py`)

Learn how to load configuration files with environment variable overrides.

**Demonstrates:**
- Loading YAML/JSON configuration files
- Environment variable overrides using `CORE_` prefix
- Configuration validation and type coercion
- Config merging patterns

```bash
# Run the config example
python examples/core/example_config.py
```

**Output:** `output/examples/core/config_example.json`

### File I/O Patterns (`example_io.py`)

Master file input/output operations for biological data formats.

**Demonstrates:**
- Reading/writing JSON, CSV, and JSONL files
- Gzip compression support
- Atomic file writing to prevent corruption
- Iterator-based processing for large files

```bash
# Run the I/O example
python examples/core/example_io.py
```

**Output:** `output/examples/core/io_example.{json,csv,jsonl}`

### Logging Setup (`example_logging.py`)

Set up structured logging for bioinformatics workflows.

**Demonstrates:**
- Logger configuration with different levels
- Structured logging with metadata
- Log formatting for different outputs
- Best practices for bioinformatics logging

```bash
# Run the logging example
python examples/core/example_logging.py
```

**Output:** `output/examples/core/logging_example.log`

### Path Management (`example_paths.py`)

Handle file paths safely and efficiently.

**Demonstrates:**
- Path expansion and resolution
- Containment validation to prevent directory traversal
- Safe filename generation
- Directory creation and validation

```bash
# Run the paths example
python examples/core/example_paths.py
```

**Output:** `output/examples/core/paths_example.json`

### Workflow Orchestration (`example_workflow.py`)

Learn basic workflow orchestration patterns.

**Demonstrates:**
- Simple config-based workflow execution
- Error handling and validation
- Progress tracking
- Result aggregation

```bash
# Run the workflow example
python examples/core/example_workflow.py
```

**Output:** `output/examples/core/workflow_results.json`

## Learning Progression

1. **Start Here**: `example_config.py` - Learn configuration management
2. **File Operations**: `example_io.py` - Master data file handling
3. **Observability**: `example_logging.py` - Set up proper logging
4. **Safety**: `example_paths.py` - Handle paths securely
5. **Integration**: `example_workflow.py` - Build simple workflows

## Related Documentation

- **Core Module Docs**: [`docs/core/`](../../docs/core/) - Detailed core utilities documentation
- **Configuration Guide**: [`docs/core/config.md`](../../docs/core/config.md) - Advanced config patterns
- **I/O Guide**: [`docs/core/io.md`](../../docs/core/io.md) - File operation patterns
- **Logging Guide**: [`docs/core/logging.md`](../../docs/core/logging.md) - Logging best practices
