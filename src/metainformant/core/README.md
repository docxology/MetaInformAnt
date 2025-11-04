# Core Utilities Module

The `core` module provides essential shared utilities and infrastructure used across all METAINFORMANT domains. These utilities handle common operations like configuration management, I/O operations, logging, and parallel processing.

## Overview

This module contains foundational components that are used throughout the METAINFORMANT package. Each submodule provides specific functionality while maintaining consistency and reusability across the entire codebase.

## Submodules

### Configuration (`config.py`)
Centralized configuration management with environment variable overrides.

**Key Features:**
- YAML and TOML configuration file support (optional dependencies)
- Environment variable integration and validation
- Database connection configuration (PostgreSQL)
- Configuration inheritance and merging
- Runtime configuration validation

**Usage:**
```python
from metainformant.core import config

# Load configuration from file with environment overrides
cfg = config.load_config("config.yaml")

# Access configuration values
db_config = cfg.database
threads = cfg.compute.threads

# Validate configuration
config.validate_config(cfg)
```

### I/O Operations (`io.py`)
Robust file I/O utilities with support for multiple formats and compression.

**Key Features:**
- JSON/JSONL reading and writing
- CSV/TSV processing with type inference
- Gzip compression support
- Atomic write operations
- Directory creation and path validation

**Usage:**
```python
from metainformant.core import io

# JSON operations
data = io.read_json("data.json")
io.write_json(data, "output.json")

# CSV operations with type inference
df = io.read_csv("data.csv")
io.write_csv(df, "output.csv")

# Compressed file handling
io.write_json_gz(data, "compressed.json.gz")
```

### Logging (`logging.py`)
Consistent, structured logging across all modules.

**Key Features:**
- Structured logging with consistent formatting
- Multiple output targets (console, file, etc.)
- Log level management and filtering
- Context-aware logging for distributed operations

**Usage:**
```python
from metainformant.core import logging

# Setup logger
logger = logging.setup_logger("my_module")

# Log messages with context
logger.info("Processing started", extra={"items": 100})
logger.error("Processing failed", extra={"error": str(exc)})
```

### Parallel Processing (`parallel.py`)
Thread-based parallel execution with order preservation.

**Key Features:**
- Thread pool execution with result ordering
- Progress tracking and cancellation
- Memory-efficient processing of large datasets
- Exception handling and aggregation

**Usage:**
```python
from metainformant.core import parallel

# Parallel map with progress tracking
def process_item(item):
    return item * 2

results = parallel.thread_map(process_item, items, max_workers=4)
```

### Path Management (`paths.py`)
Path expansion, resolution, and containment validation.

**Key Features:**
- Path expansion with environment variables and user home
- Path resolution and normalization
- Containment checks for security
- Cross-platform path handling

**Usage:**
```python
from metainformant.core import paths

# Expand paths with environment variables
expanded = paths.expand_path("~/data/{species}")
resolved = paths.resolve_path("./relative/path")

# Security containment checks
if paths.is_contained(resolved, "/allowed/directory"):
    # Safe to use path
    pass
```

### Caching (`cache.py`)
JSON-based caching with TTL (Time To Live) support.

**Key Features:**
- Disk-based JSON caching
- TTL-based expiration
- Cache size management
- Thread-safe operations

**Usage:**
```python
from metainformant.core import cache

# Create cache with TTL
cache_obj = cache.JSONCache("cache_dir", ttl_seconds=3600)

# Cache operations
cache_obj.set("key", {"data": "value"})
value = cache_obj.get("key")
```

### Text Processing (`text.py`)
Text normalization and processing utilities.

**Key Features:**
- Unicode normalization and cleaning
- Case conversion with locale awareness
- Whitespace normalization
- Text encoding detection and conversion

**Usage:**
```python
from metainformant.core import text

# Text normalization
clean_text = text.normalize_text("  Mixed\tCase\nText  ")
title_case = text.title_case("species name")
```

### Hash Functions (`hash.py`)
Content and file hashing utilities.

**Key Features:**
- Multiple hash algorithms (MD5, SHA256, etc.)
- File content hashing
- String hashing with encoding support
- Hash comparison utilities

**Usage:**
```python
from metainformant.core import hash

# File hashing
file_hash = hash.hash_file("data.txt")
content_hash = hash.hash_string("content to hash")
```

### Database Integration (`db.py`) - Optional
Database client helpers for PostgreSQL integration.

**Key Features:**
- Connection pooling and management
- Query execution with parameter binding
- Transaction management
- Result set processing

**Usage:**
```python
from metainformant.core import db

# Database operations (requires database configuration)
with db.get_connection(db_config) as conn:
    results = db.execute_query(conn, "SELECT * FROM table WHERE id = %s", (123,))
```

## Architecture Principles

### Defensive Imports
Optional dependencies are imported defensively to avoid hard failures during unrelated operations. This ensures the package remains functional even when optional dependencies are missing.

### Consistency
All utilities follow consistent patterns for configuration, error handling, and return values.

### Performance
Operations are designed to be memory-efficient and suitable for large-scale biological data processing.

### Extensibility
Utilities are designed to be easily extended with new formats, algorithms, or integrations.

## Error Handling

All core utilities implement comprehensive error handling:
- Clear error messages with context
- Appropriate exception types
- Graceful degradation for optional features
- Input validation and sanitization

## Testing

Each submodule includes comprehensive tests covering:
- Normal operation scenarios
- Error conditions and edge cases
- Performance characteristics
- Integration with other modules

## Dependencies

- **Required**: Standard library only (pathlib, json, csv, etc.)
- **Optional**: PyYAML (for YAML config), tomllib (Python 3.11+ for TOML)
- **Database**: psycopg2 (optional, for PostgreSQL integration)

## Usage in Other Modules

Core utilities are designed to be imported and used throughout the METAINFORMANT package. Here are comprehensive examples showing how core utilities are used across different modules:

### Configuration Management Across Workflows

```python
from metainformant.core import config
from metainformant.rna import workflow
from metainformant.gwas import workflow as gwas_workflow

# Load configuration for RNA-seq workflow
rna_cfg = config.load_config("config/rna_config.yaml")
rna_results = workflow.execute_workflow(rna_cfg)

# Load configuration for GWAS workflow
gwas_cfg = config.load_config("config/gwas_config.yaml")
gwas_results = gwas_workflow.run_gwas_workflow(gwas_cfg)

# Environment variable overrides
import os
os.environ["METAINFORMANT_COMPUTE_THREADS"] = "16"
os.environ["METAINFORMANT_OUTPUT_DIR"] = "output/custom"
# Configuration automatically uses these overrides
cfg = config.load_config("config.yaml")
```

### I/O Operations in Different Domain Modules

```python
from metainformant.core import io, paths
from metainformant.dna import sequences
from metainformant.rna import workflow
from metainformant.multiomics import MultiOmicsData

# DNA module: Reading sequences with core I/O
seq_path = paths.expand_path("data/sequences.fasta")
seqs = sequences.read_fasta(seq_path)

# RNA module: Writing expression data
expression_data = workflow.extract_expression("expression.tsv")
output_path = paths.join_paths("output", "rna", "expression.json")
io.write_json(expression_data.to_dict(), output_path)

# Multiomics module: Reading multiple omics files
genomics_path = paths.expand_path("data/genomics.csv")
transcriptomics_path = paths.expand_path("data/transcriptomics.tsv")
genomics_df = io.read_csv(genomics_path)
transcriptomics_df = io.read_csv(transcriptomics_path, delimiter="\t")

# Create multi-omics dataset
omics_data = MultiOmicsData(genomics=genomics_df, transcriptomics=transcriptomics_df)

# Compressed I/O for large datasets
large_data = {"results": [...]}
io.write_json_gz(large_data, "output/large_results.json.gz")
loaded_data = io.read_json_gz("output/large_results.json.gz")
```

### Logging Patterns in Multi-Module Pipelines

```python
from metainformant.core import logging
from metainformant.dna import sequences
from metainformant.rna import workflow
from metainformant.protein import parse_fasta

# Setup logger for multi-module pipeline
logger = logging.setup_logger("multiomics_pipeline", level="INFO")

# DNA processing with logging
logger.info("Starting DNA sequence analysis", extra={"module": "dna"})
seqs = sequences.read_fasta("sequences.fasta")
logger.info(f"Loaded {len(seqs)} sequences", extra={"count": len(seqs)})

# RNA processing with logging
logger.info("Starting RNA expression analysis", extra={"module": "rna"})
expression_data = workflow.extract_expression("expression.tsv")
logger.info(f"Extracted expression for {len(expression_data)} genes", 
            extra={"genes": len(expression_data)})

# Protein processing with logging
logger.info("Starting protein analysis", extra={"module": "protein"})
proteins = parse_fasta(Path("proteome.fasta"))
logger.info(f"Loaded {len(proteins)} proteins", extra={"count": len(proteins)})

# Error logging with context
try:
    results = process_data(data)
except Exception as e:
    logger.error("Processing failed", extra={"error": str(e), "module": "processing"})
    raise
```

### Parallel Processing in Cross-Domain Analyses

```python
from metainformant.core import parallel
from metainformant.dna import sequences
from metainformant.protein import parse_fasta
from metainformant.ml import BiologicalClassifier

# Parallel sequence processing
def process_sequence(seq_data):
    seq_id, seq = seq_data
    # Process sequence (e.g., calculate features)
    features = extract_features(seq)
    return seq_id, features

sequences_list = list(sequences.read_fasta("sequences.fasta").items())
# Process sequences in parallel
results = parallel.thread_map(process_sequence, sequences_list, max_workers=4)

# Parallel protein analysis
def analyze_protein(protein_data):
    prot_id, seq = protein_data
    composition = calculate_aa_composition(seq)
    return prot_id, composition

proteins_list = list(parse_fasta(Path("proteome.fasta")).items())
protein_results = parallel.thread_map(analyze_protein, proteins_list, max_workers=8)

# Parallel ML model training
def train_model(fold_data):
    X_train, X_test, y_train, y_test = fold_data
    classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
    classifier.fit(X_train, y_train)
    score = classifier.score(X_test, y_test)
    return score

# Cross-validation folds
cv_folds = [(X_train_fold, X_test_fold, y_train_fold, y_test_fold) 
            for fold in range(5)]
cv_scores = parallel.thread_map(train_model, cv_folds, max_workers=5)
```

### Path Management Across Modules

```python
from metainformant.core import paths
from metainformant.dna import sequences
from metainformant.rna import workflow
from metainformant.quality import analyze_fastq_quality

# Ensure output directory exists
output_dir = paths.join_paths("output", "multiomics", "analysis")
paths.ensure_dir(output_dir)

# Safe path handling for different modules
dna_input = paths.expand_path("data/dna/sequences.fasta")
rna_input = paths.expand_path("data/rna/expression.tsv")
protein_input = paths.join_paths("data", "protein", "proteome.fasta")

# Validate paths before processing
if paths.file_exists(dna_input):
    seqs = sequences.read_fasta(dna_input)
else:
    raise FileNotFoundError(f"DNA input not found: {dna_input}")

# Write outputs with proper path handling
dna_output = paths.join_paths(output_dir, "dna_results.json")
rna_output = paths.join_paths(output_dir, "rna_results.json")
io.write_json(dna_results, dna_output)
io.write_json(rna_results, rna_output)
```

### Complete Multi-Module Workflow Example

```python
from metainformant.core import config, io, logging, paths, parallel
from metainformant.dna import sequences
from metainformant.rna import workflow
from metainformant.multiomics import MultiOmicsData

# 1. Configuration
cfg = config.load_config("config/multiomics_workflow.yaml")
logger = logging.setup_logger("multiomics", level=cfg.logging.level)

# 2. Path setup
output_dir = paths.join_paths(cfg.output.base_dir, "multiomics")
paths.ensure_dir(output_dir)

# 3. Data loading with logging
logger.info("Loading genomic data")
dna_seqs = sequences.read_fasta(paths.expand_path(cfg.input.dna))

logger.info("Loading transcriptomic data")
rna_data = workflow.extract_expression(paths.expand_path(cfg.input.rna))

# 4. Parallel processing
def process_dna_batch(batch):
    return [extract_features(seq) for seq in batch]

dna_batches = create_batches(list(dna_seqs.values()), batch_size=100)
dna_features = parallel.thread_map(process_dna_batch, dna_batches, 
                                    max_workers=cfg.compute.threads)

# 5. Integration
omics_data = MultiOmicsData(
    genomics=pd.DataFrame(dna_features),
    transcriptomics=rna_data
)

# 6. Save results
results_path = paths.join_paths(output_dir, "results.json")
io.write_json(omics_data.to_dict(), results_path)
logger.info("Workflow completed", extra={"output": results_path})
```

### Config-Based Processing (`__init__.py`)
High-level utilities for configuration-driven data processing workflows.

**Key Features:**
- Download and process data based on configuration files
- Config file validation and error reporting
- Sample configuration file generation
- End-to-end workflow orchestration

**Usage:**
```python
from metainformant.core import download_and_process_data, validate_config_file, create_sample_config

# Process data based on configuration
results = download_and_process_data(config_dict)

# Validate configuration files
is_valid, errors = validate_config_file("config.json")
if not is_valid:
    print("Config errors:", errors)

# Create sample configurations
create_sample_config("sample_config.json", "scientific")
```

This ensures consistency and reduces code duplication across the entire codebase.
