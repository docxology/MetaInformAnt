# Core: Logging Framework

The `logging` module provides consistent, structured logging across all METAINFORMANT modules with support for console and file output.

## Functions

### Logger Creation
- **`get_logger(name)`** → `logging.Logger`
  - Get or create a logger with console output
  - Automatic handler setup if none exists
  - Default INFO level formatting

- **`setup_logger(name, log_file=None, level="INFO")`** → `logging.Logger`
  - Create logger with console and optional file output
  - Configurable logging level
  - Automatic directory creation for log files

### Structured Logging
- **`log_with_metadata(logger, message, metadata)`** → `None`
  - Log messages with structured JSON metadata
  - Useful for adding context to log entries

## Usage Examples

### Basic Logging
```python
from metainformant.core import logging

# Simple logger creation
logger = logging.get_logger("metainformant.dna")
logger.info("Starting sequence analysis")
logger.warning("Low quality sequence detected")
logger.error("Failed to parse FASTA file")
```

### Advanced Logger Setup
```python
from metainformant.core import logging

# Logger with file output
logger = logging.setup_logger(
    name="metainformant.rna.pipeline",
    log_file="output/logs/rna_pipeline.log",
    level="DEBUG"
)

logger.info("Pipeline started")
logger.debug("Processing 1000 samples")
logger.warning("Some samples failed quality control")
```

### Structured Logging with Metadata
```python
from metainformant.core import logging

logger = logging.setup_logger("metainformant.analysis")

# Log with structured metadata
analysis_metadata = {
    "sample_count": 500,
    "genes_analyzed": 25000,
    "computation_time": 45.2,
    "memory_peak": "2.1GB"
}

logging.log_with_metadata(
    logger,
    "Differential expression analysis completed",
    analysis_metadata
)
```

## Log Format

All log messages use a consistent format:
```
2024-01-15 14:30:25 | INFO | metainformant.dna.alignment | Alignment completed for 100 sequences
```

Format components:
- **Timestamp**: YYYY-MM-DD HH:MM:SS format
- **Level**: INFO, WARNING, ERROR, DEBUG
- **Logger Name**: Hierarchical module naming
- **Message**: Log message content

## Best Practices

### Logger Naming
Use hierarchical naming that reflects the module structure:
```python
# Good
logger = logging.get_logger("metainformant.dna.alignment")

# Avoid
logger = logging.get_logger("my_script")
```

### Log Levels
- **DEBUG**: Detailed diagnostic information
- **INFO**: General information about program execution
- **WARNING**: Potential issues that don't stop execution
- **ERROR**: Serious problems that may affect results

### Structured Metadata
Use `log_with_metadata` for machine-readable log entries:
```python
# Good for monitoring/analysis
logging.log_with_metadata(logger, "Batch processed", {
    "batch_id": "batch_001",
    "records_processed": 1000,
    "processing_time": 15.3,
    "success_rate": 0.95
})
```

## File Output

When using `setup_logger` with a log file:
- Parent directories are created automatically
- Log files are appended (not overwritten)
- Both console and file get the same formatted output

```python
# Logs to both console and file
logger = logging.setup_logger(
    "metainformant.workflow",
    log_file="output/workflow.log"
)
```

## Integration with Other Modules

The logging module is designed to integrate seamlessly:
```python
from metainformant.core import logging, io

logger = logging.setup_logger("metainformant.export")

try:
    # Some operation
    data = io.load_json("input/data.json")
    logger.info(f"Loaded {len(data)} records")

    # Export operation
    io.dump_json(data, "output/processed.json")
    logger.info("Export completed successfully")

except Exception as e:
    logger.error(f"Export failed: {e}")
    raise
```

## Dependencies

- **Required**: Standard library `logging` module
- **Optional**: None (all functionality uses stdlib)
