# Core Infrastructure: Getting Started

This 5-minute guide demonstrates how to use the core infrastructure components together in a realistic bioinformatics pipeline scenario: downloading sequence data, caching results, processing with parallel execution, and tracking progress with structured logging.

## Prerequisites

```bash
# Install METAINFORMANT in development mode
cd /path/to/MetaInformAnt
uv pip install -e .

# Optional dependencies for full functionality
uv pip install requests pandas pyarrow psycopg2-binary pyyaml tqdm
```

## Scenario: Species Data Processing Pipeline

You need to:
1. Download gene annotation files for multiple species
2. Cache the downloads to avoid re-fetching
3. Parse and validate the files in parallel
4. Store results in a database
5. Generate a summary report

## Complete Example

```python
#!/usr/bin/env python3
"""
Complete pipeline using core infrastructure components.
Run time: ~5 minutes to download and process 3 small files.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from metainformant.core import io, cache, config, paths, logging as core_logging
from metainformant.core.execution import parallel
from metainformant.core.data import db

# ----------------------------------------------------------------------
# Step 1: Initialize logging and directories
# ----------------------------------------------------------------------
logger = core_logging.get_logger("pipeline")
logger.info("Starting species data processing pipeline")

# Ensure output directories exist
output_dir = paths.ensure_directory(Path("output/species_data"))
cache_dir = paths.ensure_directory(Path("cache/downloads"))
log_dir = paths.ensure_directory(Path("logs"))

# ----------------------------------------------------------------------
# Step 2: Define data sources (simulated URLs for demo)
# ----------------------------------------------------------------------
SPECIES_DATA = {
    "homo_sapiens": {
        "url": "https://raw.githubusercontent.com/greenelab/adage/refs/heads/master/adage/data/gene_annotations.tsv",
        "description": "Human gene annotations (ENSG)",
    },
    "mus_musculus": {
        "url": "https://raw.githubusercontent.com/条形码/barcode/refs/heads/main/data/mouse_genes.tsv",
        "description": "Mouse gene annotations",
    },
    "danio_rerio": {
        "url": "https://raw.githubusercontent.com/zf誠信/zebrafish/refs/heads/main/data/zf_genes.tsv",
        "description": "Zebrafish gene annotations",
    },
}

# ----------------------------------------------------------------------
# Step 3: Download with caching (using core.cache)
# ----------------------------------------------------------------------
def download_species_data(species: str, info: dict[str, Any]) -> dict[str, Any]:
    """Download and cache annotation file for a species."""
    logger.info(f"Processing {species}: {info['description']}")

    # Generate cache key from species name
    cache_key = f"annotation/{species}"

    # Check cache first (TTL: 24 hours = 86400 seconds)
    cached = cache.load_cached_json(cache_dir, cache_key, ttl_seconds=86400)
    if cached is not None:
        logger.info(f"✓ Cache hit for {species}")
        return {"species": species, "cached": True, "data": cached}

    # Cache miss: download fresh data
    logger.info(f"↻ Cache miss, downloading {species}...")
    dest_path = cache_dir / f"{species}_annotations.tsv"

    try:
        # Use core.io.download_file (which wraps download_with_progress)
        success = io.download_file(info["url"], dest_path, chunk_size=8192, timeout=30)
        if not success:
            raise RuntimeError(f"Download failed for {species}")

        # Read and parse the file
        with io.open_text_auto(dest_path) as fh:
            lines = fh.readlines()

        # Simple parsing: assume TSV with gene_id and gene_symbol columns
        genes = []
        for line in lines[1:]:  # Skip header
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                genes.append({"id": parts[0], "symbol": parts[1]})

        result_data = {
            "species": species,
            "gene_count": len(genes),
            "genes": genes[:5],  # Store first 5 as sample
        }

        # Cache the result for 24 hours
        cache.cache_json(cache_dir, cache_key, result_data, ttl_seconds=86400)
        logger.info(f"✓ Downloaded and cached {species} ({len(genes)} genes)")

        return {"species": species, "cached": False, "data": result_data}

    except Exception as e:
        logger.error(f"✗ Failed to process {species}: {e}")
        return {"species": species, "error": str(e)}

# ----------------------------------------------------------------------
# Step 4: Parallel execution (using core.parallel)
# ----------------------------------------------------------------------
def run_parallel_downloads() -> list[dict[str, Any]]:
    """Download all species data in parallel (I/O-bound)."""
    logger.info("=" * 60)
    logger.info("PHASE 1: Parallel download with caching")
    logger.info("=" * 60)

    # Prepare items for parallel mapping
    items = list(SPECIES_DATA.items())

    # Use thread_map for I/O-bound downloads (4 workers)
    results = parallel.thread_map(
        lambda item: download_species_data(item[0], item[1]),
        items,
        max_workers=4,
        ordered=True,  # Preserve input order
    )

    # Summarize
    successful = [r for r in results if "error" not in r]
    cached_count = sum(1 for r in successful if r.get("cached"))
    logger.info(f"Completed: {len(successful)}/{len(items)} species processed")
    logger.info(f"  - {cached_count} from cache")
    logger.info(f"  - {len(successful) - cached_count} freshly downloaded")

    return results

# ----------------------------------------------------------------------
# Step 5: Data validation (simulated)
# ----------------------------------------------------------------------
def validate_species_data(result: dict[str, Any]) -> dict[str, Any]:
    """Validate downloaded gene annotation data."""
    if "error" in result:
        return result

    species = result["species"]
    data = result["data"]

    # Check required fields
    if not data.get("genes"):
        return {"species": species, "valid": False, "error": "No genes found"}

    # Check gene ID format (simple heuristic)
    sample_gene = data["genes"][0]
    if not sample_gene.get("id") or not sample_gene.get("symbol"):
        return {"species": species, "valid": False, "error": "Missing gene fields"}

    result["valid"] = True
    return result

def validate_all(data: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Validate all species in parallel."""
    logger.info("=" * 60)
    logger.info("PHASE 2: Parallel validation")
    logger.info("=" * 60)

    validated = parallel.thread_map(validate_species_data, data, max_workers=4)
    valid_count = sum(1 for v in validated if v.get("valid"))
    logger.info(f"Validation: {valid_count}/{len(validated)} passed")

    return validated

# ----------------------------------------------------------------------
# Step 6: Database storage (using core.db)
# ----------------------------------------------------------------------
def store_in_database(data: list[dict[str, Any]]) -> None:
    """Store validated data in PostgreSQL database."""
    logger.info("=" * 60)
    logger.info("PHASE 3: Database storage")
    logger.info("=" * 60)

    try:
        # Configure database connection via environment variables:
        # export PG_HOST="localhost"
        # export PG_PORT="5432"
        # export PG_DATABASE="metainformant"
        # export PG_USER="postgres"
        # export PG_PASSWORD="your_password"

        with db.get_connection() as conn:
            with conn.connect() as db_conn:
                # Create table if not exists
                schema = {
                    "species": "VARCHAR(100) PRIMARY KEY",
                    "gene_count": "INTEGER",
                    "validation_status": "BOOLEAN",
                    "processed_at": "TIMESTAMP DEFAULT NOW()",
                }
                conn.create_table(db_conn, "species_annotations", schema, indexes=["species"])

                # Insert data
                valid_data = [d for d in data if d.get("valid")]
                for item in valid_data:
                    query = """
                        INSERT INTO species_annotations (species, gene_count, validation_status)
                        VALUES (%s, %s, %s)
                        ON CONFLICT (species) DO UPDATE SET
                            gene_count = EXCLUDED.gene_count,
                            validation_status = EXCLUDED.validation_status,
                            processed_at = NOW()
                    """
                    params = (item["species"], item["data"]["gene_count"], True)
                    conn.execute_query(db_conn, query, params, fetch=False)

                logger.info(f"✓ Stored {len(valid_data)} species in database")

    except ImportError:
        logger.warning("psycopg2 not installed – skipping database storage")
    except Exception as e:
        logger.error(f"Database storage failed: {e}")

# ----------------------------------------------------------------------
# Step 7: Generate summary report
# ----------------------------------------------------------------------
def generate_report(
    results: list[dict[str, Any]],
    validated: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Generate JSON summary report of the pipeline run."""
    logger.info("=" * 60)
    logger.info("PHASE 4: Report generation")
    logger.info("=" * 60)

    successful = [r for r in results if "error" not in r]
    errors = [r for r in results if "error" in r]
    valid = [v for v in validated if v.get("valid")]

    report = {
        "pipeline": "species_data_processing",
        "timestamp": time.time(),
        "summary": {
            "total_species": len(SPECIES_DATA),
            "successful_downloads": len(successful),
            "failed_downloads": len(errors),
            "validated": len(valid),
            "from_cache": sum(1 for r in successful if r.get("cached")),
        },
        "details": {
            "errors": errors,
            "validated_species": [v["species"] for v in valid],
        },
    }

    # Atomic write (core.io.dump_json uses atomic=True by default)
    io.dump_json(report, output_path, indent=2)
    logger.info(f"✓ Report saved to {output_path}")

    # Also print to stdout
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    print(f"Total species: {len(SPECIES_DATA)}")
    print(f"Successful: {len(successful)}")
    print(f"From cache: {report['summary']['from_cache']}")
    print(f"Validated: {len(valid)}")
    print(f"Errors: {len(errors)}")
    if errors:
        print("\nFailed species:")
        for err in errors:
            print(f"  - {err['species']}: {err['error']}")
    print(f"\nReport: {output_path}")
    print("=" * 60)

# ----------------------------------------------------------------------
# Main pipeline orchestration
# ----------------------------------------------------------------------
def main() -> None:
    """Run the complete pipeline."""
    start_time = time.time()

    try:
        # Phase 1: Parallel downloads with caching
        download_results = run_parallel_downloads()

        # Phase 2: Validation
        validated = validate_all(download_results)

        # Phase 3: Database storage (optional)
        store_in_database(validated)

        # Phase 4: Report
        report_path = output_dir / "pipeline_report.json"
        generate_report(download_results, validated, report_path)

        elapsed = time.time() - start_time
        logger.info(f"✓ Pipeline completed in {elapsed:.1f}s")

    except KeyboardInterrupt:
        logger.warning("Pipeline interrupted by user")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise

if __name__ == "__main__":
    # Configure logging from environment (CORE_LOG_LEVEL=DEBUG for verbose)
    core_logging.configure_logging_from_env(default_level="INFO")
    main()
```

## Key Components Demonstrated

| Component | Functions Used | Purpose |
|-----------|---------------|---------|
| `core.io` | `download_file()`, `open_text_auto()`, `dump_json()` | File I/O and downloads |
| `core.cache` | `cache_json()`, `load_cached_json()` | TTL-based caching |
| `core.paths` | `ensure_directory()` | Directory creation |
| `core.utils.logging` | `get_logger()`, `configure_logging_from_env()` | Structured logging |
| `core.parallel` | `thread_map()` | Parallel I/O execution |
| `core.data.db` | `get_connection()`, `execute_query()` | Database operations |
| `core.config` | (Not shown) `load_mapping_from_file()` | Config file loading |

## Configuration Files

You can extend the pipeline with external config:

**config/pipeline.yaml**:
```yaml
pipeline:
  name: "Species Data Processing"
  cache_ttl_seconds: 86400
  max_download_workers: 4

data_sources:
  homo_sapiens:
    url: "https://example.com/human_genes.tsv"
    description: "Human gene annotations"
  mus_musculus:
    url: "https://example.com/mouse_genes.tsv"
    description: "Mouse gene annotations"
```

Load with:
```python
from metainformant.core import config

pipeline_config = config.load_mapping_from_file("config/pipeline.yaml")
max_workers = pipeline_config["pipeline"]["max_download_workers"]
```

## Environment Variables

Set these in your shell or `.env` file:

```bash
# Logging
export CORE_LOG_LEVEL=INFO   # DEBUG, INFO, WARNING, ERROR

# Database (optional)
export PG_HOST=localhost
export PG_PORT=5432
export PG_DATABASE=metainformant
export PG_USER=postgres
export PG_PASSWORD=secret

# General overrides
export AK_THREADS=4          # Override thread pool size
export AK_WORK_DIR=/data     # Working directory
export AK_LOG_DIR=/var/log   # Log directory
```

## Running the Example

```bash
# 1. Make the script executable
chmod +x pipeline_example.py

# 2. Run with INFO logging
./pipeline_example.py

# 3. Or with DEBUG logging for verbose output
CORE_LOG_LEVEL=DEBUG ./pipeline_example.py

# 4. View logs
tail -f logs/pipeline.log

# 5. Check cache size
du -sh cache/downloads/

# 6. Force re-download (clear cache)
rm -rf cache/downloads/
./pipeline_example.py
```

## Expected Output

```
2026-04-26 19:15:01 | INFO | pipeline | Starting species data processing pipeline
2026-04-26 19:15:01 | INFO | pipeline | ============================================================
2026-04-26 19:15:01 | INFO | pipeline | PHASE 1: Parallel download with caching
2026-04-26 19:15:01 | INFO | pipeline | ============================================================
2026-04-26 19:15:02 | INFO | pipeline | ✓ Cache miss, downloading homo_sapiens...
2026-04-26 19:15:03 | INFO | pipeline | ✓ Downloaded and cached homo_sapiens (20456 genes)
...
2026-04-26 19:15:08 | INFO | pipeline | ✓ Pipeline completed in 7.2s
```

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'metainformant'"
**Fix**: Install in editable mode: `uv pip install -e .`

### Issue: Database connection fails
**Fix**: Set `PG_HOST`, `PG_DATABASE`, `PG_USER`, `PG_PASSWORD` environment variables, or skip database phase by catching `ImportError`.

### Issue: Downloads are slow
**Fix**: Increase `max_workers` in `parallel.thread_map()`, but respect server rate limits. Use `rate_limited_map()` if needed.

### Issue: Cache not expiring
**Fix**: The `ttl_seconds` parameter controls expiration. Use `cache.clear_cache_dir()` to manually clear.

### Issue: "psycopg2 not available"
**Fix**: Install: `uv pip install psycopg2-binary` (development) or `psycopg2` (production).

## Next Steps

- Read detailed component guides in `[I/O Operations](./io.md)`, `[Caching](./cache.md)`, `[Parallel](./parallel.md)`
- Review [Architecture](./ARCHITECTURE.md) for system design
- Explore [API Reference](./SPEC.md) for complete function signatures
- See real-world usage in `examples/core/` directory
