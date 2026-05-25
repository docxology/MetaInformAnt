# Core: Database Integration

The `db` module provides PostgreSQL database connectivity with secure credential management, connection pooling, transaction handling, and robust error recovery. It serves as the primary data persistence layer for METAINFORMANT workflows.

## Purpose

Bioinformatics analyses generate and require large-scale structured data storage:
- **Sample metadata**: Thousands to millions of records with complex relationships
- **Analysis results**: Quantification data, variant calls, association statistics
- **Reference annotations**: Gene models, pathway memberships, ontology terms
- **Workflow provenance**: Execution logs, parameter tracking, data lineage

The `db` module provides production-ready database integration with connection pooling, batch operations, and transaction management to handle these demanding workloads efficiently and securely.

## Design Principles

### 1. **Connection Pooling**
Uses `psycopg2.pool.ThreadedConnectionPool` to maintain a pool of reusable connections. This avoids the overhead of establishing a new TCP connection and PostgreSQL backend process for each query, which can take 100-500ms per connection. Pools can be tuned with `min_connections` and `max_connections` to match workload characteristics.

### 2. **Context Manager Support**
The `PostgresConnection` class and `get_connection()` context manager ensure connections are always properly released back to the pool, even when exceptions occur. This prevents connection leaks that would otherwise exhaust the pool and block subsequent queries.

### 3. **Parameterized Queries Only**
All query execution uses Python parameter placeholders (`%s`) with separate parameter tuples. This prevents SQL injection and allows PostgreSQL to cache query plans. String interpolation or concatenation for query building is explicitly avoided.

### 4. **Transaction Integrity**
`execute_transaction()` groups multiple queries into an atomic unit with automatic rollback on failure. Bulk operations use `executemany()` within a single transaction for both correctness and performance (reducing WAL overhead).

### 5. **Graceful Degradation**
The module detects optional `psycopg2` availability and raises informative `ImportError` messages only when database functionality is actually invoked. This allows METAINFORMANT to be installed in environments without PostgreSQL for other uses.

### 6. **Security First**
Connection parameters are sanitized to strip dangerous SQL keywords and characters. Passwords are never logged. All database credentials must come from environment variables, never from code or config files that might be committed to version control.

### 7. **Typed Results**
Queries return `List[Dict[str, Any]]` using `RealDictCursor` so column names are accessible as dictionary keys. This is more maintainable than positional tuple access, especially for queries with 10+ columns.

## Module Organization

The `db` module is in `src/metainformant/core/data/db.py`. It exposes:

**Classes**:
- `PostgresConnection` — Main connection pool wrapper with query execution methods
- `DownloadResult` (re-export from download module for convenience)

**Functions**:
- `get_db_client()` — Create a `PostgresConnection` instance with environment-based config
- `build_postgres_url()` — Build connection URL string (for SQLAlchemy, etc.)
- `sanitize_connection_params()` — Clean user-provided connection parameters
- `get_connection()` — Context manager for automatic connection cleanup

**Protocols**:
- `DownloadHandler` — Interface for pluggable download handlers (HTTP, FTP, file)

## Configuration

### Environment Variables

All database credentials are read from environment variables. The module checks two variable naming schemes:

**Long-form PG_* variables** (preferred):
```bash
export PG_HOST="db.example.com"
export PG_PORT="5432"
export PG_DATABASE="metainformant"
export PG_USER="analysis_user"
export PG_PASSWORD="secure_password_here"
```

**Short-form DB_* variables** (fallback):
```bash
export DB_NAME="metainformant"
export DB_USER="analysis_user"
export DB_PASSWORD="secure_password_here"
# PG_HOST still used for host
```

If both are set, PG_* takes precedence. The `get_db_client()` function reads these variables at call time.

### Optional Configuration

The `PostgresConnection` constructor accepts additional parameters:
- `min_connections` (default: 1) — Minimum number of connections maintained in pool
- `max_connections` (default: 10) — Maximum number of connections allowed
- All standard PostgreSQL connection parameters (host, port, database, user, password)

## API Reference

### Classes

#### `PostgresConnection`

Wrapper around `psycopg2.pool.ThreadedConnectionPool` providing high-level database operations.

**Constructor**:
```text
PostgresConnection(
    host: str = "localhost",
    port: int = 5432,
    database: str = "metainformant",
    user: str = "postgres",
    password: str = "",
    min_connections: int = 1,
    max_connections: int = 10,
) -> PostgresConnection
```

**Raises**:
- `ImportError`: If `psycopg2` is not installed
- `RuntimeError`: If connection pool cannot be initialized

**Methods**:

##### `connect() -> psycopg2.extensions.connection`

Get a connection from the pool.

**Returns**: Database connection object

**Raises**:
- `RuntimeError`: If pool is exhausted or connection fails

**Example**:
```python
conn_wrapper = PostgresConnection(host="localhost", database="mydb")
conn = conn_wrapper.connect()
try:
    # Use conn
    pass
finally:
    conn_wrapper.release_connection(conn)
```

##### `release_connection(conn: psycopg2.extensions.connection) -> None`

Return a connection to the pool. This is called automatically by context managers. You generally don't call this directly—use `with conn.connect() as db_conn:` pattern instead.

##### `execute_query(
    conn: psycopg2.extensions.connection,
    query: str,
    params: Tuple = (),
    fetch: bool = True,
) -> List[Dict[str, Any]]`

Execute a single SQL query and optionally fetch results.

**Parameters**:
- `conn`: Database connection from `connect()`
- `query`: SQL query string with `%s` placeholders for parameters
- `params`: Tuple of parameters to substitute (empty tuple for no parameters)
- `fetch`: If `True`, fetch and return all results; if `False`, execute without fetch (for INSERT/UPDATE/DELETE)

**Returns**: List of dictionaries (each dict = one row) if `fetch=True`; empty list otherwise

**Raises**:
- `RuntimeError`: If query execution fails (with automatic rollback)

**Example**:
```python
conn_wrapper = get_db_client()
with conn_wrapper.connect() as conn:
    # Query with parameters
    results = conn_wrapper.execute_query(
        conn,
        "SELECT * FROM samples WHERE species = %s AND quality > %s",
        ("human", 20.0)
    )
    for row in results:
        print(row["sample_id"], row["species"])

    # Insert without fetch
    conn_wrapper.execute_query(
        conn,
        "INSERT INTO logs (message, level) VALUES (%s, %s)",
        ("Pipeline started", "INFO"),
        fetch=False
    )
```

##### `execute_transaction(
    conn: psycopg2.extensions.connection,
    queries: List[Tuple[str, Tuple]],
) -> List[List[Dict[str, Any]]]`

Execute multiple queries atomically within a single transaction.

**Parameters**:
- `conn`: Database connection from `connect()`
- `queries`: List of `(query_string, params_tuple)` pairs

**Returns**: List of result lists—each element is the result list from the corresponding query

**Raises**:
- `RuntimeError`: If any query fails (triggers rollback of entire transaction)

**Example**:
```python
queries = [
    ("INSERT INTO samples (id, name) VALUES (%s, %s)", ("S001", "Sample 1")),
    ("INSERT INTO measurements (sample_id, value) VALUES (%s, %s)", ("S001", 42.5)),
    ("UPDATE pipeline_status SET status = 'completed' WHERE run_id = %s", ("run_123",)),
]

with get_connection(database="mydb") as conn_wrapper:
    with conn_wrapper.connect() as conn:
        results = conn_wrapper.execute_transaction(conn, queries)
        # All three queries succeed or all are rolled back
```

##### `bulk_insert(
    conn: psycopg2.extensions.connection,
    table: str,
    columns: List[str],
    data: List[Tuple],
    batch_size: int = 1000,
) -> int`

Efficiently insert many rows using batched `executemany()`.

**Parameters**:
- `conn`: Database connection from `connect()`
- `table`: Target table name (must be trusted/sanitized—no user input)
- `columns`: List of column names in order
- `data`: List of tuples, each tuple containing values for one row
- `batch_size`: Number of rows per batch (default 1000; tune for your workload)

**Returns**: Total number of rows inserted

**Raises**:
- `RuntimeError`: If insert fails (with rollback)

**Performance notes**:
- Batch size of 1000–5000 is typically optimal for PostgreSQL
- Larger batches can cause memory pressure; smaller batches increase round-trips
- Uses a single prepared statement reused across all batches in the call

**Example**:
```python
# Prepare 100,000 variant records
columns = ["chrom", "pos", "ref", "alt", "sample_id", "quality"]
rows = []
for variant in variants:
    rows.append((
        variant.chrom,
        variant.pos,
        variant.ref,
        variant.alt,
        variant.sample_id,
        variant.quality,
    ))

with get_connection() as conn_wrapper:
    with conn_wrapper.connect() as conn:
        inserted = conn_wrapper.bulk_insert(
            conn,
            table="variants",
            columns=columns,
            data=rows,
            batch_size=5000
        )
        print(f"Inserted {inserted} variant records")
```

##### `create_table(
    conn: psycopg2.extensions.connection,
    table_name: str,
    schema: Dict[str, str],
    indexes: Optional[List[str]] = None,
    if_not_exists: bool = True,
) -> None`

Create a database table with optional indexes.

**Parameters**:
- `conn`: Database connection from `connect()`
- `table_name`: Name of table to create
- `schema`: Dictionary mapping column names to PostgreSQL type strings (e.g., `"SERIAL PRIMARY KEY"`, `"VARCHAR(255)"`, `"FLOAT"`, `"TIMESTAMP"`)
- `indexes`: Optional list of column names to create B-tree indexes on
- `if_not_exists`: If `True`, use `CREATE TABLE IF NOT EXISTS` to avoid error on re-creation

**Raises**:
- `RuntimeError`: If table creation fails

**Example**:
```python
schema = {
    "id": "SERIAL PRIMARY KEY",
    "sample_id": "VARCHAR(100) NOT NULL",
    "species": "VARCHAR(50)",
    "expression_mean": "FLOAT",
    "created_at": "TIMESTAMP DEFAULT NOW()",
}

indexes = ["sample_id", "species"]

with get_connection() as conn_wrapper:
    with conn_wrapper.connect() as conn:
        conn_wrapper.create_table(
            conn,
            table_name="gene_expression",
            schema=schema,
            indexes=indexes,
            if_not_exists=True
        )
```

##### `table_exists(
    conn: psycopg2.extensions.connection,
    table_name: str,
) -> bool`

Check whether a table exists in the current database.

**Parameters**:
- `conn`: Database connection from `connect()`
- `table_name`: Table name to check

**Returns**: `True` if table exists, `False` otherwise

**Example**:
```python
with get_connection() as conn_wrapper:
    with conn_wrapper.connect() as conn:
        if not conn_wrapper.table_exists(conn, "samples"):
            print("Table not found—creating...")
            # Create table
```

##### `close() -> None`

Close the connection pool and release all connections. Typically called automatically by context manager, but can be called manually for cleanup.

**Example**:
```python
pool = get_db_client()
# ... use pool ...
pool.close()  # Clean shutdown
```

### Standalone Functions

#### `build_postgres_url(host, port, database, user, password) -> str`

Construct a PostgreSQL connection URL in the format:
```
postgresql://user:***@host:port/database
```

Password is masked as `***` in the returned string for safe logging.

**Parameters**:
- `host`: Database hostname or IP
- `port`: Port number (typically 5432)
- `database`: Database name
- `user`: Username for authentication
- `password`: Password (included but masked in output)

**Returns**: Connection URL string with password hidden

**Example**:
```python
url = build_postgres_url(
    host="db.cluster.example.com",
    port=5432,
    database="omics",
    user="analyst",
    password="s3cr3t!"
)
print(url)  # postgresql://analyst:***@db.cluster.example.com:5432/omics
# Use with SQLAlchemy:
# engine = create_engine(url)
```

#### `get_db_client(**kwargs) -> PostgresConnection`

Factory function that reads environment variables and creates a `PostgresConnection` instance.

**Parameters**: Any `PostgresConnection` constructor parameters can be passed as keyword arguments to override environment variables.

**Returns**: Configured `PostgresConnection` instance

**Raises**:
- `ImportError`: If psycopg2 unavailable
- `RuntimeError`: If required environment variables missing

**Environment variable priority**:
1. Explicit kwargs passed to function
2. `PG_*` environment variables
3. `DB_*` environment variables (for backward compatibility)

**Example**:
```python
# Uses environment variables PG_HOST, PG_PORT, PG_DATABASE, PG_USER, PG_PASSWORD
client = get_db_client()

# Override with explicit parameters
client = get_db_client(
    host="override.host",
    database="override_db",
    max_connections=20
)
```

#### `sanitize_connection_params(params: Dict[str, Any]) -> Dict[str, Any]`

Remove potentially dangerous characters and SQL keywords from connection parameter values to mitigate injection attacks. This is a basic sanitization layer; parameterized queries remain the primary defense.

**Parameters**:
- `params`: Dictionary of connection parameters (e.g., `{"host": "...", "user": "..."}`)

**Returns**: New dictionary with sanitized string values

**Sanitization rules**:
- Single quotes (`'`) removed
- Semicolons (`;`) removed (statement separator)
- Double-dash (`--`) comment sequence removed
- SQL keywords removed from values: `DROP`, `DELETE`, `UPDATE`, `INSERT`, `ALTER`, `CREATE`, `TRUNCATE`, `EXEC`, `EXECUTE`
- Non-string values passed through unchanged

**Important**: This function is a belt-and-suspenders measure. The real protection is using parameterized queries (`%s` placeholders) everywhere.

**Example**:
```python
raw = {
    "host": "localhost; DROP TABLE users;--",
    "user": "admin'; DELETE FROM samples;--",
    "password": "p@ss",
}
safe = sanitize_connection_params(raw)
print(safe["host"])  # "localhost DROP TABLE users"
print(safe["user"])   # "admin DELETE FROM samples"
# Dangerous characters removed, but still use parameterized queries!
```

#### `get_connection(**kwargs) -> Iterator[PostgresConnection]`

Context manager that yields a `PostgresConnection` and ensures cleanup.

**Parameters**: Same as `get_db_client()`

**Yields**: `PostgresConnection` instance

**Example**:
```python
from metainformant.core import db

with db.get_connection(database="analysis") as conn_wrapper:
    with conn_wrapper.connect() as conn:
        results = conn_wrapper.execute_query(conn, "SELECT COUNT(*) FROM genes")
        count = results[0]["count"]
        print(f"Database has {count} genes")
# Connection automatically returned to pool
```

### Protocols

#### `DownloadHandler` (Protocol)

Interface for pluggable download strategy implementations.

**Methods**:
- `get_file_size(url: str, *, timeout: int = 60) -> int | None`
- `download(url: str, dest_path: str | Path, **kwargs) -> DownloadResult`

Built-in implementations:
- `HTTPDownloadHandler` — HTTP/HTTPS with Range resume support
- `FTPDownloadHandler` — FTP (no resume)
- `FileDownloadHandler` — Local file:// URLs

#### `ProgressCallback` (Protocol)

Callback signature for progress reporting during downloads:
```python
def progress_cb(*, bytes_written: int, total_bytes: int | None) -> None: ...
```

## Usage Examples

### Basic Database Connection and Query

```python
from metainformant.core import db

def count_samples_by_species() -> dict[str, int]:
    """Count samples grouped by species."""
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            results = conn_wrapper.execute_query(
                conn,
                """
                SELECT species, COUNT(*) as count
                FROM samples
                GROUP BY species
                ORDER BY count DESC
                """
            )
            return {row["species"]: row["count"] for row in results}

counts = count_samples_by_species()
for species, count in counts.items():
    print(f"{species}: {count} samples")
```

### Transaction with Multiple Statements

```python
from metainformant.core import db

def record_analysis_run(samples: list[str], parameters: dict) -> int:
    """Record a complete analysis run with samples and parameters."""
    queries = [
        # Insert run metadata
        (
            "INSERT INTO analysis_runs (started_at, parameters) "
            "VALUES (NOW(), %s) RETURNING run_id",
            (json.dumps(parameters),)
        ),
        # Insert sample associations
        (
            "INSERT INTO run_samples (run_id, sample_id) "
            "SELECT %s, s.id FROM samples s WHERE s.sample_id = ANY(%s)",
            # run_id filled in after first query, handled in real code with RETURNING
            # This example simplified for clarity—see real implementation using CTEs
        ),
    ]

    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            results = conn_wrapper.execute_transaction(conn, queries)
            run_id = results[0][0]["run_id"]
    return run_id
```

### Bulk Insert with Error Handling

```python
from metainformant.core import db
import json

def store_variant_calls(vcf_path: str, run_id: int) -> int:
    """Parse VCF and bulk-insert variant calls into database."""
    variants = []

    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            variants.append((
                chrom, int(pos), vid, ref, alt,
                float(qual) if qual != "." else None,
                filt, info, run_id
            ))

    columns = ["chrom", "pos", "variant_id", "ref", "alt",
               "quality", "filter", "info", "run_id"]

    try:
        with db.get_connection() as conn_wrapper:
            with conn_wrapper.connect() as conn:
                inserted = conn_wrapper.bulk_insert(
                    conn,
                    table="variants",
                    columns=columns,
                    data=variants,
                    batch_size=5000
                )
                print(f"Inserted {inserted} variants")
                return inserted
    except Exception as e:
        print(f"Failed to insert variants: {e}")
        raise
```

### Connection Pool Tuning

```python
from metainformant.core import db

# High-throughput pipeline: larger pool
pool = db.PostgresConnection(
    host="db.prod.example.com",
    database="omics",
    user="pipeline",
    password=os.environ["PIPELINE_DB_PASSWORD"],
    min_connections=5,
    max_connections=50,  # Handle many concurrent workers
)

# Use throughout pipeline
# ...
pool.close()
```

### Safe Connection with Environment Overrides

```python
import os
from metainformant.core import db, config

# Load database config from environment
pg_config = config.load_postgres_config_from_env(prefix="PG")
if pg_config:
    conn_wrapper = db.PostgresConnection(
        host=pg_config.host,
        port=pg_config.port,
        database=pg_config.database,
        user=pg_config.user,
        password=pg_config.password,
        max_connections=20,
    )
else:
    raise RuntimeError(
        "PostgreSQL configuration not found. Set PG_HOST, PG_DATABASE, "
        "PG_USER, PG_PASSWORD environment variables."
    )
```

### Context Manager Pattern for Reuse

```python
from contextlib import contextmanager
from metainformant.core import db

@contextmanager
def get_db_cursor():
    """Yield a cursor for quick one-liners."""
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                yield cur

# Usage
with get_db_cursor() as cur:
    cur.execute("SELECT * FROM samples WHERE id = %s", (sample_id,))
    row = cur.fetchone()
    # Connection returned to pool on context exit
```

## Integration with Core Modules

### With I/O Module

```python
from metainformant.core import db, io
import json

def export_query_to_json(query: str, output_path: str | Path) -> int:
    """Execute query and write results to JSON file."""
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            results = conn_wrapper.execute_query(conn, query)

    # Use io module for atomic JSON write
    io.dump_json(results, output_path, atomic=True, indent=2)
    return len(results)

# Usage
count = export_query_to_json(
    "SELECT * FROM gene_expressions WHERE p_value < 0.05",
    "output/significant_genes.json"
)
print(f"Exported {count} records")
```

### With Configuration Module

```python
from metainformant.core import db, config

def setup_database_from_config(config_path: str) -> None:
    """Create tables based on YAML configuration."""
    cfg = config.load_mapping_from_file(config_path)

    db_config = cfg["database"]
    schema_config = cfg["schema"]

    # Create connection
    with db.get_connection(
        host=db_config["host"],
        database=db_config["name"],
        user=db_config["user"],
        password=db_config["password"],
    ) as conn_wrapper:
        with conn_wrapper.connect() as conn:
            # Create each table from config
            for table_name, table_def in schema_config["tables"].items():
                conn_wrapper.create_table(
                    conn,
                    table_name=table_name,
                    schema=table_def["columns"],
                    indexes=table_def.get("indexes"),
                    if_not_exists=True
                )
            print(f"Created {len(schema_config['tables'])} tables")
```

### With Logging Module

```python
from metainformant.core import db
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

def query_with_logging(query: str, params: tuple = ()) -> list[dict]:
    """Execute query with detailed logging for audit trail."""
    start_time = time.time()

    try:
        with db.get_connection() as conn_wrapper:
            with conn_wrapper.connect() as conn:
                results = conn_wrapper.execute_query(conn, query, params)

        elapsed = time.time() - start_time
        logger.info(
            "Query executed",
            extra={
                "query": query[:100],  # Truncate long queries
                "params": str(params)[:100],
                "rows_returned": len(results),
                "elapsed_seconds": elapsed,
            }
        )
        return results
    except Exception as e:
        logger.error("Query failed", extra={"query": query, "error": str(e)})
        raise
```

### In a Workflow Pipeline

```python
from metainformant.core import io
from metainformant.core.data import db
from metainformant.core.execution import parallel
from typing import List

def process_and_store_batch(
    batch: List[dict],
    conn_wrapper: db.PostgresConnection
) -> int:
    """Process a batch of records and store in database."""
    # Transform batch
    processed = [
        (r["sample_id"], r["gene_id"], r["expression"])
        for r in batch
    ]

    with conn_wrapper.connect() as conn:
        return conn_wrapper.bulk_insert(
            conn,
            table="expression_values",
            columns=["sample_id", "gene_id", "expression"],
            data=processed,
            batch_size=1000
        )

def bulk_import_pipeline(
    input_file: str,
    max_workers: int = 4
) -> None:
    """Parallel processing with database batch inserts."""
    # Read input data
    records = io.read_jsonl(input_file)

    with db.get_connection(max_connections=10) as conn_wrapper:
        # Process in parallel batches
        def process_batch(batch: List[dict]) -> int:
            return process_and_store_batch(batch, conn_wrapper)

        # Use parallel.batch to chunk and process
        results = parallel.parallel_batch(
            process_batch,
            records,
            batch_size=500,
            max_workers=max_workers
        )
        total = sum(results)
        print(f"Total records imported: {total}")
```

## Error Handling

### Connection Failures

**Symptom**: `RuntimeError: Database connection failed: ...` when calling `get_db_client()` or `conn.connect()`.

**Common causes and fixes**:

1. **PostgreSQL server not running**
   ```bash
   sudo systemctl status postgresql
   sudo systemctl start postgresql
   ```

2. **Network/firewall blocking connection**
   ```bash
   telnet db.example.com 5432
   # Or:
   nc -zv db.example.com 5432
   ```

3. **Invalid credentials**
   ```bash
   # Test with psql first
   PGPASSWORD=yourpass psql -h localhost -U youruser -d yourdb -c "SELECT 1"
   ```

4. **Database doesn't exist**
   ```bash
   # Connect to postgres database to create new one
   PGPASSWORD=pass psql -h localhost -U postgres -d postgres -c "CREATE DATABASE mydb"
   ```

5. **User lacks CONNECT privilege**
   ```sql
   GRANT CONNECT ON DATABASE mydb TO analysis_user;
   ```

### Query Failures

**Symptom**: `RuntimeError: Database query failed: ...` from `execute_query()`.

**Common causes**:

1. **Syntax error** — Check query with `EXPLAIN` or `psql` first
2. **Table/column doesn't exist** — Verify schema with `\dt` and `\d table_name` in psql
3. **Type mismatch** — Ensure parameter types match column types
4. **Constraint violation** — Check NOT NULL, UNIQUE, FOREIGN KEY constraints
5. **Deadlock** — Queries may need retry logic or reordering

**Recovery pattern**:
```python
import time
from metainformant.core import db

def robust_query(query: str, params: tuple, max_retries: int = 3):
    for attempt in range(max_retries):
        try:
            with db.get_connection() as conn_wrapper:
                with conn_wrapper.connect() as conn:
                    return conn_wrapper.execute_query(conn, query, params)
        except RuntimeError as e:
            if "deadlock" in str(e).lower() and attempt < max_retries - 1:
                wait = 0.5 * (2 ** attempt)  # Exponential backoff
                print(f"Deadlock detected, retrying in {wait}s...")
                time.sleep(wait)
            else:
                raise
```

### Transaction Rollbacks

**Symptom**: `psycopg2.InternalError: current transaction is aborted, commands ignored until end of transaction block`.

**Cause**: An error occurred within a transaction, marking it as failed. All subsequent commands in that transaction will fail until `ROLLBACK`.

**Fix**: Ensure you're using `execute_transaction()` for multi-statement operations, or manually rollback on error:
```python
try:
    with conn_wrapper.connect() as conn:
        conn_wrapper.execute_query(conn, "INSERT INTO a VALUES (1)")
        conn_wrapper.execute_query(conn, "INSERT INTO b VALUES (2)")  # If this fails
except Exception:
    conn.rollback()  # Explicit rollback
    raise
```

### Pool Exhaustion

**Symptom**: `psycopg2.pool.PoolError: connection pool exhausted`.

**Cause**: All connections in pool are in use and `max_connections` reached.

**Fixes**:
1. Increase `max_connections` in pool configuration
2. Ensure connections are released (context managers or explicit `release_connection()`)
3. Check for long-running queries—consider indexes or query optimization
4. Reduce `max_workers` in parallel processing if each worker holds a DB connection

**Monitoring**:
```python
pool = get_db_client(max_connections=20)
print(f"Pool stats: {pool._pool._used} of {pool._pool.maxconn} connections in use")
```

### Missing psycopg2

**Symptom**: `ImportError: psycopg2 is required for database operations`.

**Fix**:
```bash
# Install with uv
uv add psycopg2

# Or for binary package (no compilation needed):
uv add psycopg2-binary
```

**Note**: `psycopg2-binary` is fine for development and many production uses. For maximum performance and control, build from source with `psycopg2` (requires libpq and Python headers).

## Performance Considerations

### Connection Pool Sizing

**Too small** (e.g., 2 connections for 8 workers):
- Workers block waiting for connections
- Underutilizes CPU and I/O

**Too large** (e.g., 100 connections for 4 workers):
- Wastes memory (each connection ~10MB)
- Increases context switching overhead
- Can overwhelm PostgreSQL (each connection is a backend process)

**Rule of thumb**:
- For parallel processing: `pool_max = workers * 1.5` (workers may hold connections briefly)
- For long-lived connection pools: `pool_max = min(2 * CPU_cores + 1, 50)`
- Monitor: `SELECT COUNT(*) FROM pg_stat_activity;`

### Bulk Operations

**INSERT patterns**:

| Pattern | Performance | When to use |
|---------|-------------|-------------|
| Single `INSERT` per row | Very slow (1000+ rows) | Never for bulk |
| `executemany()` in transaction | Good (1000–10K rows) | Moderate batch sizes |
| `COPY FROM STDIN` | Excellent (100K+ rows) | Very large imports |

Example using `COPY`:
```python
import io
import csv

def copy_from_csv(conn, table: str, csv_path: str) -> int:
    with open(csv_path) as f:
        with conn.cursor() as cur:
            cur.copy_expert(
                f"COPY {table} FROM STDIN WITH CSV HEADER",
                f
            )
    # Note: COPY is faster but bypasses ORM and triggers
```

### Query Optimization

**Indexes**: Essential for WHERE, JOIN, ORDER BY:
```sql
CREATE INDEX idx_samples_species ON samples(species);
CREATE INDEX idx_variants_chrom_pos ON variants(chrom, pos);
```

**EXPLAIN ANALYZE**: Profile slow queries:
```python
results = conn_wrapper.execute_query(
    conn,
    "EXPLAIN ANALYZE SELECT * FROM variants WHERE chrom = 'chr1' AND pos > 100000"
)
print(results)
```

**Connection-side caching**: Consider caching reference data:
```python
from functools import lru_cache

@lru_cache(maxsize=1)
def get_species_list():
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            results = conn_wrapper.execute_query(conn, "SELECT DISTINCT species FROM samples")
            return [row["species"] for row in results]
```

### Prepared Statements

If you execute the same query thousands of times with different parameters, prepare it once:
```python
with db.get_connection() as conn_wrapper:
    with conn_wrapper.connect() as conn:
        cur = conn.cursor()
        cur.execute("PREPARE myplan AS SELECT * FROM samples WHERE species = $1")
        for species in species_list:
            cur.execute("EXECUTE myplan", (species,))
            # ...
```

### Keep Transactions Short

Long-running transactions:
- Prevent vacuum from cleaning up old rows
- Increase deadlock probability
- Hold locks longer

**Pattern**:
```python
# BAD: Long transaction spanning user interaction
conn = pool.getconn()
cur = conn.cursor()
cur.execute("SELECT ...")
# User thinks for 5 minutes...
cur.execute("UPDATE ...")  # Locks held for 5+ minutes

# GOOD: Short, focused transaction
with conn_wrapper.connect() as conn:
    results = conn_wrapper.execute_query(conn, "SELECT ...")
# Transaction committed/rolled back immediately
# Do user thinking outside transaction
```

### Monitoring

Use PostgreSQL's built-in views:
```python
# Active queries
results = conn_wrapper.execute_query(
    conn,
    "SELECT pid, query, state, query_start FROM pg_stat_activity WHERE state = 'active'"
)

# Table sizes
results = conn_wrapper.execute_query(
    conn,
    """
    SELECT
        schemaname, tablename,
        pg_size_pretty(pg_total_relation_size(schemaname||'.'||tablename)) as size
    FROM pg_tables
    WHERE schemaname NOT IN ('pg_catalog', 'information_schema')
    ORDER BY pg_total_relation_size(schemaname||'.'||tablename) DESC
    """
)
```

## Security Notes

### Credential Management

**NEVER**:
- Commit credentials to version control
- Hard-code passwords in Python files
- Store in plain-text config files in the repository

**ALWAYS**:
- Use environment variables
- Consider a secrets manager (HashiCorp Vault, AWS Secrets Manager) for production
- Rotate credentials periodically
- Use separate database users with minimal required privileges

**Example with python-dotenv** (local development):
```bash
# .env file (NOT committed)
PG_HOST=localhost
PG_PORT=5432
PG_DATABASE=metainformant_dev
PG_USER=dev_user
PG_PASSWORD=dev_pass
```

```python-snippet
from dotenv import load_dotenv
load_dotenv()  # Loads .env into os.environ

from metainformant.core import db
# Now get_db_client() picks up loaded env vars
```

### Principle of Least Privilege

Don't use superuser for application queries:
```sql
-- Create restricted user
CREATE USER analysis_user WITH PASSWORD 'secure_password';

-- Grant specific privileges only
GRANT CONNECT ON DATABASE metainformant TO analysis_user;
GRANT USAGE ON SCHEMA public TO analysis_user;
GRANT SELECT, INSERT, UPDATE ON samples TO analysis_user;
GRANT SELECT, INSERT ON variants TO analysis_user;

-- Don't grant:
-- - SUPERUSER
-- - CREATE on database (unless you need to create tables)
-- - DROP, ALTER, TRUNCATE (unless absolutely needed)
```

### Parameterized Queries

**ALWAYS** use `%s` placeholders:
```python
# GOOD — parameterized
cur.execute("SELECT * FROM samples WHERE id = %s", (sample_id,))

# BAD — string interpolation (SQL injection risk)
cur.execute(f"SELECT * FROM samples WHERE id = {sample_id}")
cur.execute("SELECT * FROM samples WHERE id = '%s'" % sample_id)
```

### SSL/TLS Connections

For production deployments across networks, use SSL:
```python
client = PostgresConnection(
    host="db.example.com",
    database="mydb",
    user="user",
    password="secret",
    sslmode="require",  # In connection string or as parameter
)
# Or via connection parameter:
conn = psycopg2.connect(
    "postgresql://...",
    sslmode="verify-full",
    sslrootcert="/path/to/ca.crt"
)
```

### Audit Logging

All database operations should be logged for security auditing:
```python
import logging
from metainformant.core.utils.logging import get_logger

logger = get_logger("db.audit")

def audited_query(query: str, params: tuple = ()):
    start = time.time()
    try:
        with db.get_connection() as conn_wrapper:
            with conn_wrapper.connect() as conn:
                results = conn_wrapper.execute_query(conn, query, params)
        logger.info(
            "DB_QUERY_SUCCESS",
            extra={
                "query_hash": hashlib.sha256(query.encode()).hexdigest()[:8],
                "duration_s": time.time() - start,
                "rows": len(results),
            }
        )
        return results
    except Exception as e:
        logger.warning(
            "DB_QUERY_FAILED",
            extra={"query": query[:200], "error": str(e)}
        )
        raise
```

## Testing Guidelines

### Test Database Setup

Use a dedicated test database:
```bash
# Create test database
createdb metainformant_test

# Or via SQL:
psql -c "CREATE DATABASE metainformant_test OWNER = your_user;"
```

Configure tests via environment:
```bash
export PG_HOST=localhost
export PG_PORT=5432
export PG_DATABASE=metainformant_test
export PG_USER=test_user
export PG_PASSWORD=test_pass
export TEST_DATABASE_AVAILABLE=1
```

### Test Fixtures

```python
# tests/conftest.py
import pytest
from metainformant.core import db

@pytest.fixture(scope="session")
def db_pool():
    """Create a connection pool for the test session."""
    pool = db.get_db_client()
    yield pool
    pool.close()

@pytest.fixture
def db_conn(db_pool):
    """Provide a connection per test with automatic cleanup."""
    with db_pool.connect() as conn:
        yield conn
        # Optional: clean up data after each test
        conn.rollback()  # Or truncate tables
```

### Real Implementation Policy

Follow METAINFORMANT's real-implementation policy—tests must use real PostgreSQL:
```python
# CORRECT: Real database test
@pytest.mark.network
def test_bulk_insert_real():
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            count = conn_wrapper.bulk_insert(
                conn, "test_table", ["col1", "col2"],
                data=[(1, "a"), (2, "b")]
            )
            assert count == 2

# WRONG: Replacing the database connection layer with a test double hides
# connection, transaction, and serialization behavior that this project tests
# through real adapters.
```

### Temporary Test Tables

```python
def test_query_with_temp_table():
    with db.get_connection() as conn_wrapper:
        with conn_wrapper.connect() as conn:
            # Create temp table (auto-dropped at session end)
            conn_wrapper.execute_query(
                conn,
                "CREATE TEMP TABLE temp_samples (id INT, name TEXT)"
            )
            # Insert test data
            conn_wrapper.execute_query(
                conn,
                "INSERT INTO temp_samples VALUES (1, 'test')"
            )
            # Query it
            results = conn_wrapper.execute_query(
                conn,
                "SELECT * FROM temp_samples WHERE id = 1"
            )
            assert len(results) == 1
```

### Slow Test Marker

Database tests can be slow; use appropriate markers:
```python
import pytest

@pytest.mark.slow
def test_large_bulk_insert():
    # This takes >30 seconds
    pass

@pytest.mark.network
def test_remote_database_connection():
    # Requires network access
    pass
```

## Dependencies

### Required
- `psycopg2` — PostgreSQL adapter (or `psycopg2-binary` for easier install)

### Optional
- None

Install:
```bash
uv add psycopg2
# or:
uv add psycopg2-binary
```

### System Dependencies

For `psycopg2` from source (not binary):
- `libpq-dev` (Debian/Ubuntu) or `postgresql-devel` (RHEL/CentOS)
- Python header files (`python3-dev` or `python3-devel`)

## Repository Conventions

Following METAINFORMANT data flow conventions:
- **Input data**: May be sourced from databases (read operations)
- **Output data**: Query exports go under `output/`
- **Configuration**: Database credentials via environment variables ONLY—never checked into repository
- **Security**: No database passwords in config files or source code

## Troubleshooting Cheat Sheet

| Symptom | Likely Cause | Quick Fix |
|---------|--------------|-----------|
| `PoolError: exhausted` | Too few connections for parallel workers | Increase `max_connections` or decrease `max_workers` |
| `query timeout` or hangs | Missing index on filtered/sorted column | Run `EXPLAIN ANALYZE`; add index |
| `duplicate key violates` | Bulk insert with overlapping data | Use `ON CONFLICT DO NOTHING` or pre-filter data |
| `connection reset by peer` | Connection timeout/idle | Configure `pool_recycle` or reduce `max_connections` |
| `no pg_hba.conf entry` | pg_hba.conf disallows host/user | Update pg_hba.conf to add `host all all <IP> md5` line |
| `permission denied` | DB user missing GRANTs | Ask DBA to grant required privileges |

## Further Reading

- PostgreSQL documentation: https://www.postgresql.org/docs/
- psycopg2 documentation: https://www.psycopg.org/docs/
- SQLAlchemy Core (if using ORM layer): https://docs.sqlalchemy.org/en/20/core/
