# Core: Database Integration

The `db` module provides PostgreSQL database connectivity with secure credential management and connection handling.

## Functions

### Connection Management
- **`get_db_client()`** → `Tuple[connection, cursor]`
  - Establish connection to PostgreSQL database
  - Load configuration from environment variables
  - Returns connection and cursor objects
  - Performs basic connectivity test (SELECT 1)

### URL Construction
- **`build_postgres_url(host, port, database, user, password)`** → `str`
  - Build properly formatted PostgreSQL connection URL
  - Useful for SQLAlchemy and other database libraries
  - Handles URL encoding automatically

### Security
- **`sanitize_connection_params(params)`** → `Dict`
  - Sanitize database connection parameters
  - Remove SQL injection attempts and dangerous characters
  - Filter out dangerous SQL keywords

## Usage Examples

### Basic Connection
```python
from metainformant.core import db

try:
    conn, cur = db.get_db_client()

    # Execute queries
    cur.execute("SELECT COUNT(*) FROM samples")
    count = cur.fetchone()[0]
    print(f"Found {count} samples")

    # Always close connections
    cur.close()
    conn.close()

except RuntimeError as e:
    print(f"Database connection failed: {e}")
```

### Environment Configuration
Set these environment variables before using database functions:

```bash
export PG_HOST="localhost"
export PG_PORT="5432"
export PG_DATABASE="metainformant"
export PG_USER="analysis_user"
export PG_PASSWORD="secure_password"
```

Or using the shorter aliases:
```bash
export DB_NAME="metainformant"
export DB_USER="analysis_user"
export DB_PASSWORD="secure_password"
```

### Connection URL Building
```python
from metainformant.core import db

# Build connection URL for SQLAlchemy or other tools
url = db.build_postgres_url(
    host="db.example.com",
    port=5432,
    database="biology_db",
    user="researcher",
    password="secure123"
)
print(f"Connection URL: {url}")
# postgresql://researcher:secure123@db.example.com:5432/biology_db
```

### Secure Parameter Handling
```python
from metainformant.core import db

# Sanitize user input before using in connection
raw_params = {
    "host": "localhost; DROP TABLE users;--",
    "database": "analysis_db",
    "user": "researcher",
    "password": "pass123"
}

safe_params = db.sanitize_connection_params(raw_params)
# Dangerous characters and SQL keywords removed
```

## Integration with Core Modules

```python
from metainformant.core import db, io, config

def export_query_results(query, output_file):
    """Execute query and export results to JSON."""
    conn, cur = db.get_db_client()

    try:
        cur.execute(query)
        columns = [desc[0] for desc in cur.description]
        rows = cur.fetchall()

        # Convert to list of dicts
        results = [dict(zip(columns, row)) for row in rows]

        # Export to output directory
        io.dump_json(results, f"output/{output_file}")

        return len(results)

    finally:
        cur.close()
        conn.close()

# Usage
count = export_query_results(
    "SELECT * FROM gene_expressions WHERE species = 'human'",
    "human_expressions.json"
)
print(f"Exported {count} records")
```

## Error Handling

Database operations include robust error handling:
- Connection failures raise clear RuntimeError messages
- Invalid credentials or unreachable hosts handled gracefully
- Parameter sanitization prevents injection attacks
- Automatic cleanup with context managers recommended

## Best Practices

### Connection Management
```python
# Recommended: Use context managers for automatic cleanup
from contextlib import contextmanager

@contextmanager
def get_db_connection():
    conn, cur = db.get_db_client()
    try:
        yield conn, cur
    finally:
        cur.close()
        conn.close()

# Usage
with get_db_connection() as (conn, cur):
    cur.execute("SELECT * FROM experiments")
    results = cur.fetchall()
```

### Security Considerations
- **Never store credentials in code or config files**
- **Use environment variables for all database credentials**
- **Apply parameter sanitization to user inputs**
- **Limit database user permissions to required operations**
- **Use SSL/TLS connections in production**

### Performance Optimization
- **Reuse connections when possible** (avoid opening/closing frequently)
- **Use appropriate indexing** on frequently queried columns
- **Batch operations** for bulk inserts/updates
- **Connection pooling** for high-throughput applications

## Dependencies

- **Required**: `psycopg2` (PostgreSQL adapter for Python)
- **Optional**: None

Install with:
```bash
uv add psycopg2
```

## Repository Conventions

Following the METAINFORMANT data flow conventions:
- **Input data**: May be sourced from databases (read operations)
- **Output data**: Query results and exports go under `output/`
- **Configuration**: Database credentials via environment variables only
- **Security**: No secrets stored in `config/` directory
