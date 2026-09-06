# Data

Input validation utilities and PostgreSQL database connectivity for the core data layer.

## Contents

| File | Purpose |
|------|---------|
| `validation.py` | Type, range, path, and schema validators with decorator support |
| `db.py` | PostgreSQL connection management, query execution, and URL building |

## Key Functions and Classes

| Symbol | Description |
|--------|-------------|
| `validate_type()` | Assert a value matches an expected type |
| `validate_range()` | Assert a numeric value falls within min/max bounds |
| `validate_path_exists()` | Confirm a filesystem path exists and return a resolved Path |
| `validate_path_is_file()` | Confirm a path exists and is a file |
| `validate_path_is_dir()` | Confirm a path exists and is a directory |
| `validate_path_within()` | Ensure a path stays within a parent directory (security) |
| `validate_not_empty()` | Reject empty strings, lists, or dicts |
| `validate_not_none()` | Reject None values |
| `validate_schema()` | Check a dict against a required-keys schema |
| `validate_json_schema()` | Validate data against a JSON Schema file |
| `validator()` | Decorator that converts a bool predicate into a raising validator |
| `PostgresConnection` | Pooled PostgreSQL connection wrapper with query helpers |
| `build_postgres_url()` | Construct a connection URL from host, port, database, and credentials |
| `get_db_client()` | Factory that returns a configured PostgresConnection |
| `get_connection()` | Context manager that closes the pool on exit |
| `sanitize_connection_params()` | Strip unsafe characters from connection parameters |

## Usage

```python
from metainformant.core.data.validation import validate_type, validate_path_within
from metainformant.core.data.db import PostgresConnection, build_postgres_url, get_connection

validate_type(value, int, name="count")
validate_path_within("/safe/base", user_path)

url = build_postgres_url(host="localhost", port=5432, database="bio", user="postgres", password="secret")
client = PostgresConnection(host="localhost", database="bio", user="postgres", password="secret")
with get_connection(
    host="localhost", port=5432, database="bio", user="postgres", password="secret"
) as client:
    with client.connect() as conn:
        rows = client.execute_query(conn, "SELECT * FROM samples")
```
