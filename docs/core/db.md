# Core: db

Function: `get_db_client`

Requires environment variables for Postgres (`PG_HOST`, `PG_DATABASE`, `PG_USER`, `PG_PASSWORD`, `PG_PORT`).

```python
from metainformant.core import db

conn, cur = db.get_db_client()
cur.execute("SELECT 1")
```

Notes

- Credentials are sourced from environment variables (see `core.config`). Do not store secrets under `config/`.
- Query results and diagnostics written to disk should go under `output/`.
