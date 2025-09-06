### Core: config

Classes: `PostgresConfig`

Functions: `load_postgres_config_from_env`, `load_typed_env`

```python
from metainformant.core import config

pg = config.load_postgres_config_from_env()
typed = config.load_typed_env(prefix="APP", keys={"PORT": int, "DEBUG": bool})
```
