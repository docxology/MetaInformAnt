# Core: Configuration Management

The `config` module provides comprehensive configuration management with support for multiple file formats, environment variable overrides, and type-safe loading.

## Classes

### `PostgresConfig`
Dataclass for PostgreSQL database configuration.

```python
@dataclass(frozen=True)
class PostgresConfig:
    host: str
    port: int
    database: str
    user: str
    password: str
```

## Functions

### Database Configuration
- **`load_postgres_config_from_env(prefix="PG")`** → `PostgresConfig | None`
  - Load PostgreSQL config from environment variables
  - Supports custom prefix for variable names
  - Returns None if required variables are missing

### File-Based Configuration
- **`load_config_file(config_path)`** → `Dict[str, Any]`
  - Load configuration from YAML, TOML, or JSON files
  - Supports .yaml/.yml, .toml, and .json formats
  - Raises ValueError for unsupported formats

- **`load_mapping_from_file(config_path)`** → `Dict[str, Any]`
  - Load structured mapping from config files
  - Same format support as load_config_file
  - Includes fallback YAML parser for simple configs

### Environment Variable Handling
- **`get_env_or_default(env_var, default)`** → `str`
  - Get environment variable with fallback default

- **`load_typed_env(prefix, keys)`** → `Dict[str, Any]`
  - Load typed environment variables with prefix
  - Supports int, float, bool, and str types
  - Automatic type coercion and validation

- **`apply_env_overrides(config, prefix="AK")`** → `Dict[str, Any]`
  - Apply environment overrides to config mapping
  - Supports threads, work_dir, and log_dir overrides

## Usage Examples

### PostgreSQL Configuration
```python
from metainformant.core import config

# Load from environment variables
pg_config = config.load_postgres_config_from_env()
if pg_config:
    print(f"Connected to {pg_config.database} at {pg_config.host}:{pg_config.port}")
```

### Typed Environment Loading
```python
from metainformant.core import config

# Load typed environment variables
app_config = config.load_typed_env(
    prefix="APP",
    keys={"PORT": int, "DEBUG": bool, "TIMEOUT": float}
)
print(f"Port: {app_config.get('PORT', 8080)}")
```

### File-Based Configuration
```python
from metainformant.core import config

# Load from YAML/TOML/JSON
try:
    cfg = config.load_mapping_from_file("config.yaml")
    cfg_with_overrides = config.apply_env_overrides(cfg)
except FileNotFoundError:
    print("Config file not found")
```

### Environment Overrides
```python
from metainformant.core import config

# Apply environment overrides to config
base_config = {"threads": 4, "work_dir": "/tmp/work"}
final_config = config.apply_env_overrides(base_config)
# Can be overridden with AK_THREADS, AK_WORK_DIR env vars
```

## Supported File Formats

### YAML (.yaml, .yml)
```yaml
database:
  host: localhost
  port: 5432
  name: mydb

compute:
  threads: 8
  memory_gb: 16
```

### TOML (.toml)
```toml
[database]
host = "localhost"
port = 5432
name = "mydb"

[compute]
threads = 8
memory_gb = 16
```

### JSON (.json)
```json
{
  "database": {
    "host": "localhost",
    "port": 5432,
    "name": "mydb"
  },
  "compute": {
    "threads": 8,
    "memory_gb": 16
  }
}
```

## Environment Variables

### Database Configuration
- `PG_HOST` - Database host
- `DB_NAME` or `PG_DATABASE` - Database name
- `DB_USER` or `PG_USER` - Database user
- `DB_PASSWORD` or `PG_PASSWORD` - Database password
- `PG_PORT` - Database port (default: 5432)

### Application Overrides
- `AK_THREADS` - Number of threads
- `AK_WORK_DIR` - Working directory path
- `AK_LOG_DIR` - Logging directory path

## Error Handling

The config module provides robust error handling:
- Clear error messages for missing files
- Graceful handling of optional dependencies (PyYAML, tomllib)
- Type validation for environment variables
- Fallback parsers for simple configurations

## Dependencies

- **Required**: Standard library only
- **Optional**: PyYAML (for YAML), tomllib (Python 3.11+ for TOML)
- **Database**: psycopg2 (optional, for PostgreSQL integration)
