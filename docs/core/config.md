# Core: Configuration Management

The `config` module provides flexible configuration loading with support for multiple formats (YAML, TOML, JSON), environment variable overrides, deep merging, and type coercion. It's the central configuration system for all METAINFORMANT domain modules.

## Purpose

Bioinformatics pipelines require configuration across many dimensions:
- **File paths**: Input datasets, reference genomes, output directories
- **Parameters**: Algorithm parameters (e.g., minimum mapping quality, FDR thresholds)
- **Compute resources**: Thread counts, memory limits
- **Environment-specific overrides**: Development vs production credentials

The `config` module unifies these into a simple, unified API.

## Design Principles

### 1. **Multiple Format Support**
Supports YAML (human-friendly, hierarchical), TOML (config-as-code), and JSON (machine-generated). Falls back to a minimal YAML parser if PyYAML is unavailable.

### 2. **Environment Overrides**
Configuration can be overridden via environment variables without editing files. This follows the [Twelve-Factor App](https://12factor.net/config) principle.

### 3. **Deep Merging**
Base configuration files can be overridden by environment-specific configs through recursive dictionary merging, not shallow replacement.

### 4. **Type Coercion**
Environment variables are strings; the module automatically converts them to `int`, `float`, `bool`, or `str` based on a type map.

### 5. **Graceful Degradation**
Optional dependencies (PyYAML, tomllib) are handled gracefully with informative error messages.

### 6. **No Global State**
All functions are pure or accept explicit parameters. No global config singleton.

## Module Organization

The `config` module is in `src/metainformant/core/utils/config.py`. It exposes:

- **Core loaders**: `load_mapping_from_file()`, `load_config_file()`
- **Merging**: `merge_configs()`
- **Environment**: `apply_env_overrides()`, `load_typed_env()`, `load_postgres_config_from_env()`
- **Type coercion**: `coerce_config_types()`
- **Discovery**: `discover_config_files()`, `get_config_schema()`

## Configuration File Formats

### YAML (Recommended for Human-Edited Configs)

```yaml
# config/rna/workflow.yaml
pipeline:
  name: "RNA-seq Analysis"
  threads: 8

input:
  fastq_dir: "/data/fastq"
  sample_sheet: "samples.tsv"

reference:
  genome: "/data/genome/hg38.fa"
  gtf: "/data/annotation/hg38.gtf"

steps:
  - name: "quality_control"
    tool: "fastqc"
    enabled: true
  - name: "quantification"
    tool: "kallisto"
    bootstrap_samples: 100
```

Supported YAML features:
- Mappings (key: value)
- Nested mappings with indentation
- Inline dicts: `{a: 1, b: 2}`
- Lists: `[item1, item2]` or block style with `-`
- Integers, floats, booleans, strings

**No support for** YAML anchors/aliases (`&`, `*`) or complex types (dates, binary). These are rejected by minimal parser and require full PyYAML.

### TOML (Recommended for Config-as-Code)

```toml
# config/pipeline.toml
[pipeline]
name = "GWAS Analysis"
threads = 16

[input.fastq]
directory = "/data/fastq"
pattern = "*.fastq.gz"

[reference.genome]
path = "/data/genome/hg38.fa"
index = "/data/genome/hg38.kidx"

[[steps]]
name = "alignment"
tool = "bwa"
enabled = true
```

Supported via Python 3.11+ `tomllib` (stdlib) or `tomli` backport.

### JSON (Machine-Generated Configs)

```json
{
  "pipeline": {
    "name": "Metagenomics",
    "threads": 4
  },
  "input": {
    "fastq_dir": "/data"
  }
}
```

Supported via `json` module.

## API Reference

### Loading Configuration Files

#### `load_mapping_from_file(config_path: str | Path) -> dict[str, Any]`

Load a configuration mapping from YAML, TOML, or JSON.

**Parameters**:
- `config_path`: Path to configuration file

**Returns**: Dictionary with configuration data

**Raises**:
- `FileNotFoundError` if file doesn't exist
- `ValueError` if format unsupported or parsing fails
- `RuntimeError` if required parser (PyYAML/tomllib) missing

**Example**:
```python
from metainformant.core.utils import config

cfg = config.load_mapping_from_file("config/workflow.yaml")
pipeline_name = cfg["pipeline"]["name"]
```

**Auto-detection**: Based on file suffix (`.yaml`, `.yml`, `.toml`, `.json`).

#### `load_config_file(config_path: Path) -> dict[str, Any]`

Lower-level loader; similar to `load_mapping_from_file()` but returns raw dict (older name; use `load_mapping_from_file()` for new code).

### Environment Variable Overrides

#### `apply_env_overrides(config: Mapping[str, Any], *, prefix: str = "AK") -> dict[str, Any]`

Apply environment variable overrides to a configuration dictionary. Shallow merge only (no deep merging).

**Recognized Variables** (with prefix `AK` by default):

| Variable | Type | Description |
|----------|------|-------------|
| `{prefix}_THREADS` | `int` | Override thread count |
| `{prefix}_WORK_DIR` | `str` | Working/output directory |
| `{prefix}_LOG_DIR` | `str` | Log directory |

**Example**:
```bash
export AK_THREADS=16
export AK_WORK_DIR=/mnt/raid/output
```

```python
base_cfg = {"threads": 4, "work_dir": "output"}
override_cfg = config.apply_env_overrides(base_cfg)
# override_cfg["threads"] == 16 if AK_THREADS set
```

#### `load_typed_env(prefix: str, keys: Mapping[str, type]) -> dict[str, Any]`

Load a typed set of environment variables into a dictionary.

**Parameters**:
- `prefix`: Environment variable prefix (e.g., `"APP"`)
- `keys`: Mapping of variable names (without prefix) to expected types

**Supported types**: `int`, `float`, `bool`, `str`

**Example**:
```python
cfg = config.load_typed_env(
    prefix="MYAPP",
    keys={
        "PORT": int,
        "DEBUG": bool,
        "TIMEOUT": float,
        "HOST": str,
    }
)
# Reads: MYAPP_PORT, MYAPP_DEBUG, MYAPP_TIMEOUT, MYAPP_HOST
# Converts "True" → True, "8080" → 8080, "30.5" → 30.5
```

#### `load_postgres_config_from_env(prefix: str = "PG") -> PostgresConfig | None`

Load PostgreSQL database configuration from environment variables.

**Environment variables** (checked in order):

| Variable | Fallback |
|----------|----------|
| `{prefix}_HOST` | — |
| `{prefix}_PORT` | `"5432"` |
| `{prefix}_DATABASE` | `DB_NAME` |
| `{prefix}_USER` | `DB_USER` |
| `{prefix}_PASSWORD` | `DB_PASSWORD` |

**Returns**: `PostgresConfig` dataclass if all required vars present, `None` otherwise.

**Example**:
```python
pg_cfg = config.load_postgres_config_from_env()
if pg_cfg:
    print(f"Connecting to {pg_cfg.host}:{pg_cfg.port}/{pg_cfg.database}")
```

### Merging and Coercion

#### `merge_configs(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]`

Deep-merge two configuration dictionaries. Nested dictionaries are merged recursively; non-dict values are replaced.

**Algorithm**:
```python
def merge_configs(base, override):
    result = dict(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge_configs(result[key], value)  # Recurse
        else:
            result[key] = value  # Override
    return result
```

**Example**:
```python
base = {
    "threads": 4,
    "database": {"host": "localhost", "port": 5432},
    "output": "results/"
}
override = {
    "threads": 8,
    "database": {"port": 5433},  # host preserved
    "input": {"path": "/data"}    # new key added
}

merged = merge_configs(base, override)
# Result:
# {
#   "threads": 8,
#   "database": {"host": "localhost", "port": 5433},
#   "output": "results/",
#   "input": {"path": "/data"}
# }
```

**Use case**: Layered configuration (defaults → environment → user overrides).

#### `coerce_config_types(config: dict[str, Any], type_map: dict[str, type]) -> dict[str, Any]`

Convert configuration values to expected types in-place.

**Supported types**: `bool`, `int`, `float`, `str`

**Boolean coercion** (case-insensitive):
```text
_bool_map = {"1", "true", "yes", "y", "on"} → True
_bool_map = {"0", "false", "no", "n", "off"} → False
```

**Example**:
```python
config = {"threshold": "0.05", "enabled": "true", "max_iter": "1000"}
type_map = {"threshold": float, "enabled": bool, "max_iter": int}

typed = coerce_config_types(config, type_map)
# typed == {"threshold": 0.05, "enabled": True, "max_iter": 1000}
```

**Raises**: `ValueError` if coercion fails (e.g., `int("abc")`).

### Configuration Discovery

#### `discover_config_files(repo_root: str | Path, domain: str | None = None) -> list[dict[str, Any]]`

Discover all configuration files under `config/` directory.

**Parameters**:
- `repo_root`: Repository root directory
- `domain`: Filter by domain (e.g., `"rna"`, `"gwas"`)

**Returns**: List of config info dicts:
```python
[
    {
        "path": "/abs/path/config/rna/workflow.yaml",
        "domain": "rna",
        "format": "yaml",
        "size": 2048,
        "modified_time": 1703023456.78,
        "relative_path": "rna/workflow.yaml"
    },
    ...
]
```

**Use case**: Discover available pipelines programmatically.

#### `get_config_schema(config_path: str | Path) -> dict[str, Any]`

Infer the structure of a configuration file without fully validating.

**Returns**:
```text
{
    "top_level_keys": ["pipeline", "input", "reference"],
    "nested_structure": {
        "pipeline": {"threads": "int", "name": "str"},
        # ...
    },
    "types": {"pipeline.threads": int, ...}
}
```

**Use case**: Auto-generate documentation or validation schemas.

### PostgreSQL Configuration

#### `PostgresConfig`

Dataclass holding database credentials:

```python
@dataclass(frozen=True)
class PostgresConfig:
    host: str
    port: int
    database: str
    user: str
    password: str
```

#### `get_db_client(**kwargs) -> PostgresConnection`

Convenience function that returns a `PostgresConnection` (from `core.data.db`). See [Database](./db.md) for full details.

## Configuration File Patterns

### Standard Project Layout

```
MetaInformAnt/
├── config/
│   ├── default.yaml          # Base configuration
│   ├── development.yaml      # Dev overrides
│   ├── production.yaml       # Production overrides
│   ├── rna/
│   │   ├── workflow.yaml
│   │   └── samples.tsv
│   ├── gwas/
│   │   └── analysis.yaml
│   └── schemas/
│       └── workflow.json
├── src/metainformant/
└── output/
```

### Domain-Specific Configs

Each domain (`rna`, `gwas`, `dna`, etc.) may have its own config schema. The `discover_config_files()` function indexes them.

**Typical domain config sections**:
- `pipeline`: Global pipeline settings
- `input`: Data sources, sample sheets, file patterns
- `reference`: Genome/annotation files, indices
- `parameters`: Algorithm-specific parameters
- `output`: Output directory, naming patterns
- `resources`: Compute resources (threads, memory)

### Config Loading Pattern

```python
from metainformant.core.io import paths
from metainformant.core.utils import config

def load_pipeline_config(domain: str, env: str = "production") -> dict:
    """Load layered configuration for a domain."""
    repo_root = paths.get_project_root()
    config_dir = repo_root / "config"

    # 1. Base config
    base_path = config_dir / "default.yaml"
    base_cfg = config.load_mapping_from_file(base_path)

    # 2. Environment override
    env_path = config_dir / f"{env}.yaml"
    env_cfg = config.load_mapping_from_file(env_path) if env_path.exists() else {}

    # 3. Domain-specific config
    domain_path = config_dir / domain / "workflow.yaml"
    domain_cfg = config.load_mapping_from_file(domain_path) if domain_path.exists() else {}

    # 4. Merge: base ← env ← domain
    merged = config.merge_configs(base_cfg, env_cfg)
    merged = config.merge_configs(merged, domain_cfg)

    # 5. Apply env vars
    merged = config.apply_env_overrides(merged)

    return merged
```

## Environment Variables: Complete Reference

### Standard METAINFORMANT Overrides

All prefixed with `AK_` (AmalgKit legacy), configurable via `prefix` parameter:

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `AK_THREADS` | `int` | CPU count | Number of worker threads |
| `AK_WORK_DIR` | `str` | `"output"` | Output directory for results |
| `AK_LOG_DIR` | `str` | `"logs"` | Directory for log files |

**Example**:
```bash
export AK_THREADS=32
export AK_WORK_DIR=/mnt/raid/results
export AK_LOG_DIR=/var/log/metainformant
```

### PostgreSQL Variables

Either `PG_*` or `DB_*` forms:

| Variable | Meaning |
|----------|---------|
| `PG_HOST` / `DB_HOST` | Database server hostname |
| `PG_PORT` / `DB_PORT` | Port (default: 5432) |
| `PG_DATABASE` / `DB_NAME` | Database name |
| `PG_USER` / `DB_USER` | Username |
| `PG_PASSWORD` / `DB_PASSWORD` | Password |

Priority: `PG_*` checked first, then `DB_*` as fallback for some keys.

### Logging

| Variable | Type | Default |
|----------|------|---------|
| `CORE_LOG_LEVEL` | `str` (DEBUG/INFO/WARNING/ERROR) | `INFO` |

## Common Pitfalls and Troubleshooting

### Pitfall 1: YAML Parsing Fails Without PyYAML

**Symptom**: `RuntimeError: PyYAML not available for YAML config files`

**Cause**: You're using `.yaml` file but `PyYAML` package not installed.

**Fix**: Install PyYAML or use `.json` config:
```bash
uv pip install pyyaml
```

Alternatively, use the minimal YAML subset (simple key: value mappings only). The built-in `_parse_simple_yaml_mapping()` handles:
- Key-value pairs (strings, integers)
- Empty dicts `{}` and lists `[]`
- Inline dicts `{a: 1}`
- One-level nested mappings and block lists

Complex YAML features (anchors `&`, aliases `*`, multi-line strings `|`, `>`) require full PyYAML.

### Pitfall 2: Deep Merge Not Happening

**Symptom**: Override config's nested dict completely replaces base config's nested dict instead of merging.

**Cause**: Using simple `dict.update()` instead of `merge_configs()`.

**Fix**:
```python
# WRONG
config = {**base, **override}  # Shallow merge; nested dict overwritten

# CORRECT
config = config.merge_configs(base, override)  # Deep merge
```

### Pitfall 3: Environment Variable Type Errors

**Symptom**: `ValueError: invalid literal for int() with base 10: 'abc'`

**Cause**: Non-numeric string in numeric env var (`AK_THREADS=abc`).

**Fix**: Validate or use defaults. `apply_env_overrides()` catches `ValueError` and logs warning, falling back to default.

### Pitfall 4: Missing Config File

**Symptom**: `FileNotFoundError: Config file not found: ...`

**Fix**:
- Check path relative to current working directory
- Use absolute path or `Path(__file__).parent / "config.yaml"`
- Ensure file committed to repository

### Pitfall 5: Sensitive Data in Config Files

**Symptom**: Database passwords in plain text YAML checked into Git.

**Fix**: Use environment variables for secrets:
```yaml
# config/database.yaml
database:
  host: ${PG_HOST}
  user: ${PG_USER}
  password: ${PG_PASSWORD}  # Not in file; injected at runtime
```

The module doesn't interpolate `${VAR}` automatically; use `os.getenv()` in code:
```python
db_cfg = config.load_mapping_from_file("config/db.yaml")
db_cfg["password"] = os.getenv("DB_PASSWORD", db_cfg.get("password"))
```

Or use separate `config/production-secrets.yaml` (git-ignored) to load optionally.

### Pitfall 6: TOML in Python < 3.11

**Symptom**: `ImportError: tomllib not available`

**Cause**: Using `.toml` files on Python 3.10 or earlier.

**Fix**:
- Upgrade to Python 3.11+
- Or install `tomli` backport: `uv pip install tomli`
- The module tries `importlib.import_module("tomllib")` and falls back gracefully

## Performance Considerations

### Config Loading Overhead

- **YAML + PyYAML**: ~5–20 ms for 10 KB file (parser overhead)
- **TOML + tomllib**: ~2–10 ms (stdlib, faster)
- **JSON**: ~1–5 ms (fastest)

For large configs (MBs), consider splitting or using binary formats (MessagePack).

### Caching Configs

Configs rarely change at runtime. Load once at startup:

```python
# Module-level singleton (common pattern)
CONFIG = config.load_mapping_from_file("config/pipeline.yaml")
```

No built-in caching; reload on SIGHUP if needed.

### Merging Large Configs

`merge_configs()` is recursive and O(n) in total keys. For very large configs (10k+ keys), profiling shows negligible overhead (<1 ms).

### Environment Variable Lookup

`os.getenv()` is fast. No need to cache env var values unless in tight loops (unlikely for config).

## Security Considerations

### Path Traversal in Configs

Config files may contain paths. Always sanitize:

```python
user_path = cfg["input"]["path"]
if not paths.is_safe_path(user_path):
    raise ValueError("Unsafe path in config")
abs_path = paths.expand_and_resolve(user_path)
```

### SQL Injection Risk

Configs may contain SQL queries. Never interpolate directly:

```python
# BAD
query = f"SELECT * FROM {cfg['table']}"  # Injection risk

# GOOD
query = "SELECT * FROM table WHERE column = %s"
cursor.execute(query, (cfg["value"],))
```

### Secrets in Config Files

Never commit passwords/tokens. Use environment variables or secret managers (AWS Secrets Manager, HashiCorp Vault).

## Type-checking and Validation

The config module does **not** validate schema or types automatically. You must validate:

```python
from metainformant.core.data.validation import validate_schema

schema = {
    "pipeline": {"name": str, "threads": int},
    "input": {"fastq_dir": str}
}
errors = validate_schema(cfg, schema)
if errors:
    raise ValueError(f"Invalid config: {errors}")
```

Or use `pydantic` (optional):
```python
from pydantic import BaseModel, Field

class PipelineConfig(BaseModel):
    name: str
    threads: int = Field(ge=1, le=64)

validated = PipelineConfig(**cfg["pipeline"])
```

## Integration with Other Core Modules

```python
from metainformant.core import io
from metainformant.core.io import paths
from metainformant.core.utils import config
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Load config
cfg = config.load_mapping_from_file("config/workflow.yaml")

# Apply environment overrides
cfg = config.apply_env_overrides(cfg)

# Resolve paths relative to project root
project_root = paths.get_project_root()
input_dir = paths.expand_and_resolve(cfg["input"]["dir"])
input_dir = paths.is_within(input_dir, project_root) or project_root

# Use typed values
threads = int(cfg.get("pipeline", {}).get("threads", 4))

logger.info(f"Running with {threads} threads on {input_dir}")
```

## Dependencies

- **Required**: Python stdlib `os`, `pathlib`, `json`
- **Optional**: `PyYAML` for full YAML support, `tomli`/`tomllib` for TOML (stdlib in Python 3.11+)

## Testing

Config functions are pure; no mocking needed:
```python
def test_merge_configs():
    base = {"a": 1, "b": {"x": 10}}
    override = {"b": {"y": 20}, "c": 3}
    merged = config.merge_configs(base, override)
    assert merged == {"a": 1, "b": {"x": 10, "y": 20}, "c": 3}
```

## Summary

The `config` module is the configuration backbone of METAINFORMANT. Its simplicity (load/merge/override/coerce) encourages declarative pipeline configuration while maintaining flexibility for environment-specific customization.
