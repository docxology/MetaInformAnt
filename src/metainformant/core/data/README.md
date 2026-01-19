# Core Data Module

Data structures, validation, and database integration for MetaInformAnt.

## Purpose

This module provides:
- Data validation utilities
- PostgreSQL database integration helpers
- Type-safe data containers

## Key Components

| File | Description |
|------|-------------|
| [validation.py](validation.py) | Type, range, path, and schema validators |
| [db.py](db.py) | `PostgresConnection` class for database operations |

## Usage

```python
from metainformant.core.data import validate_type, validate_range

validate_type(value, int, "count")
validate_range(value, min_val=0, max_val=100, name="percent")
```

## Related Documentation

- **Parent**: [src/metainformant/core/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **I/O Module**: [../io/README.md](../io/README.md)
