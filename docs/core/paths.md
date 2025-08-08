# Core: paths

Utilities for safe, predictable path handling across the repository.

Functions: `expand_and_resolve`, `is_within`

Conventions

- Root-level directories:
  - `config/`: configuration files and run options
  - `data/`: datasets and local databases (inputs)
  - `output/`: all outputs from tests and real runs (ephemeral)

Usage example

```python
from pathlib import Path
from metainformant.core.paths import expand_and_resolve, is_within

root = expand_and_resolve(".")
out_dir = expand_and_resolve("output/run_001")
assert is_within(out_dir, root)
Path(out_dir).mkdir(parents=True, exist_ok=True)
```

