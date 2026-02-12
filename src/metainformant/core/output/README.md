# Output

Designated output directory for core module results. All execution outputs from the core module are written here following the project convention of routing results to `output/`.

## Contents

This directory serves as the runtime output location for core operations. It may contain cached discovery data and intermediate results. No source modules are defined here.

## Convention

All metainformant modules write their outputs to `output/<module>/` directories. The core output directory is managed automatically by the IO and paths utilities:

```python
from metainformant.core.io.paths import get_project_root

output_dir = get_project_root() / "output" / "core"
output_dir.mkdir(parents=True, exist_ok=True)
```
