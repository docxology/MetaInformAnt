### RNA: amalgkit Wrapper

Utilities converting parameter dicts to CLI flags and executing `amalgkit` subcommands.

Key functions: `build_cli_args`, `build_amalgkit_command`, `check_cli_available`, `run_amalgkit`, and convenience wrappers per subcommand: `metadata`, `integrate`, `config`, `select`, `getfastq`, `quant`, `merge`, `cstmm`, `curate`, `csca`, `sanity`.

```mermaid
flowchart LR
  A[params dict] --> B[build_cli_args]
  B --> C[build_amalgkit_command]
  C --> D[run_amalgkit]
  D --> E[return code/stdout/stderr]
```

Example

```python
from metainformant.rna import amalgkit as ak

ok, help_text = ak.check_cli_available()
cmd = ak.build_amalgkit_command("metadata", {"species-list": ["Apis_mellifera"], "threads": 8})
```


