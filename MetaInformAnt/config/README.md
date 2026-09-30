# config

Configuration files for METAINFORMANT domain modules and the amalgkit
campaign. See each subfolder's README.

## Environment Variable Overrides

Configuration supports environment-variable overrides via
`metainformant.core.utils.config.apply_env_overrides` (module-specific
prefixes, e.g. `AMALGKIT_` for RNA/amalgkit settings; the NCBI client reads
`NCBI_EMAIL`). Example:

```python
from metainformant.core.utils.config import load_mapping_from_file, apply_env_overrides

config = load_mapping_from_file("config/ncbi/ncbi.yaml")
config = apply_env_overrides(config, prefix="NCBI")
```

## Quick Start

```bash
# Inspect the canonical project inventory (never omit --data-root)
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
uv run python scripts/rna/run_all_species.py   --config-dir projects/hymenoptera_amalgkit/config/amalgkit   --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Related Documentation

- [RNA Workflow Guide](../docs/rna/workflow.md)
- [GWAS Workflow Guide](../docs/gwas/workflow.md)
- [Core Config Module](../src/metainformant/core/utils/config.py)
- Amalgkit campaign configs live under
  `projects/hymenoptera_amalgkit/config/amalgkit/` (parent-repo mirror;
  28 species YAMLs plus template, test, and cross-species configs).
