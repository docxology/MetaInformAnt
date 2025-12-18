# Setup

Environment management uses `uv` with Python 3.11+.

## Quick Setup

```bash
git clone https://github.com/q/metainformant.git
cd metainformant
bash scripts/package/setup.sh --ncbi-email "you@example.com"
# Default setup installs dev + scientific dependencies + amalgkit
# You can either use uv run directly, or activate the venv for the shell:
# uv run python -V
source .venv/bin/activate  # or /tmp/metainformant_venv/bin/activate on FAT filesystems
```

**Note**: Setup scripts automatically detect FAT filesystems (exFAT, FAT32) and configure UV cache and virtual environment locations accordingly. See [UV Setup Guide](UV_SETUP.md) for details.

## FAT Filesystem Support

On FAT-formatted external drives, setup is **automatic**:

```bash
# Same command works on FAT filesystems
bash scripts/package/setup.sh

# The script automatically:
# - Detects FAT filesystem
# - Sets UV_CACHE_DIR=/tmp/uv-cache
# - Creates venv at /tmp/metainformant_venv (if needed)
```

See [UV Setup Guide](UV_SETUP.md) for complete FAT filesystem documentation.

## Verify Setup

```bash
# Standard verification
uv run python -V
uv run pytest -q
uv run metainformant --help

# Comprehensive verification (includes FAT filesystem checks)
bash scripts/package/verify_uv_setup.sh
```

Directories policy


- `config/`: repository-level configuration and options; read via `metainformant.core.config` with environment overrides.
- `data/`: datasets and local databases (inputs), organized by domain/version.
- `output/`: all outputs from tests and real runs. Ephemeral and reproducible; safe to delete.

Examples respect this policy by defaulting to `output/` when writing files.

External tools

- RNA: `amalgkit` (installed automatically by default setup script, via uv)
- Optional MSA: `muscle` or `clustalo` in PATH (MUSCLE shim installed automatically if not found)
- Optional NCBI Datasets: `ncbi-datasets-pylib` is a dependency; verify availability in runtime

Next steps
- Explore [CLI](./cli.md)
- DNA quickstart: [DNA overview](./dna/index.md)
- RNA workflow: [RNA overview](./rna/index.md)
