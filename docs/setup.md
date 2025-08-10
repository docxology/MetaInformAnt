# Setup

Environment management uses `uv` with Python 3.11+.

```bash
git clone <repo>
cd METAINFORMANT
bash scripts/setup_uv.sh --with-amalgkit --ncbi-email "you@example.com"
# You can either use uv run directly, or activate the venv for the shell:
# uv run python -V
source .venv/bin/activate
```

Verify:

```bash
uv run python -V
uv run pytest -q
uv run metainformant --help
```

Directories policy


- `config/`: repository-level configuration and options; read via `metainformant.core.config` with environment overrides.
- `data/`: datasets and local databases (inputs), organized by domain/version.
- `output/`: all outputs from tests and real runs. Ephemeral and reproducible; safe to delete.

Examples respect this policy by defaulting to `output/` when writing files.

External tools

- RNA: `amalgkit` (installed automatically by setup when `--with-amalgkit` is provided, via uv)
- Optional MSA: `muscle` or `clustalo` in PATH
- Optional NCBI Datasets: `ncbi-datasets-pylib` is a dependency; verify availability in runtime

Next steps
- Explore [CLI](./cli.md)
- DNA quickstart: [DNA overview](./dna/index.md)
- RNA workflow: [RNA overview](./rna/index.md)


