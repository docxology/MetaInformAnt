# Architecture Guide — Thin Orchestration Pattern

> [!NOTE]
> This guide explains the core design principle of all MetaInformAnt projects and how the template is structured to enforce it.

## The Thin Orchestration Pattern

In a **thin orchestrator**, each numbered script in `scripts/` does exactly four things:

1. **Parse** CLI arguments (`--config`, `--force`, etc.)
2. **Load** configuration from the YAML file
3. **Delegate** all computation to `metainformant.*` library functions
4. **Record** provenance (timing, input/output file metadata) to `logs/`

Scripts contain **no domain algorithms**. There are no hand-rolled statistical functions, no inline data-wrangling loops that belong in a library, and no hardcoded paths. If business logic is needed, it goes into the appropriate `metainformant.*` module and is imported.

### Why?

| Concern | Without Thin Orchestration | With Thin Orchestration |
| :--- | :--- | :--- |
| Reproducibility | Logic scattered in scripts, version-ambiguous | Versioned library, tested independently |
| Reuse | Copy-paste between projects | Import and call |
| Testing | Hard to unit-test embedded logic | Library functions tested; scripts tested via subprocess |
| Maintenance | Fix the same bug in N projects | Fix once in the library |

## MetaInformAnt Module Mapping

| Domain | MetaInformAnt Module | Use in this Template |
| :--- | :--- | :--- |
| I/O utilities | `metainformant.core.io` | `load_json`, `dump_json`, `read_csv`, `write_csv` |
| Path safety | `metainformant.core.paths` | `expand_and_resolve`, `ensure_directory` |
| Validation | `metainformant.core.validation` | `validate_path_exists`, `validate_type` |
| Caching | `metainformant.core.cache` | `JsonCache` for idempotency manifests |
| Logging | `metainformant.core.utils.logging` | `get_logger`, `setup_logging` |
| DNA sequences | `metainformant.dna.sequences` | When sequencing data is involved |
| RNA analysis | `metainformant.rna.amalgkit` | When transcriptomic data is involved |
| GWAS | `metainformant.gwas` | When genome-association work is involved |

> [!TIP]
> When adapting this template for a specific domain, replace `pandas` / `numpy` direct calls with the appropriate `metainformant.*` module functions for maximum reuse and testing coverage.

## Adding External Tool Dependencies

If your pipeline needs CLI tools (e.g., `bwa`, `samtools`, `bcftools`), follow this pattern:

```python
import subprocess
import shutil

def check_tool(name: str) -> None:
    if not shutil.which(name):
        raise RuntimeError(f"Required tool not found in PATH: {name}")

check_tool("samtools")
result = subprocess.run(["samtools", "sort", ...], check=True)
```

Document the tool in `doc/index.md` under **Prerequisites** and in `README.md`.

## Idempotency Pattern

Every script checks for its primary output before doing any work:

```python
if output_file.exists() and not force:
    logger.info("Output already exists, skipping (use --force to reprocess).")
    return
```

This makes re-runs safe and cheap. The `--force` flag bypasses the guard for full reprocessing.

## Extending the Pipeline

### Adding a New Stage

1. Create `scripts/NN_stage_name.py` following the docstring + argparse + idempotency template.
2. Delegate logic to `metainformant.*`.
3. Write logs to `logs/NN_stage_name.log`.
4. Register in `run.sh` (add function + case block) and `main.py` (`STAGES` dict + `ORDER` list).
5. Add doc at `doc/stages/NN_stage_name.md`.
6. Add tests in `tests/test_pipeline.py`.

### Forking for a New Project

1. Copy this template directory.
2. Rename in `pyproject.toml` (`name`, `description`).
3. Adapt `config/default.yaml` for your organism, paths, and parameters.
4. Replace or augment scripts with domain-specific stages.
5. Update `README.md`, `AGENTS.md`, `SPEC.md`, and `doc/`.
