# Development Guide for METAINFORMANT

This guide covers development environment setup, workflow, testing, debugging, and adding new modules. It's intended for contributors who plan to modify or extend METAINFORMANT.

For user-oriented quick start, see [GETTING_STARTED.md](GETTING_STARTED.md). For contribution guidelines, see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Development Environment

### Prerequisites
- **Python 3.11 or 3.12** (managed by `uv`)
- **Git**
- **uv** – fast Python package manager and resolver ([install](https://docs.astral.sh/uv/))

### Virtual Environment Setup with `uv`

`uv` replaces `virtualenv` + `pip`. All development should use a project-local virtual environment.

#### Create and Activate a Virtual Environment

```bash
# From repository root
uv venv  # creates .venv/ by default

# Activate (bash/zsh)
source .venv/bin/activate

# Or use uv run without activating (recommended)
uv run python -c "print('Hello')"
```

On FAT filesystems (external drives), the setup scripts automatically configure a temp-based venv. See [UV_SETUP.md](UV_SETUP.md#fat-filesystem-support).

#### Install Dependencies

```bash
# Editable install with all dev dependencies
uv pip install -e ".[dev,scientific,ml,networks,singlecell]"

# For full installation with all optional extras:
uv pip install -e ".[all]"

# Or use the setup script (installs amalgkit, verifies external tools):
bash scripts/package/setup.sh
```

The `dev` group includes: pytest,black,isort,mypy,pre-commit,sphinx.

#### Verify Installation

```bash
uv run python -V                   # Should be >=3.11
uv run pytest -q                   # Quick test sanity check
uv run metainformant --help         # CLI entry point
```

---

### UV Workflow in Daily Development

| Task | Command |
|------|---------|
| Run a script | `uv run python path/to/script.py` |
| Run tests | `uv run pytest` or `uv run scripts/run_tests.sh` |
| Format code | `uv run black .` |
| Sort imports | `uv run isort .` |
| Type check | `uv run mypy src/metainformant` |
| Launch REPL | `uv run python` (has access to package) |
| Lint | `uv run flake8 src/metainformant` |

Prefer `uv run` over activating the venv manually; it's consistent and avoids activation state confusion.

---

### Pre-commit Hooks

Install once after cloning:
```bash
uv run pre-commit install
```

Hooks run on every `git commit`. They check:
- Formatting (black, isort)
- Import order
- Large files
- Merge conflicts
- YAML/JSON/TOML syntax
- End-of-line whitespace
- Private keys
- Example code validity

Force-run on all files (e.g., after changing config):
```bash
uv run pre-commit run --all-files
```

Some hooks are disabled by default (mypy, flake8, bandit, pydocstyle). You can enable them locally in `.pre-commit-config.yaml` if desired.

---

## Testing

### Running Tests

METAINFORMANT uses `pytest` with a strict configuration.

**Quick run (fast tests only):**
```bash
uv run pytest -q
```

**Full suite (including slow/network tests):**
```bash
uv run pytest -v
```

**Parallel execution (much faster):**
```bash
uv run pytest -n auto     # one worker per CPU core
```

**With coverage:**
```bash
uv run pytest --cov=src/metainformant --cov-report=term-missing --cov-report=html
# HTML report: output/coverage_html/index.py
```

**Exclude slow tests:**
```bash
uv run pytest -m "not slow"
```

**Run a specific test file:**
```bash
uv run pytest tests/dna/test_composition.py -v
```

**Run the CI wrapper script** (mirrors CI environment):
```bash
bash scripts/run_tests.sh
```

### Test Organization

- `tests/` mirrors `src/metainformant/`
- `tests/conftest.py` – shared fixtures (temporary directories, sample data)
- Markers:
  - `@pytest.mark.slow` – long-running tests
  - `@pytest.mark.network` – requires internet (real API calls)
  - `@pytest.mark.external_tool` – requires external binaries (e.g., muscle, amalgkit)
  - `@pytest.mark.integration` – full pipeline integration tests
  - `@pytest.mark.no_mock` – explicitly marks tests adhering to No Mocking policy

### The No-Mocking Policy

**Zero tolerance for mocks, fakes, or placeholders.** All functions must perform real computations or genuine I/O. Tests must use real implementations:

- **Good**: actual file I/O, real API calls to UniProt/NCBI, genuine algorithms
- **Bad**: returning hardcoded data, skipping computations with `return None`, using `unittest.mock.patch` to replace real functions with stubs

See [docs/NO_MOCKING_POLICY.md](NO_MOCKING_POLICY.md) for rationale, examples, and enforcement details.

### Writing New Tests

- Place tests in `tests/<module>/test_<feature>.py`
- Follow **one assertion per test** principle (one behavior per test)
- Use fixtures from `conftest.py` for reusable resources
- Keep tests fast; mark slow ones with `@pytest.mark.slow`
- Real data only: generate synthetic data with `numpy.random` or use small example files stored under `tests/data/`

Example test:
```python
def test_gc_content_basic() -> None:
    from metainformant.dna.composition import gc_content
    assert gc_content("ATCG") == 0.5
    assert gc_content("AAAA") == 0.0
    assert gc_content("") == 0.0
```

### Debugging Failing Tests

Use `pdb` or `breakpoint()`:
```bash
uv run pytest tests/dna/test_composition.py::test_gc_content_basic --pdb
```

Or add `breakpoint()` in test or source code; `uv run pytest` will drop into debugger on failure.

For print debugging: `uv run pytest -s` to see stdout/stderr.

---

## Adding New Modules

A typical module lives in `src/metainformant/<module>/` and has the following structure:

```
src/metainformant/<module>/
├── __init__.py          # Public API exports
├── README.md           # User-facing documentation (overview, examples)
├── AGENTS.md           # AI agent guidance (this file pattern)
├── SPEC.md             # Technical specification (API signatures, data structures)
├── <submodule1>/       # Functional subpackages
├── <submodule2>/
└── tests/
    └── test_<module>.py  # Mirror tests in top-level tests/
```

**But note**: METAINFORMANT places module tests in the top-level `tests/` mirroring the package structure, not inside each module. So you'd add `tests/<module>/test_*.py`.

### Steps to Add a New Module

1. **Create directory** under `src/metainformant/`: `mkdir -p src/metainformant/mymodule`
2. **Add `__init__.py`** to define public API:
   ```python
   """Brief module docstring."""
   from . import submodule1, submodule2
   __all__ = ["submodule1", "submodule2"]
   ```
3. **Add `AGENTS.md`** – AI agent directives (see [docs/AGENTS.md](../docs/AGENTS.md) for template)
4. **Add `SPEC.md`** – technical specification: functions, classes, types, error codes
5. **Add `README.md`** – user guide with examples (follow existing module READMEs as templates)
6. **Write code** in subpackages (e.g., `analysis/`, `io/`, `utils/`)
7. **Write documentation** in `docs/<domain>/`:
   - Add or update `docs/<domain>/index.md`
   - Possibly create new domain directory if module introduces a new domain
8. **Add tests** in `tests/<module>/test_*.py`
9. **Update `pyproject.toml`** if necessary (dependencies, optional extras)
10. **Update root documentation**:
    - `README.md` module table if new top-level module
    - `docs/index.md` module list if applicable
11. **Update CI** (`.github/workflows/`) only if new test matrix needed (rare)

### Submodule Conventions

- Use small, focused subpackages (e.g., `dna/sequence/`, `dna/alignment/`, `visualization/plots/`)
- Each submodule should have a clear purpose and public API (define in its `__init__.py`)
- Follow the **no mocking** policy: every function must do real work
- Use `metainformant.core` utilities:
  - `metainformant.core.io` for file I/O
  - `metainformant.core.utils.logging` for logging
  - `metainformant.core.paths` for path manipulation

### Dependency Management

- Core dependencies (numpy, pandas, matplotlib, biopython, scipy, requests, pyyaml, networkx, scikit-learn, statsmodels, seaborn, psutil, rich) are always installed.
- Optional dependencies should be declared in `[project.optional-dependencies]` in `pyproject.toml`:
  - `scientific` – heavy numeric packages (numba, umap-learn, scanpy, anndata)
  - `ml` – ML stacks (optuna, xgboost, lightgbm)
  - `networks` – graph analysis (python-louvain, cdlib)
  - `singlecell` – single-cell analysis packages
  - `visualization` – extra plotting (plotly, bokeh, altair, graphviz)

Users install extras as needed: `uv pip install -e ".[ml,singlecell]"`.

---

## Debugging

### Using pdb / breakpoint()

Python 3.11+ includes built-in `breakpoint()`:

```python
def my_function(x: int) -> int:
    breakpoint()  # Execution will stop here
    return x * 2
```

Run with `uv run python script.py`; when hit, you'll enter PDB.

For pytest:
```bash
uv run pytest --pdb   # Drop into pdb on test failure
uv run pytest --trace  # Start trace on first test
```

### Logging Debug Info

Set logging level to DEBUG at program start:
```python
from metainformant.core.utils.logging import setup_logging
setup_logging(level="DEBUG", file_path="output/debug.log")
```

Or via environment variable (if supported):
```bash
export METAINFORMANT_LOG_LEVEL=DEBUG
uv run python script.py
```

### Profiling Performance

For slow code:
```bash
# CPU profiling with cProfile
uv run python -m cProfile -o output/profile.prof script.py

# Visualize with snakeviz
uv pip install snakeviz
uv run snakeviz output/profile.prof
```

Or use `pytest-benchmark`:
```bash
uv run pip install pytest-benchmark
uv run pytest --benchmark-only tests/my_module/
```

---

## Git Workflow

- **main** – always stable, deployable
- Feature branches: `git checkout -b feat/<module>-<description>`
- Bugfix branches: `git checkout -b fix/<module>-<shortdesc>`
- Docs branches: `git checkout -b docs/<area>-update`

Commit messages: see [Conventional Commits](CONTRIBUTING.md#commit-messages).

Before pushing:
```bash
git fetch origin
git rebase origin/main          # rebase onto latest main
uv run pre-commit run --all-files
uv run scripts/run_tests.sh
```

---

## CI/CD

GitHub Actions runs on every push/PR:
- **Python 3.11, 3.12** test matrix
- **Lint**: black, isort, flake8
- **Type check**: mypy
- **Tests**: pytest with coverage
- **Docs build**: Sphinx must build without errors

CI status appears on the PR. All checks must pass before merge.

---

## Finding Your Way

- **Architecture overview**: [docs/architecture.md](architecture.md)
- **Module index**: [docs/index.md](index.md)
- **Tutorials**: [docs/TUTORIALS.md](TUTORIALS.md)
- **Task reference**: [docs/tasks/](../docs/tasks/)
- **FAQ**: [docs/FAQ.md](FAQ.md)
- **Troubleshooting**: [docs/TROUBLESHOOTING.md](TROUBLESHOOTING.md)

---

## Asking for Help

- **Documentation**: Edit any `.md` file directly and PR — no gatekeeping
- **Code questions**: [GitHub Discussions](https://github.com/docxology/metainformant/discussions)
- **Chat**: `#metainformant:matrix.org`
- **Issues**: [GitHub Issues](https://github.com/docxology/metainformant/issues) - use bug template

Welcome aboard! Your contributions make METAINFORMANT better for everyone.
