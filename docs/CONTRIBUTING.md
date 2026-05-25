# Contributing to METAINFORMANT

See the root-level [CONTRIBUTING.md](../CONTRIBUTING.md) for complete contribution guidelines, coding standards, and PR process.

Quick links:
- [Development setup](SETUP.md#development-mode)
- [Testing guide](testing.md) (comprehensive testing documentation)
- [Documentation guide](DOCUMENTATION_GUIDE.md)
- [Real Implementation policy](REAL_IMPLEMENTATION_POLICY.md)

All contributions require documentation updates alongside code changes.

---

## Comprehensive Contributor Guide

This guide provides detailed information for contributing to METAINFORMANT, including code style, testing, PR process, documentation standards, and community guidelines. For a quick start, see the first-contribution steps in the [root CONTRIBUTING.md](../CONTRIBUTING.md).

---

### Code Style and Formatting

METAINFORMANT enforces strict code quality standards to maintain a clean, maintainable codebase.

#### Python Version
- **Python 3.11+ only** (no backward compatibility with 3.10)
- All code must use `from __future__ import annotations` for forward reference support

#### Formatting: Black
- Code formatting is handled by [Black](https://github.com/psf/black) with line length 120
- Black is automatically run via pre-commit hooks
- Configuration in `pyproject.toml`: `line-length = 120`, target versions `py311` and `py312`
- Manual invocation: `black .` or `uv run black .`

#### Import Sorting: isort
- Imports must be sorted using [isort](https://github.com/pycqa/isort) with the Black profile
- Configuration: `profile = "black"`, `line_length = 120`
- Pre-commit hook runs isort automatically
- Manual: `isort .` or `uv run isort .`

#### Type Hints: mypy
- **All public functions must have type hints** (PEP 484 style)
- mypy configuration in `pyproject.toml` is strict:
  - `disallow_untyped_defs = true`
  - `disallow_incomplete_defs = true`
  - `check_untyped_defs = true`
  - `warn_unused_ignores = true`
  - `strict_equality = true`
- Some third-party libraries (Bio, matplotlib, pandas, numpy, sklearn) have missing imports ignored per-module
- Run manually: `mypy src/metainformant` or `uv run mypy src/metainformant`

#### Linting: flake8 (optional but encouraged)
- flake8 is configured with `max-line-length = 120`
- Some errors are ignored (E203, W503)
- Pre-commit includes flake8 hook but it may be disabled in the current configuration

#### Docstrings
- **Google style** docstrings for all public functions, classes, and modules
- Include: `Args:`, `Returns:`, `Raises:` sections as applicable
- Example:
```python
def gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as a fraction (0.0 to 1.0)
    """
    ...
```

---

### Pre-commit Hooks

METAINFORMANT uses [pre-commit](https://pre-commit.com) to automate code quality checks. The configuration is in `.pre-commit-config.yaml` and includes:

- **black** – code formatting (Python 3.12, line-length 120)
- **isort** – import sorting (Black profile)
- **flake8** – linting (currently disabled due to config issues)
- **mypy** – type checking (currently disabled due to many issues, but encouraged locally)
- **General file checks** (via `pre-commit-hooks`):
  - `check-added-large-files` (max 1000 KB)
  - `check-merge-conflict`
  - `check-yaml`, `check-json`, `check-toml`
  - `end-of-file-fixer`
  - `trailing-whitespace` (excludes `.md`, `.rst`)
  - `detect-private-key`
- **Local hooks**:
  - `validate-examples` – validates example Python files
  - `detect-aws-credentials`
  - `mixed-line-ending` – enforce LF
- **pygrep-hooks**: Python-specific checks (blanket noqa, type annotations, no log warnings, RST syntax)
- **Disabled**: bandit (security), pydocstyle (documentation style) – can be enabled locally if desired

#### Using Pre-commit

Install pre-commit (if not already): `uv pip install pre-commit`

Run once to install git hooks: `pre-commit install`

Run manually on all files: `pre-commit run --all-files`

Hooks run automatically on `git commit`. The `fail_fast` setting is `false`, so all hooks execute even if one fails.

---

### Testing Requirements

METAINFORMANT has a **strict real-implementation policy** and emphasizes real-data, integration-style testing.

#### Test Framework: pytest
- Tests reside in `tests/` mirroring the `src/metainformant/` structure
- Test files named `test_<module>.py`
- Run tests via `scripts/run_tests.sh` (CI-parity wrapper) or directly: `uv run pytest`
- Add options: `-v` (verbose), `-x` (stop on first failure), `-m "not slow"` (exclude slow tests), `--timeout=300` (per-test timeout)

#### Real Implementation Policy
**Read [REAL_IMPLEMENTATION_POLICY.md](REAL_IMPLEMENTATION_POLICY.md)** for the complete policy. Key points:

- **ALL functions must perform real computations or genuine external calls** (UniProt, PDB, NCBI, etc.)
- **Zero tolerance for placeholders**: no `return [[0,1,2]*100]`, no dummy data, no `return None` inert plots
- **Tests use real data or generated synthetic data** (via `numpy.random`, actual algorithms)
- The purpose: ensure every function does real work; avoid false confidence

#### Coverage
- Coverage is collected when running `pytest --cov=src/metainformant`
- Configuration in `pyproject.toml` under `[tool.coverage]`
- Aim to maintain high coverage (>90% recommended). PRs should not reduce coverage.
- Coverage reports: terminal (`--cov-report=term-missing`), HTML (`--cov-report=html`), XML (`--cov-report=xml` for CI)

#### Test Organization
- `tests/conftest.py` – shared fixtures
- Markers: `slow`, `network` (real API calls), `external_tool` (requires external binaries), `integration`, `real_impl`
- Parallel execution via `pytest-xdist`: `uv run pytest -n auto`

---

### Pull Request Process

Follow these steps for a smooth PR experience:

1. **Create a branch** from `main`: `git checkout -b my-feature`
   - Use descriptive names: `feat/dna-msa-improvements`, `fix/gwas-plot-bug`, `docs/update-tutorials`
2. **Make changes** following [Coding Standards](#code-style-and-formatting)
3. **Run locally**:
   - `scripts/run_tests.sh` – ensures all tests pass
   - `pre-commit run --all-files` – ensures lint/formatting is clean
4. **Commit** with a clear message following [Conventional Commits](https://www.conventionalcommits.org/):
   ```
   feat(dna): add progressive alignment algorithm
   fix(visualization): correct axis limits in manhattan plot
   docs(rna): update amalgkit workflow documentation
   ```
   - Format: `<type>(<scope>): <subject>`
   - Types: `feat`, `fix`, `docs`, `test`, `refactor`, `chore`, `build`, `ci`
   - Scope: module name (`dna`, `rna`, `gwas`, `core`, `visualization`)
5. **Push** to your fork and open a **Draft PR** early to get feedback.
6. **CI pipeline** runs automatically:
   - Tests on multiple Python versions (3.11, 3.12)
   - Linting (black, isort, flake8)
   - Type checking (mypy)
   - Documentation build check
7. **Address review comments** from maintainers within 48 hours.
8. **Update CHANGELOG.md** under `[Unreleased]` if your change is user-facing.
9. **Request review** from the `@maintainers` team once CI passes.
10. **Merge** only after at least one approval and CI green.

#### PR Checklist
Before marking your PR ready, verify:
- [ ] Tests added (or existing tests updated)
- [ ] `scripts/run_tests.sh` passes locally
- [ ] No new `TODO` comments (or `TODO(#issue)` referencing an issue)
- [ ] Documentation updated (README in module, `docs/<domain>/`, `docs/TUTORIALS.md` if applicable)
- [ ] Docstrings complete (Google style)
- [ ] Type hints present (`mypy src/metainformant` clean)
- [ ] Pre-commit hooks pass (`pre-commit run --all-files`)
- [ ] No dead code (`vulture` optional but encouraged)
- [ ] Branch rebased onto latest `main`

---

### Documentation Standards

Every code contribution **must** include documentation updates. Incomplete documentation PRs will not be merged.

#### Required Documentation Updates

1. **Module README.md** (in `src/metainformant/<module>/README.md`)
   - Add usage examples for new functions
   - Update the "Key Functions" table with signature and brief description
   - Add "See Also" cross-references to related modules

2. **docs/<domain>/index.md** or **docs/<domain>/topic.md**
   - Write 1–2 paragraphs explaining the feature
   - Include a runnable code snippet showing typical usage
   - Cross-link to the module README and related topics

3. **docs/TUTORIALS.md** (optional but strongly preferred for user-facing features)
   - Add a section with step-by-step walkthrough
   - Include expected output (console output and figures if applicable)

4. **Docstrings** in source code (Google style, see above)

#### Documentation Style
- Clear, technical writing; avoid marketing fluff
- Code examples must be **real and runnable** (copy-paste tested)
- Use consistent Markdown formatting
- Keep cross-references up to date (`[link text](relative/path.md)`)
- Place new documentation in the appropriate `docs/<domain>/` subdirectory

#### Building Documentation
- Sphinx builds the docs from `docs/` using `conf.py`
- To build: `bash scripts/package/uv_docs.sh` or `uv run sphinx-build -b html docs/ docs/_build/html`
- View: open `docs/_build/html/index.html` in a browser
- When adding new pages, update the toctree in the parent `index.md` or relevant `README.md`

---

### Issue Triage and Labels

METAINFORMANT uses GitHub Issues for bug reports and feature requests. Labels help organize and prioritize:

- **bug** – Something is broken; requires a fix
- **enhancement** – New feature or improvement to existing functionality
- **help wanted** – Community contributions welcome (good first issues)
- **good first issue** – Suitable for new contributors
- **documentation** – Documentation improvements (no code changes)
- **question** – Usage questions (use GitHub Discussions instead)
- **duplicate** – Already reported elsewhere
- **wontfix** – Won't be addressed (closed)
- **invalid** – Not a valid issue (closed)

#### Triage Workflow
- Maintainers label new issues within 48 hours
- `help wanted` and `good first issue` labels indicate opportunities for community contribution
- Before working on an issue, comment expressing interest; assign if approved

---

### Community Guidelines

METAINFORMANT is an open, welcoming community. We adhere to the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/). Please read it.

#### Principles
- **Be respectful and inclusive** – no harassment, discrimination, or hostile behavior
- **Constructive collaboration** – provide helpful feedback, accept critique gracefully
- **Focus on the code** – debates should be technical and evidence-based
- **Help newcomers** – guide new contributors patiently

#### Communication Channels
- **GitHub Issues** – bug reports, feature requests (tracked)
- **GitHub Discussions** – Q&A, design discussions, community chat
- **Matrix Room**: `#metainformant:matrix.org` – real-time help and chat
- **Documentation** – any `.md` file can be edited via PR; no gatekeeping

---

### Quick Reference

| Topic | Reference |
|-------|-----------|
| Full contribution guide | [Root CONTRIBUTING.md](../CONTRIBUTING.md) |
| Setup instructions | [docs/SETUP.md](SETUP.md) |
| Testing guide | [docs/testing.md](testing.md) |
| Documentation guide | [docs/DOCUMENTATION_GUIDE.md](DOCUMENTATION_GUIDE.md) |
| Real-implementation policy | [docs/REAL_IMPLEMENTATION_POLICY.md](REAL_IMPLEMENTATION_POLICY.md) |
| CI/CD workflows | `.github/workflows/` |
| Code of Conduct | [Contributor Covenant](https://www.contributor-covenant.org/) |
| Community chat | `#metainformant:matrix.org` |

---

**Questions?** Open a [GitHub Discussion](https://github.com/docxology/metainformant/discussions) or ping `@maintainers`.
