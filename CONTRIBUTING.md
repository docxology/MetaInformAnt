# Contributing to METAINFORMANT

Thank you for your interest in improving METAINFORMANT! This guide covers contribution workflows, standards, and review process.

## Quick Start (First Contribution)

1. **Fork** the repository on GitHub
2. **Clone** your fork locally
3. **Create branch**: `git checkout -b my-feature`
4. **Install dev dependencies**: `uv pip install -e ".[dev,test,docs]"`
5. **Make changes** following [Project Standards](#coding-standards)
6. **Run tests**: `scripts/run_tests.sh` (must pass)
7. **Run lint**: `pre-commit run --all-files`
8. **Commit** (clear message → [conventional commits](https://www.conventionalcommits.org/))
9. **Push** → open PR on GitHub

## Areas of Contribution

| Area | What's Needed | Contact |
|------|---------------|---------|
| **New modules** | Bioinformatics domain (epigenomics, proteomics, metabolomics) | @maintainers |
| **Documentation** | Tutorials, example code, module guides | @docs-team |
| **Performance** | NumPy/SciPy optimization, GPU support | @perf-team |
| **Cloud** | AWS/Azure support, Kubernetes, workflow engines | @cloud-team |
| **LLM Integration** | New MCP tools, Claude/Cursor plugins | @ml-team |
| **Bug fixes** | Issues tagged `bug` or `help wanted` | Any maintainer |

## Coding Standards

### Python

- **Python 3.11+ only** (no 3.10 compatibility)
- **Type hints required** for all public functions (PEP 484)
- **Docstrings**: Google style — parameters, returns, raises
- **Line length**: 100 chars (not 79) — [black](https://github.com/psf/black) will format
- **Imports**: Standard lib → third-party → local, sorted by isort

```python
# Good
from typing import List
import numpy as np
from metainformant.dna import sequences

# Bad
import numpy as np, from metainformant.dna import sequences
```

### Tests

- **Location**: `tests/<module>/test_<feature>.py`
- **One assertion per test** (one concern per test)
- **Fixtures** in `tests/conftest.py` (reusable)
- **Run**: `scripts/run_tests.sh` (CI-parity wrapper)
- **No mocks** for internal functions — use real code ([NO_MOCKING_POLICY.md](NO_MOCKING_POLICY.md))

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat(rna): add differential expression for multi-factor designs

fix(gwas): correct allele frequency calculation in plink reader

docs(cloud): improve gcloud installation instructions

refactor(core): replace pandas with polars in io module

test(dna): add coalescent simulation coverage

chore: cleanup obsolete import in __init__.py
```

**Format:** `<type>(<scope>): <subject>`
- `type`: `feat`, `fix`, `docs`, `test`, `refactor`, `chore`, `build`, `ci`
- `scope`: module name (`rna`, `gwas`, `core`, `cloud`)
- `subject`: imperative mood, no period, max 72 chars

## Documentation Requirements

Every feature contribution MUST include documentation update:

1. **Module README.md** (if touching module code)
 - Add usage examples for new function
 - Update "Key Functions" table with signature
 - Add "See Also" cross-reference if related to other modules

2. **docs/<module>/index.md** or **docs/<domain>/index.md**
 - Add 1–2 paragraph overview of feature
 - Include code snippet showing typical usage
3. **docs/TUTORIALS.md** (optional but preferred for user-facing features)
 - Add section with step-by-step walkthrough
 - Include expected output (console + figures)
4. **docstring** in source code (Google style)

**No docs → no merge** (enforced by CI).

## Pull Request Process

1. **Open draft PR** early — get feedback before code complete
2. **CI must pass** — all tests green, lint clean
3. **Address review comments** within 48h
4. **Update CHANGELOG.md** (if user-facing change) — bullet under `[Unreleased]`
5. **Request review** from @maintainers team
6. **Merge** only after at least 1 approval

### PR Checklist

- [ ] Tests added (or existing tests updated)
- [ ] `scripts/run_tests.sh` passes locally
- [ ] No new `TODO` comments (or `TODO(#issue)` with issue number)
- [ ] Documentation updated (README + docs/)
- [ ] Docstrings complete (Google style)
- [ ] Type hints present (mypy clean)
- [ ] Pre-commit hooks pass (`pre-commit run --all-files`)
- [ ] No dead code (`vulture` optional but encouraged)
- [ ] Branch rebased onto latest `main`

## Reporting Bugs

**Before filing:** Check existing [issues](https://github.com/docxology/metainformant/issues) and [FAQ](docs/FAQ.md).

Use the bug template — include:
1. **Python version**, `metainformant.__version__`
2. **OS** (Linux/macOS/Windows + distro)
3. **Minimal reproducible code** (ideally <20 lines)
4. **Error traceback** (full, not partial)
5. **Expected vs actual behavior**

## Feature Requests

- **Search existing requests** first (don't duplicate)
- **Provide use case**: "I need X because Y"
- **Include example input/output** if applicable
- **Consider scope**: Smaller features get faster review

## Community

- **Discussions**: GitHub Discussions (Q&A, design questions)
- **Chat**: Matrix room #metainformant:matrix.org (real-time help)
- **Issues**: GitHub Issues (bug reports, feature requests)
- **Docs**: Edit any `.md` file and PR — no gatekeeping

## License

By contributing, you agree your work is licensed under the project's [Apache 2.0](LICENSE) license.

---

**Questions?** Open a [Discussion](https://github.com/docxology/metainformant/discussions) or ping @maintainers.

**Thanks** to all 200+ contributors! 
