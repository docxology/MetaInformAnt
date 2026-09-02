# scripts/test — Per-module test infrastructure

## Purpose

Thin orchestrators for the METAINFORMANT test suite. Business logic for
running tests stays minimal here; the suite itself lives under `tests/`.

## Entry points

### run_module_tests.py — per-module runner

Runs each module's test directory in its **own pytest subprocess**, enforcing
the repo rule that a single pytest invocation must never span multiple
test directories with colliding `conftest` plugins.

```bash
# One module
uv run python scripts/test/run_module_tests.py --module rna

# Several modules, each in its own pytest process
uv run python scripts/test/run_module_tests.py --module rna --module dna

# With per-module coverage of src/metainformant/<module> (cov-append merges)
uv run python scripts/test/run_module_tests.py --module rna --coverage

# Every discovered module, one pytest per module
uv run python scripts/test/run_module_tests.py --all

# Machine-readable results for CI
uv run python scripts/test/run_module_tests.py --module rna --json results.json
```

Exit codes: `0` all passed, `1` any module failed/timed out, `2` usage error.

Per-module timeout defaults to 1800 s (`--timeout` to override); a timeout is
reported as `timed_out: true` with exit code 124 in the JSON payload.

`--list` prints the discovered module names (directories under `tests/` that
contain `test_*.py` and are not shared infrastructure: `_support`, `data`,
`fixtures`, `integration`, `infrastructure`, `other`, `scripts`).

## See also

- `tests/README.md` — test suite map (generated; do not hand-edit the table)
- `tests/AGENTS.md` — test standards and zero-mocks policy
- `scripts/quality/` — mypy budget and other quality gates
