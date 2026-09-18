# Security and dependency check of the release surface — 2026-09-17

## Verified outcome

**PASS — gate established; one medium-severity finding dispositioned to the
baseline; no known dependency vulnerabilities.** This record documents the
first run of the blocking security and dependency gate
(`scripts/package/security_checks.sh`, wired as the blocking
`security-checks` job in `.github/workflows/security.yml`). All scanners ran
read-only via `uvx` against the repository at commit `3f48bb2d7` (clean tree,
branch `main`); no environment, `pyproject.toml`, or `uv.lock` change was made
or required. The live RNA/Amalgkit producer was not touched.

This record is an engineering audit ledger, not a scientific result. It makes
no significance, causal, or inference claims beyond gate mechanics.

## Scope and tooling

| Check | Tool invocation | Surface |
| --- | --- | --- |
| Static security analysis | `uvx bandit -q -r src/metainformant -f json` | `src/metainformant` (tests excluded; bandit default excludes) |
| Dependency vulnerabilities | `uv export --frozen --no-emit-project --format requirements-txt --no-hashes` then `uvx pip-audit -r <export>` | All 52 pinned distributions derivable from `uv.lock` |
| Shell-execution / path-traversal pattern audit | `grep -rnE` over `src/metainformant` (`shell=True`, `os.system(`, `os.popen(`, bare `eval(`/`exec(`, string literals beginning `../`) | Maintained `src/metainformant` Python sources, `__pycache__`/`tests` excluded |

## Observed findings

### Bandit (365 findings: 1 MEDIUM, 364 LOW)

| Severity | Count | Blocking under gate policy |
| --- | --- | --- |
| MEDIUM | 1 | Yes — dispositioned below |
| LOW | 364 | No — reported for review, never blocks |

The single medium-or-higher finding:

| Severity | ID | Location | Observation | Disposition |
| --- | --- | --- | --- | --- |
| MEDIUM | B310 (`urllib_urlopen` blacklist) | `src/metainformant/gwas/workflow/amellifera_pipeline.py:193` | `urllib.request.urlretrieve(url, gz_path)` where `url` is a hardcoded `https://ftp.ncbi.nlm.nih.gov/...` constant. B310 guards against user-controlled URL schemes; no user input reaches `urlretrieve` here, so the guarded risk does not apply. | Baselined in `scripts/package/security_baseline.txt` (dated 2026-09-17, owner `gwas-maintainers`). Recommended remediation, non-blocking: migrate the call to the hardened download helper in `src/metainformant/core/io/download.py`, which already carries a justified `# nosec B310` with explicit scheme handling. |

The task's in-scope fix area (`src/metainformant/rna/engine` subprocess
argument handling and temp-file handling) required **no changes**: all
observed `subprocess.run`/`Popen` calls there pass argument lists (no
`shell=True`, no string commands), and gzip/kallisto/amalgkit invocations use
constant or validated argv entries (e.g. `src/metainformant/rna/engine/sra_extraction.py:140,166`,
`streaming_orchestrator.py:1440`, `workflow_planning.py:681,725`).

Low-severity distribution by rule (non-blocking, retained for review):

| Rule | Count | Note |
| --- | --- | --- |
| B311 `blacklist` (standard pseudo-random) | 121 | Statistical/simulation code; not a security control |
| B603 `subprocess_without_shell_equals_true` | 104 | All argument-list calls; house style already prefers argv lists |
| B607 `start_process_with_partial_path` | 44 | e.g. `gzip`, `kallisto` resolved via PATH |
| B110 `try_except_pass` | 34 | Error-suppression hygiene, not security |
| B404 `blacklist` (import subprocess) | 26 | Advisory only |
| B101 `assert_used` | 26 | Test/contract assertions |
| B112 `try_except_continue` | 4 | Hygiene |
| B107 `hardcoded_password_default` | 3 | Non-secret default parameter values |
| B403 / B105 `blacklist` / `hardcoded_password_string` | 1 + 1 | Advisory |

### pip-audit (0 known vulnerabilities)

`pip-audit` reported **zero known vulnerabilities across all 52 pinned
distributions** exported from `uv.lock` (`--frozen`, project package excluded).
The dependency surface needs no bump, override, or `uv.lock` modification.

### Pattern audit (shell execution / path traversal)


**Shell-execution patterns (blocking): zero hits** in `src/metainformant`
after filtering pure comment lines (one hit was a comment —
`src/metainformant/gwas/analysis/calling.py:101`, "no shell=True needed" — and
cannot execute). No `shell=True`, no `os.system`/`os.popen`, and no bare
`eval(`/`exec(` on executable code lines were found. ZIP-extraction traversal
hardening from the 2026-08-13 review (`RESEARCH_SOFTWARE_REVIEW_2026-08-13.md`)
remains in place.

**Path-traversal string literals (informational, non-blocking): 10 review
candidates**, all reviewed and assessed benign:

| Location | Observation | Assessment |
| --- | --- | --- |
| `src/metainformant/core/io/paths.py:58` | `dangerous_patterns = ["..//", "..\\", ";", "|", "&", "$"]` | Defensive sanitizer denylist, not traversal use |
| `src/metainformant/gwas/reporting/html_report.py:216-228` (8 hits) | HTML template strings with relative `src="../../*.png"` plot links | Static relative links inside a generated report directory; not filesystem path construction |
| `src/metainformant/ontology/workflow/run_ontology.py:49` | `PROJECT_ROOT / "../../src"` | Repo-root-relative path constant anchored to a computed project root; no external input |

These remain printed on every gate run as `REVIEW (non-blocking)` lines and
archived in `output/security_pattern_traversal.txt`; the blocking shell-exec
class is what fails the gate.
  medium+ bandit findings, any pip-audit vulnerability, and any blocking
  shell-execution pattern-audit hit fail the gate unless covered by
  `scripts/package/security_baseline.txt`.
- `scripts/package/security_baseline.txt` (new): dated, commented, owner-column
  baseline of accepted findings. Currently one entry (B310 above). Adding an
  entry requires a dated justification comment above the
  `SEVERITY|OWNER|ID|LOCATION` line.
- `.github/workflows/security.yml` (new): `Blocking security and dependency
  checks` job on ubuntu-latest, running the script on push and pull request to
  `main` plus `workflow_dispatch`; uploads the security reports (bandit JSON,
  pip-audit JSON, pattern-audit and path-traversal TXT) as artifacts on every
  run. Existing workflows were not modified.

## Boundary

- All scanner invocations were read-only data collection; no gate verdict was
  altered and no source file was edited in response to findings.
- No environment or dependency change: `uvx` ephemeral tools, `uv export`
  to a temp file, no `uv sync`/`install`/`add`, no `pyproject.toml` or
  `uv.lock` edit.
- Nested repositories under `projects/` and `projects/hymenoptera_amalgkit`
  were not read for modification purposes and were not edited.
