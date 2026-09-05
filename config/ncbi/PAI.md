# PAI - config/ncbi

## Context & Intent

`config/ncbi/` holds NCBI Entrez API configuration shared by every module that accesses NCBI services. Per `AGENTS.md`: "NCBI API and data retrieval configuration (email, API keys, rate limiting)." It holds one file, `ncbi.yaml`, with the email, rate limiting, and retry settings loaded through `src/metainformant/core/utils/config.py`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer).
- Cross-cutting scope: unlike sibling dirs, `ncbi.yaml` is not module-specific — it serves all domains hitting NCBI (SRA fetch, assembly lookup, etc.) via the sequential-failover pattern in root `SPEC.md` (local sources before remote NCBI/SRA acquisition).
- Downstream consumer: `metainformant.core.utils.config` (explicitly linked in `README.md`).

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the NCBI namespace documented by the module."

From `README.md`: `NCBI_EMAIL` overrides the `email` setting — real settings in `ncbi.yaml` are `email: ""` (set your own), `rate_limit_delay: 0.34` (NCBI caps at 3 requests/second), `max_retries: 3`, `retry_delay: 1.0`.

## AI Workflows

- **Editing `ncbi.yaml`**: preserve the four documented keys; never commit a personal email — rely on `NCBI_EMAIL` for user-specific values.
- **Verification**: load via `metainformant.core.utils.config` under `uv run python` and confirm the mapping keys; targeted core-config tests only, never repo-wide (shared-instruction USB constraint).
- **Scope discipline**: config data only; the retry/rate-limit logic lives in `core/` and must not be modified from here.
