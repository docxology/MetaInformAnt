# AGENTS.md — `MetaInformAnt/.github/`

GitHub configuration for METAINFORMANT (verified 2026-11-13): `workflows/` (CI:
`test.yml`, `build.yml`, `release.yml`, `security.yml`, plus a workflows
`README.md`). The Python floor is 3.14; the workflows README lists the Python
versions of each job.

## Gotchas
- CI runs the repo's no-mocks test policy; do not add mock-based tests.
Repo-wide policy: see the repository-root `AGENTS.md`.