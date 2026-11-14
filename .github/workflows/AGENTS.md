# AGENTS.md — `MetaInformAnt/.github/workflows/`

CI workflows (verified 2026-11-13): `test.yml`, `build.yml`, `release.yml`,
`security.yml`, plus `README.md`. Single-version jobs run Python 3.14; the test
fast lane, build and release install smoke also run 3.15. Changes here affect
GitHub Actions runs only — local doc passes never trigger them.
Repo-wide policy: see the repository-root `AGENTS.md`.
