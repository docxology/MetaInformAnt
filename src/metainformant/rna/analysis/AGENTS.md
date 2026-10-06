# Agent guidance: RNA analysis methods

Own reusable matrix/profile, orthology, normalization, comparative-statistics
and analysis-contract methods. Project adapters retain cohort paths, manifests,
artifact assembly and project readiness gates. Parent methods import package APIs,
not nested project scripts.

Before changing an estimand or input contract, trace callers and update real
numerical/refusal tests, [README](README.md) and [SPEC](SPEC.md). Keep fingerprint,
ortholog mean-profile and sample-aligned gene distances distinct. Unavailable
pairs remain explicit; API availability does not establish executed project results.
Use canonical absolute package imports and the parent frozen `uv` environment.
Repo-wide guidance is in the repository-root `AGENTS.md`.
