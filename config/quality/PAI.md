# PAI - config/quality

## Context & Intent

`config/quality/` is the configuration slot for the quality module (QC metrics, contamination detection). Its job is to hold YAML/JSON config that the quality module loads through `metainformant.core.utils.config`. Today the directory is documentation-only: no YAML config file exists yet; the only non-doc file is `mypy_error_budget.txt`, a small counter maintained by typed-CI tooling.

## Virtual Hierarchy

- Parent: `config/` — the repo-wide YAML configuration layer described in root `SPEC.md`.
- Sibling config dirs: `eqtl/`, `singlecell/`, `networks/`, `multiomics/`, `longread/`, `ncbi/`, `gwas/`, `amalgkit/`, among others.
- Downstream consumer: `src/metainformant/quality/` (business logic) — config here must never contain logic, only data.

## Maintenance Notes

From `AGENTS.md` (this directory):

- Repo-wide policy lives in the repository-root `AGENTS.md`; follow it.
- "YAML/JSON config for the quality module" is the directory's entire role — verified 2026-08-29 as 1 file, 0 subdirs.
- Root rules that apply here: dependency management via `uv` only; REAL IMPLEMENTATION policy (no mocks); all outputs go to `output/` (root `PAI.md`).
- Do not add Python files here; `config/` holds data only.

## AI Workflows

- **Adding a config**: create `config/quality/<name>.yaml` following the section style of sibling templates (e.g. `config/singlecell/singlecell_template.yaml`), validate against the quality module's schema before committing, then update `README.md` and `SPEC.md` in the same change.
- **Verification**: load it with `uv run python` and `metainformant.core.utils.config.load_mapping_from_file`, confirming sorted keys — the same pattern the sibling `config/eqtl/README.md` demonstrates.
- **Testing**: quality-module tests live under `tests/quality/`; run them targeted, never repo-wide (USB-volume constraint in shared instructions).
