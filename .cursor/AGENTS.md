# AGENTS.md — `MetaInformAnt/.cursor/`

Cursor editor configuration for METAINFORMANT. Contains `skills/` — 156
generated per-directory skill folders plus a `README.md`. Skills are
**generated** by `scripts/package/generate_cursor_skills.py` from the nested
`AGENTS.md` files (verify with `--check`); never hand-edit a skill folder.

## Children
- `skills/` — generated skill folders `metainformant-<path-slug>` mapping
  1:1 to repo directories (e.g. `metainformant-scripts-eqtl` → `scripts/eqtl`),
  each holding a `SKILL.md` that links to the target directory's AGENTS.md.
  Three hash-named folders are stale generated artifacts for amalgkit doc
  paths (`doc/01_infrastructure`, `doc/03_data_management`,
  `doc/04_troubleshooting`).
Repo-wide policy: see `/Volumes/external_drive/Git/template/projects/ongoing/AGENTS.md`.