# AGENTS.md — `MetaInformAnt/.cursor/skills/metainformant-projects-hymenoptera-amalgkit-doc-00-setup`

Generated METAINFORMANT Cursor skill wrapper (one per repo folder that has an `AGENTS.md`; produced by `scripts/package/generate_cursor_skills.py` — verified from disk, 2026-08-30 doc-realization pass). Regenerate after adding/moving any `AGENTS.md`, then run the parity check: `uv run python scripts/package/generate_cursor_skills.py --check`.

## Contents
- `SKILL.md` — the skill: YAML frontmatter (`name`, `description`) plus instructions to read the linked folder docs before editing that subtree.
- `AGENTS.md`, `README.md` — fleet skeletons (this pass replaces them with this file's content).

## Skill target
- Skill name: `metainformant-projects-hymenoptera-amalgkit-doc-00-setup`
- Points to the real `AGENTS.md` at repo path `projects/hymenoptera_amalgkit/doc/00_setup/AGENTS.md` (verified present).
- Related overview: `README.md` link inside `SKILL.md` → `projects/hymenoptera_amalgkit/doc/00_setup/README.md` (unverified).

## Gotchas
- Do not hand-edit skill wrappers for content changes — edit the target folder's `AGENTS.md`/`README.md` and regenerate.
- Skill is a thin pointer only; real rules live in the target folder docs and root `CLAUDE.md` (uv only, `output/` for outputs, `.tmp/` for temp, real implementations).
