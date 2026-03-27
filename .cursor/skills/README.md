# Project Cursor skills

This directory holds **Agent Skills** for Cursor: each subfolder contains a `SKILL.md` with YAML front matter (`name`, `description`) and a short body.

- **Source of truth** for technical detail remains each folder’s `AGENTS.md` (and `README.md` where present). Skills are thin entry points so the agent loads the right subtree rules.
- **Regenerate** after adding or moving `AGENTS.md`:

  ```bash
  uv run python scripts/package/generate_cursor_skills.py
  ```

- **Verify** skills match every `AGENTS.md`:

  ```bash
  uv run python scripts/package/generate_cursor_skills.py --check
  ```

Do not copy long documentation into `SKILL.md`; link to `AGENTS.md` and repo guides instead.
