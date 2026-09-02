# Cursor / Hermes skills surface

`.cursor/skills/` holds 419 generated wrapper skills, one per `AGENTS.md`
location in the repository. They are produced deterministically by
[`scripts/package/generate_cursor_skills.py`](../scripts/package/generate_cursor_skills.py)
— never hand-edit a wrapper; regenerate instead.

## Regenerate and verify

```bash
uv run python scripts/package/generate_cursor_skills.py          # write
uv run python scripts/package/generate_cursor_skills.py --check  # exit 1 on drift
uv run pytest tests/scripts/test_generate_cursor_skills.py -q    # regression tests
```

`--check` verifies, for every live `AGENTS.md`:

- a matching `SKILL.md` exists with correct frontmatter (`name` == slug),
- byte-exact content matches regeneration (no drift),
- all canonical links resolve,
- no orphan skill directories (a `SKILL.md` without an `AGENTS.md` source),
- module-surface claims are true: for `src/metainformant/<module>` wrappers,
  every backtick-quoted submodule exists on disk and every relative import in
  the package `__init__.py` resolves statically (AST-based, import-free).

## Naming

Slugs are derived from the source path (`src/metainformant/rna` →
`metainformant-src-metainformant-rna`). Paths whose readable slug would exceed
64 characters are compressed by dropping leading segments with prefix-aware
deduplication (e.g. `src/metainformant/structural_variants/visualization` →
`metainformant-structural-variants-visualization`). A sha256 digest slug is the
last-resort fallback only when even a single segment cannot fit; a regression
test fails if any digest-named skill directory appears.

## Enriched module skills

`src/metainformant/<module>` wrappers carry a generated **Module surface**
section: purpose (from the package docstring), public submodules (from
`__all__` / relative imports, verified on disk), the canonical import, and the
test entry point (`uv run pytest tests/<module> -q`, one pytest directory per
invocation).

## How the set surfaces in Hermes

Hermes-style agent runtimes discover skills by scanning skill roots for
`SKILL.md` frontmatter (`name`, `description`) — the same convention as Cursor.
For this repo the discovery root is `.cursor/skills/`; each generated wrapper
appears as one entry in `list_skills`-style inventories keyed by its slug.

To sync a fresh checkout into such a runtime:

1. `uv run python scripts/package/generate_cursor_skills.py` (regenerate if
   `--check` reports drift),
2. point the runtime's skill-search roots at this checkout's `.cursor/skills/`
   (Hermes loads skills from configured roots per profile; adding this path to
   the profile's skill roots exposes all 419 wrappers),
3. re-run the runtime's skill listing (Hermes: ask for the skill list or run
   the profile's `list_skills` surface) and confirm the count matches
   `--check` output (`check ok: N skills match AGENTS.md locations`).

Because every wrapper is byte-deterministic, regeneration produces identical
output on any machine at the same commit; syncing is therefore a pure
directory mount, not a conversion step.

## Related

- [`scripts/package/AGENTS.md`](../scripts/package/AGENTS.md)
- [`docs/REAL_IMPLEMENTATION_POLICY.md`](REAL_IMPLEMENTATION_POLICY.md)
