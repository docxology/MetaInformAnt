# Testing

Policy

- All tests must write artifacts only under `output/`. Avoid networked writes and other side effects unless explicitly required by an integration test.
- Prefer deterministic seeds and stable filenames to enable reproducible CI diffs.

Running

```bash
pytest -q
```

Conventions in tests

- Use `metainformant.core.io.ensure_directory` to create subfolders under `output/`.
- Keep input fixtures under `tests/data/` or `data/` as appropriate.

Run tests via CLI wrapper or directly with pytest.

```bash
metainformant tests -q
# or
pytest -q
```

Structure
- Unit tests cover core utilities, DNA modules (sequences/alignment/msa/phylogeny/popgen), RNA configs/workflow, simulation, visualization, and domain availability.

Write tests first (TDD), then implement or update functions. See `tests/` for examples.

Related: [CLI](./cli.md), [Core](./core.md)


