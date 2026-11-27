# Offline table reader

The **MetaInformAnt Offline Reader** covers local table I/O in version 0.4.0. It checks the schema and values of a synthetic TSV, writes a second TSV with the package API, compares both parsed tables, and charts the fixture. The [historical manuscript status](../MANUSCRIPT_STATUS.md) records the earlier project context. The [documentation index](../../index.md) covers the broader package.

The [reader configuration](config.yaml) orders the [scope](01_scope.md) and [table check](02_table_check.md). The [bibliography](references.bib) and [synthetic TSV](synthetic_counts.tsv) are maintained beside the prose.

## Reproduce the check and render

Use Python 3.14 with `matplotlib==3.11.2`, `psutil==7.2.2`, `PyYAML==6.0.3`, and `pytest==9.1.1` available in the active environment. From the MetaInformAnt checkout, set `READER_ROOT` to a new absolute path ending in `MetaInformAnt`:

```bash
export PYTHONPATH="$PWD/src"
python docs/manuscript/reader/export.py "$READER_ROOT"
python -m pytest --noconftest -o addopts= -q \
  tests/core/test_core_io.py::test_tsv_roundtrip \
  tests/core/test_core_io.py::test_json_roundtrip \
  tests/core/test_core_io.py::test_jsonl_roundtrip
```

The exporter writes the manuscript sources under `docs/manuscript/` in the new directory, the round-tripped table and its numeric JSON summary under `output/data/`, and the chart under `output/figures/`. It uses only local files. The test selection covers three local storage round trips.

From a separate [Template checkout at the verified pin](https://github.com/docxology/template/tree/87704feae0a4490ef7dad9c931d268e8c30273c2) with its Python environment initialized and Pandoc, pandoc-crossref, and LaTeX on PATH, register the export and run:

```bash
unset PYTHONPATH
git -C "$READER_ROOT" init --quiet
git -C "$READER_ROOT" add docs/manuscript output/data output/figures
ln -s "$READER_ROOT" projects/working/MetaInformAnt
uv run --no-sync python scripts/pipeline/stage_03_render.py --project working/MetaInformAnt --skip-manuscript-hydration
uv run --no-sync python scripts/pipeline/stage_04_validate.py --project working/MetaInformAnt
```

The combined PDF, HTML, and validation report appear in the exported project’s `output/pdf/`, `output/web/`, and `output/reports/`. This check covers the synthetic reader and the three named I/O tests.
