# Software-methods manuscript

**MetaInformAnt Core Table and Metadata Methods** documents public table I/O and validation interfaces using independent checks and deterministic controls. It is a software paper with synthetic evidence. The [historical manuscript status](MANUSCRIPT_STATUS.md) describes prior project context. The [offline reader](reader/README.md) remains available for a shorter worked example.

The canonical [configuration](config.yaml) orders the [abstract](00_abstract.md), [introduction](01_introduction.md), [architecture](02_architecture.md), [methods](03_methods.md), [results](04_results.md), [discussion](05_discussion.md), and [availability statement](06_availability.md). Citations are in [references.bib](references.bib). The [documentation index](../index.md) links to package guides.

## Reproduce the paper

Use Python 3.14 with `matplotlib==3.11.2`, `psutil==7.2.2`, `PyYAML==6.0.3`, and `pytest==9.1.1` available in the active environment. From a clean MetaInformAnt checkout, set `PAPER_ROOT` to a new absolute path ending in `MetaInformAnt-methods`:

```bash
export PYTHONPATH="$PWD/src"
python docs/manuscript/export_methods.py "$PAPER_ROOT"
python -m pytest --noconftest -o addopts= -q \
  tests/core/test_core_io.py::test_tsv_roundtrip \
  tests/core/test_core_io.py::test_json_roundtrip \
  tests/core/test_core_io.py::test_jsonl_roundtrip
```

The export copies the ordered sources and runs the benchmark. The generated `output/data/` holds all fifteen case tables and `table_contract.json`, which records expected and observed outcomes, valid counts, artifact sizes, and SHA-256 hashes. The chart is written to `output/figures/table_contract.png`.

In a separate [Template checkout at the verified pin](https://github.com/docxology/template/tree/87704feae0a4490ef7dad9c931d268e8c30273c2) with its Python environment initialized and Pandoc, pandoc-crossref, and LaTeX on PATH, register the export using an unused working alias and run:

```bash
unset PYTHONPATH
git -C "$PAPER_ROOT" init --quiet
git -C "$PAPER_ROOT" add docs/manuscript output/data output/figures
ln -s "$PAPER_ROOT" projects/working/MetaInformAnt-methods
uv run --no-sync python scripts/pipeline/stage_03_render.py --project working/MetaInformAnt-methods --skip-manuscript-hydration
uv run --no-sync python scripts/pipeline/stage_04_validate.py --project working/MetaInformAnt-methods
```

The PDF, HTML, validation report, and rendered provenance receipt appear under the export’s `output/` directory. The checks establish the behavior of the table and metadata interfaces exercised here.
