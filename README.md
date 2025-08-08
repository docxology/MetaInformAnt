# METAINFORMANT: Integrated Multi-Omic Biological Systems Modeling

METAINFORMANT is a unified, modular toolkit for integrated multi-omic analysis across domains:
- DNA: genome retrieval, metadata, and storage, at population and phylogenetic levels
- RNA: transcriptomic metadata, download, quantification, and curation
- Protein: proteome retrieval and descriptive analytics
- Epigenome: DNA modifications (e.g., chromatin/methylation) and other epigenetics pipelines
- Ontology: functional annotation (e.g., GO) and taxonomic utilities
- Phenotype: morphological and behavioral phenotype curation and web scraping (e.g., AntWiki)
- Ecology: ecology metadata, download, curation

This repository organizes practical pipelines and utilities that can be composed within and across omic layers to support systems-level biological modeling.

### Highlights
- Cohesive project layout with clear domain modules
- uv-based environment and dependency management
- Practical pipelines runnable per-domain, designed for orchestration across domains
- Extensible: add new organisms, assays, and analyses incrementally

---

### Project Structure

```
METAINFORMANT/
  src/
    metainformant/
      __init__.py
      __main__.py                 # `python -m metainformant` CLI entry
      core/                       # Shared utilities, I/O, logging, config
        __init__.py
        config.py
        io.py
        logging.py
      dna/
        __init__.py
        genomes.py                # Fetch + parse + register genomes
      rna/
        __init__.py
        pipeline.py               # Metadata → download → quant → curate
      protein/
        __init__.py
        proteomes.py
      epigenome/
        __init__.py
      ontology/
        __init__.py
      phenotype/
        __init__.py
      ecology/
        __init__.py
      legacy/                     # Prior scripts retained for reference
        README.md
  tests/
    test_repo_structure.py        # Asserts new package layout
  metadata/
    metadata.tsv                  # Canonical metadata sink
  scripts/
    setup_uv.sh                   # One-shot uv + venv + sync
  pyproject.toml                  # uv-native dependencies + tool config
```

---

### Installation and Environment (uv)

uv is a fast, modern Python package and environment manager. We provide a one-shot setup script.

Prereqs: Linux/macOS; Python 3.11+ recommended (3.12 preferred).

1. Clone and enter the repo
```bash
git clone <your-fork-or-repo-url>
cd METAINFORMANT
```

2. Run setup (installs uv if missing, creates venv, installs deps)
```bash
bash scripts/setup_uv.sh
```

3. Activate the environment (subsequent shells)
```bash
source .venv/bin/activate
```

4. Sanity check
```bash
uv run python -V
pytest -q
```

Notes
- All dependencies are declared in `pyproject.toml`. Use `uv add <pkg>` to add more.
- Pr

### License
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at `https://www.apache.org/licenses/LICENSE-2.0`.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the `LICENSE` file for the specific language governing permissions and limitations under the License.


