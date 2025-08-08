# METAINFORMANT: Integrated Multi-Omic Biological Systems Modeling

METAINFORMANT is a unified, modular toolkit for integrated multi-omic analysis across domains:
- DNA: genome retrieval, metadata, storage, plus alignment, MSA, phylogeny, and population genetics
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
        genomes.py                # Accession validation
        entrez.py                 # Entrez fetch helpers
        ncbi.py                   # NCBI Datasets wrappers
        sequences.py              # FASTA I/O
        alignment.py              # Pairwise alignment
        msa.py                    # Lightweight MSA (progressive)
        phylogeny.py              # NJ tree + Newick export
        population.py             # Popgen stats (alleles, pi, Tajima's D, Fst)
      rna/
        __init__.py
         pipeline.py               # Utility helpers
         amalgkit.py               # Thin modular wrapper around amalgkit CLI
         workflow.py               # Plan/execute complete amalgkit workflows
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

### DNA quickstart

```python
from metainformant.dna import sequences, alignment, phylogeny, population

seqs = sequences.read_fasta("tests/data/dna/toy.fasta")
aln = alignment.global_align(seqs["A"], seqs["B"])  # .score, .aligned_seq1, .aligned_seq2
tree = phylogeny.neighbor_joining_tree(seqs)
print(phylogeny.to_newick(tree))

pi = population.nucleotide_diversity(["AAAA", "AAAT"])  # 0.25
```

### RNA: Transcriptomic meta-analysis via amalgkit

Install amalgkit (external dependency):

```bash
pip install git+https://github.com/kfuku52/amalgkit
amalgkit -h  # verify
```

Programmatic usage in Python:

```python
from pathlib import Path
from metainformant.rna import (
    check_cli_available,
    metadata, integrate, config, select, getfastq, quant, merge, cstmm, curate, csca, sanity,
    AmalgkitWorkflowConfig, plan_workflow, execute_workflow,
)

ok, help_text = check_cli_available()
assert ok, help_text

# Run a single step
result = metadata({"threads": 4})
print(result.returncode, result.stdout[:200])

# Plan and execute a full workflow
cfg = AmalgkitWorkflowConfig(work_dir=Path("/path/to/work"), threads=8, species_list=["Apis_mellifera"])
steps = plan_workflow(cfg)
codes = execute_workflow(cfg)
```

CLI usage:

```bash
python -m metainformant rna plan --work-dir /tmp/amg --threads 4 --species Apis_mellifera
python -m metainformant rna run --work-dir /tmp/amg --threads 4 --species Apis_mellifera --check
```

### License
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at `https://www.apache.org/licenses/LICENSE-2.0`.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the `LICENSE` file for the specific language governing permissions and limitations under the License.


