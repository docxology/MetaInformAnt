# METAINFORMANT: Integrated Multi-Omic Biological Systems Modeling

METAINFORMANT is a unified, modular toolkit for integrated multi-omic analysis across domains:
- DNA: genome retrieval, metadata, storage, plus alignment, MSA, phylogeny, and population genetics
- RNA: transcriptomic metadata, download, quantification, and curation
- Protein: proteome retrieval and descriptive analytics
- Epigenome: DNA modifications (e.g., chromatin/methylation) and other epigenetics pipelines
- Ontology: functional annotation (e.g., GO) and taxonomic utilities
- Phenotype: morphological and behavioral phenotype curation and web scraping (e.g., AntWiki)
- Ecology: ecology metadata, download, curation
- Visualization: cohesive plotting and animation utilities for figures and trees

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
      __main__.py                 # `uv run python -m metainformant` CLI entry
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
      visualization/
        __init__.py                 # Unified plotting/animation API
        plots.py                    # Line plots, heatmaps, pairplots (matplotlib/seaborn)
        animations.py               # Simple time-series animations
        trees.py                    # Biopython Phylo tree -> matplotlib
      math/                       # Theoretical & quantitative biology
        __init__.py               # Price eq., kin/multilevel selection, DDM
      simulation/                 # Synthetic data and agent-based simulators
        __init__.py               # DNA/RNA/protein generators, grid agents
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

2. Run setup (uses uv; creates venv, installs deps)
```bash
bash scripts/setup_uv.sh
```

3. Activate the environment (subsequent shells) or use uv run directly
```bash
source .venv/bin/activate
```

4. Sanity check
```bash
uv run python -V
uv run pytest -q
```

Notes
- All dependencies are declared in `pyproject.toml`. Prefer `uv add <pkg>` for adding new dependencies, or `uv pip install <pkg>` for one-off installations.
- Prefer `uv run <cmd>` over activating the environment when scripting/CI.

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
# Install amalgkit (external dependency):
uv pip install git+https://github.com/kfuku52/amalgkit
amalgkit -h  # verify
```

### Visualization quickstart

```python
import matplotlib.pyplot as plt
from metainformant.visualization import lineplot, heatmap, animate_time_series

# Line plot
ax = lineplot(None, [0.1, 0.4, 0.2, 0.6], label="signal")
ax.set_xlabel("index"); ax.set_ylabel("value")
plt.close(ax.figure)

# Heatmap
ax = heatmap([[1, 0], [0, 1]])
plt.close(ax.figure)

# Animation (save as GIF/MP4 using matplotlib writers installed in your env)
fig, anim = animate_time_series([0, 1, 0, 1], interval_ms=100)
# Example save (requires imagemagick/ffmpeg):
# anim.save("timeseries.gif", writer="imagemagick", fps=10)
plt.close(fig)
```

## Math: Theoretical and quantitative biology

- Price equation decomposition: `metainformant.math.price_equation`
- Kin selection (Hamilton's rule): `metainformant.math.kin_selection_response`
- Multilevel selection partition: `metainformant.math.multilevel_selection_decomposition`
- Drift–Diffusion Model (DDM): `metainformant.math.ddm_analytic_accuracy`, `metainformant.math.ddm_mean_decision_time`

Example:

```python
from metainformant.math import price_equation, kin_selection_response

cov, trans, total = price_equation([1.0, 1.2, 0.9], [0.2, 0.4, 0.1], [0.25, 0.35, 0.15])
resp = kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.1)
```

## Simulation: Synthetic data and agents

- DNA/protein generators: `generate_random_dna`, `generate_random_protein`, `mutate_sequence`
- RNA counts (NB): `simulate_counts_negative_binomial`
- Minimal agent-based grid: `GridWorld`, `Agent`

Example:

```python
from metainformant.simulation import generate_random_dna, simulate_counts_negative_binomial, GridWorld

seq = generate_random_dna(1000)
counts = simulate_counts_negative_binomial(100, 6)
world = GridWorld(10, 10, num_agents=5)
world.step()
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
uv run python -m metainformant rna plan --work-dir /tmp/amg --threads 4 --species Apis_mellifera
uv run python -m metainformant rna run --work-dir /tmp/amg --threads 4 --species Apis_mellifera --check
```

### License
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at `https://www.apache.org/licenses/LICENSE-2.0`.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the `LICENSE` file for the specific language governing permissions and limitations under the License.
