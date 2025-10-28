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
        cache.py                  # JSON-based caching with TTL
        config.py                 # Configuration management
        db.py                     # Database integration (optional)
        hash.py                   # Content hashing utilities
        io.py                     # File I/O utilities
        logging.py                # Logging configuration
        parallel.py               # Parallel processing
        paths.py                  # Path handling utilities
        text.py                   # Text processing utilities
      dna/                        # DNA sequence analysis
        __init__.py
        alignment.py              # Sequence alignment algorithms
        codon.py                  # Codon usage analysis
        composition.py            # Sequence composition analysis
        consensus.py              # Consensus sequence generation
        distances.py              # Evolutionary distance calculations
        entrez.py                 # NCBI Entrez database integration
        fastq.py                  # FASTQ processing and quality
        genomes.py                # Genome data handling
        motifs.py                 # Motif discovery
        msa.py                    # Multiple sequence alignment
        mutations.py              # Mutation analysis
        ncbi.py                   # NCBI database integration
        phylogeny.py              # Phylogenetic analysis
        population.py             # Population genetics
        restriction.py            # Restriction enzyme analysis
        sequences.py              # DNA sequence I/O
        transcription.py          # DNA to RNA transcription
        translation.py            # Genetic code translation
        variants.py               # Variant analysis
      rna/                        # RNA analysis and workflows
        __init__.py
        amalgkit.py               # Amalgkit CLI wrapper
        configs.py                # Workflow configuration
        deps.py                   # Dependency management
        pipeline.py               # RNA analysis pipeline
        steps/                    # Individual workflow steps
        workflow.py               # Complete workflow orchestration
      protein/                    # Protein analysis
        __init__.py
        alignment.py              # Protein sequence alignment
        alphafold.py              # AlphaFold integration
        contacts.py               # Protein contact analysis
        interpro.py               # InterPro domain analysis
        pdb.py                    # PDB structure handling
        proteomes.py              # Proteome analysis
        secondary.py              # Secondary structure prediction
        sequences.py              # Protein sequence I/O
        structure.py              # Protein structure analysis
        structure_io.py           # Structure file I/O
        uniprot.py                # UniProt database integration
      math/                       # Mathematical biology
        __init__.py
        coalescent.py             # Coalescent theory
        ddm.py                    # Drift-diffusion models
        dynamics.py               # Population dynamics
        effective_size.py         # Effective population size
        epidemiology.py           # Disease modeling
        egt.py                    # Evolutionary game theory
        fst.py                    # F-statistics
        ld.py                     # Linkage disequilibrium
        popgen.py                 # Population genetics
        price.py                  # Price equation
        quantgen.py               # Quantitative genetics
        selection.py              # Natural selection
        selection_experiments/    # Selection experiment simulations
          __init__.py
          cli.py                  # Command-line interface
          model.py                # Selection models
          plotting.py             # Result visualization
          README.md               # Experiment documentation
      ml/                         # Machine learning
        __init__.py
        classification.py         # Supervised classification
        dimensionality.py         # Dimensionality reduction
        features.py               # Feature selection
        regression.py             # Regression models
        validation.py             # Model validation
      multiomics/                 # Multi-omics integration
        __init__.py
        integration.py            # Cross-omic data integration
      networks/                   # Network analysis
        __init__.py
        community.py              # Community detection
        graph.py                  # Graph algorithms
        pathway.py                # Pathway analysis
        ppi.py                    # Protein-protein interactions
        regulatory.py             # Regulatory networks
      ontology/                   # Ontology analysis
        __init__.py
        go.py                     # Gene Ontology
        obo.py                    # OBO format parsing
        query.py                  # Ontology querying
        types.py                  # Ontology data types
      phenotype/                  # Phenotype analysis
        __init__.py
        antwiki.py                # AntWiki integration
      ecology/                    # Ecology analysis
        __init__.py
        community.py              # Community ecology analysis
      epigenome/                  # Epigenomics
        __init__.py
        methylation.py            # DNA methylation analysis
        tracks.py                 # Genomic track processing
      quality/                    # Quality assessment
        __init__.py
        contamination.py          # Sequence contamination detection
        fastq.py                  # FASTQ quality analysis
        metrics.py                # Quality metrics and scoring
      simulation/                 # Data simulation
        __init__.py
        agents.py                 # Agent-based modeling
        rna.py                    # RNA count simulation
        sequences.py              # Sequence generation
      singlecell/                 # Single-cell analysis
        __init__.py
        clustering.py             # Single-cell clustering
        config/                   # Configuration files
        dimensionality.py         # Dimensionality reduction
        integration.py            # Data integration
        preprocessing.py          # Preprocessing pipelines
        trajectory.py             # Trajectory analysis
        visualization.py          # Single-cell visualization
      visualization/              # Data visualization
        __init__.py
        animations.py             # Animation utilities
        plots.py                  # Plotting functions
        trees.py                  # Phylogenetic tree plotting
  tests/                          # Comprehensive test suite
    test_*.py                     # Domain-specific tests
    data/                         # Test data fixtures
  config/                         # Configuration files
    amalgkit_*.yaml              # RNA workflow configurations
  scripts/                        # Utility scripts
    *.sh                         # Shell scripts for setup/testing
  docs/                           # Documentation
    */README.md                   # Module documentation
  output/                         # Output directory (ephemeral)
  pyproject.toml                  # Project configuration
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
