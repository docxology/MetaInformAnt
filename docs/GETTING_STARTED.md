# Getting Started with METAINFORMANT

**New to METAINFORMANT?** This guide will walk you through your first analysis in under 5 minutes, from installation to visualization. For detailed setup, see [SETUP.md](SETUP.md).

---

## Overview

METAINFORMANT is a comprehensive bioinformatics toolkit covering 30 modules across genomics, transcriptomics, proteomics, and systems biology. This quick start demonstrates a typical workflow: generate synthetic DNA data, analyze composition, and visualize the results.

---

## Before You Begin

Ensure you have:
- **Python 3.11+** installed
- **UV** package manager (`pip install uv` or see [UV_SETUP.md](UV_SETUP.md))
- Repository cloned: `git clone https://github.com/docxology/metainformant.git`
- Dependencies installed: `bash scripts/package/setup.sh`

Verify: `uv run python -V` (should be 3.11 or 3.12).

---

## 1. Choose the Right Module

Not sure which module to use? Refer to the **Module Selection Decision Tree** in the [documentation index](../docs/index.md#module-selection-decision-tree). It helps you select based on data type and analysis goal.

For DNA sequence analysis, we'll use the `dna` module. For visualization, we'll use `visualization.plots.basic`.

---

## 2. Create Sample Data

METAINFORMANT encourages using **synthetic data** for learning and testing. We'll create sample DNA sequences and save them using `metainformant.core.io`.

```python
# File: get_started_example.py
from metainformant.core import io
import numpy as np

# Create synthetic DNA sequences of varying lengths
sequences = {
    "gene_alpha": "ATCGATCGATCGATCGATCGATCG",  # 24 bp
    "gene_beta":  "GCTAGCTAGCTAGCTAGCTA",        # 20 bp
    "gene_gamma": "ATATATATATATATATATATATAT",    # 24 bp
    "gene_delta": "CGCGCGCGCGCGCGCGCGCG",        # 20 bp
}

# Save to FASTA format in the output directory
output_dir = io.ensure_directory("output/getting_started")
fasta_path = output_dir / "sample_sequences.fasta"

with open(fasta_path, "w") as f:
    for name, seq in sequences.items():
        f.write(f">{name}\n{seq}\n")

print(f"Saved {len(sequences)} sequences to {fasta_path}")
```

Run it:
```bash
uv run python get_started_example.py
```

---

## 3. Basic Analysis with the `dna` Module

Now analyze the DNA sequences: compute GC content, sequence lengths, and motif frequencies.

```python
from metainformant.dna.sequence import composition, core, motifs
from pathlib import Path

# Read the FASTA file we just created (returns {sequence_id: sequence})
fasta_path = Path("output/getting_started/sample_sequences.fasta")
sequences = core.read_fasta(fasta_path)

# Compute basic statistics
print("Sequence Analysis Results:")
print("-" * 40)
for name, seq in sequences.items():
    length = len(seq)
    gc = composition.gc_content(seq)
    print(f"{name}:")
    print(f"  Length: {length} bp")
    print(f"  GC content: {gc:.2%}")
    print(f"  AT content: {1 - gc:.2%}")  # sequences contain only A, C, G, T

# Find a specific motif across all sequences
motif = "ATCG"
print(f"\nMotif search ('{motif}'):")
for name, seq in sequences.items():
    positions = motifs.find_motif_positions(seq, motif)
    print(f"  {name}: found at positions {positions}")
```

**Note:** The example uses `metainformant.dna.sequence.core`, `composition` and `motifs`; consult the [DNA module documentation](dna/sequences.md) for the other sequence functions.

---

## 4. Visualize with `visualization.plots.basic`

Plot GC content across sequences using the basic plotting module.

```python
import matplotlib.pyplot as plt
from metainformant.visualization.plots import basic
import numpy as np

# Prepare data
names = list(sequences.keys())
gc_contents = [composition.gc_content(seq) for seq in sequences.values()]

# Create a bar plot
fig, ax = plt.subplots(figsize=(10, 6))
basic.bar_plot(
    np.arange(len(names)),
    np.array(gc_contents),
    ax=ax,
    color="steelblue",
    edgecolor="black",
)
ax.set_xlabel("Gene")
ax.set_ylabel("GC Content")
ax.set_title("GC Content by Gene")
ax.set_xticks(np.arange(len(names)))
ax.set_xticklabels(names, rotation=45)
plt.tight_layout()

# Save figure
output_plot = io.ensure_directory("output/getting_started") / "gc_content.png"
plt.savefig(output_plot, dpi=150, bbox_inches="tight")
print(f"Plot saved to {output_plot}")
plt.show()  # Optional: display interactively
```

---

## 5. Add Caching to Expensive Computations

If your analysis involves heavy computation (e.g., multiple sequence alignment), use the JSON-based caching system.

```python
from metainformant.core.io import cache
from pathlib import Path

# Initialize a cache (stores to output/cache by default)
cache_dir = Path("output/getting_started/cache")
my_cache = cache.JsonCache(cache_dir, ttl_seconds=3600)  # 1-hour TTL

# Example expensive computation: computing a distance matrix
def compute_distance_matrix(seq1: str, seq2: str) -> float:
    """Simulate expensive pairwise alignment distance."""
    # In reality, this could call dna.alignment.distances.jukes_cantor_distance()
    return 1.0 - composition.gc_content(seq1) * composition.gc_content(seq2)

# Cache results automatically
seq_list = list(sequences.values())
distance_matrix = np.zeros((len(seq_list), len(seq_list)))
for i, s1 in enumerate(seq_list):
    for j, s2 in enumerate(seq_list):
        key = f"dist_{i}_{j}"
        cached = my_cache.get(key)
        if cached is None:
            dist = compute_distance_matrix(s1, s2)
            my_cache.set(key, dist)
        else:
            dist = cached
        distance_matrix[i, j] = dist

print("Distance matrix computed (with caching):")
print(distance_matrix)
```

The cache persists between runs, avoiding recomputation.

---

## 6. Configure Logging

METAINFORMANT provides structured logging via `metainformant.core.utils.logging`. Use it in scripts.

```python-snippet
from metainformant.core.utils.logging import setup_logger

# Setup: can be done once at program start
logger = setup_logger(__name__, level="INFO")  # Options: DEBUG, INFO, WARNING, ERROR

logger.info("Starting GETTING_STARTED workflow")
logger.debug(f"Loaded {len(sequences)} sequences")
logger.warning("This is just a demo; real analyses require more data")
logger.error("Errors are logged with full context")
```

By default, logs go to stderr as `timestamp | level | logger name | message`. To log to file:
```python
logger = setup_logger(__name__, log_file="output/getting_started/workflow.log", level="DEBUG")
```

---

## 7. Full Pipeline Script

Putting it all together (saved as `scripts/getting_started_pipeline.py`):

```python-snippet
#!/usr/bin/env python
"""METAINFORMANT Getting Started Pipeline

A minimal end-to-end example: create synthetic DNA data, analyze composition,
generate a plot, and demonstrate caching + logging.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from metainformant.core import io
from metainformant.core.io import cache
from metainformant.core.utils.logging import setup_logger
from metainformant.dna.sequence import composition
from metainformant.visualization.plots import basic

def main() -> None:
    # Configure logging
    logger = setup_logger(__name__, level="INFO")
    logger.info("=== METAINFORMANT Getting Started ===")

    # ── Step 1: Prepare output directories ──────────────────────────────────
    output_dir = io.ensure_directory("output/getting_started")
    cache_dir = io.ensure_directory(output_dir / "cache")

    # ── Step 2: Synthetic data ─────────────────────────────────────────────
    sequences = {
        "gene_alpha": "ATCGATCGATCGATCGATCGATCG",
        "gene_beta":  "GCTAGCTAGCTAGCTAGCTA",
        "gene_gamma": "ATATATATATATATATATATATAT",
        "gene_delta": "CGCGCGCGCGCGCGCGCGCG",
    }
    logger.info(f"Created {len(sequences)} synthetic sequences")

    # Save as FASTA
    fasta_path = output_dir / "sample.fasta"
    with open(fasta_path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")
    logger.info(f"FASTA saved → {fasta_path}")

    # ── Step 3: Analysis ───────────────────────────────────────────────────
    names = list(sequences.keys())
    gc_contents = [composition.gc_content(seq) for seq in sequences.values()]
    lengths     = [len(seq) for seq in sequences.values()]

    for name, gc, length in zip(names, gc_contents, lengths):
        logger.info(f"{name}: gc={gc:.2%}, len={length}")

    # ── Step 4: Caching expensive pairwise distances ───────────────────────
    dist_cache = cache.JsonCache(cache_dir, ttl_seconds=1800)
    seq_list = list(sequences.values())
    n = len(seq_list)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            key = f"dist_{i}_{j}"
            dist = dist_cache.get(key)
            if dist is None:
                # Simulate heavy alignment
                dist = 1.0 - gc_contents[i] * gc_contents[j]
                dist_cache.set(key, dist)
            dist_matrix[i, j] = dist
    logger.info("Pairwise distance matrix computed (cached)")

    # ── Step 5: Plot GC content ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 5))
    basic.bar_plot(
        np.arange(len(names)),
        np.array(gc_contents),
        ax=ax,
        color="teal",
        edgecolor="black",
    )
    ax.set_xlabel("Gene")
    ax.set_ylabel("GC Content")
    ax.set_title("GC Content Across Synthetic Genes")
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right")
    plt.tight_layout()

    plot_path = output_dir / "gc_content.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    logger.info(f"Plot saved → {plot_path}")
    plt.close(fig)

    logger.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
```

Run it:
```bash
uv run python scripts/getting_started_pipeline.py
```
Check the `output/getting_started/` directory for results.

---

## Next Steps

You've completed a basic METAINFORMANT workflow. Continue exploring:

- **Tutorials** – Step-by-step guides for specific domains: [docs/TUTORIALS.md](TUTORIALS.md)
- **Task Reference** – Quick how-to's for common operations: [docs/tasks/](tasks/)
  - [analyze_dna.md](tasks/analyze_dna.md)
  - [visualize_results.md](tasks/visualize_results.md)
  - [run_gwas.md](tasks/run_gwas.md)
  - and more...
- **Module Documentation** – In-depth guides:
  - [DNA analysis](dna/index.md)
  - [RNA-seq & amalgkit](rna/index.md)
  - [GWAS pipeline](gwas/index.md)
  - [Visualization library](visualization/index.md)
- **Architecture** – System design and data flow: [docs/architecture.md](architecture.md)
- **CLI Reference** – All available commands: [docs/cli.md](cli.md)
- **Testing** – Writing tests and contributing quality code: [docs/testing.md](testing.md)

---

## Common Questions

**Q: How do I select which module to use?**
A: See the [Module Selection Decision Tree](../docs/index.md#module-selection-decision-tree) in the documentation index. It maps your data type and goal to the appropriate module.

**Q: Where should I put my data files?**
A: Use the `data/` directory for input datasets (gitignored). Outputs go to `output/`. See [SETUP.md](SETUP.md) for directory policies.

**Q: How do I enable debug logging?**
A: Call `logger = setup_logger(__name__, level="DEBUG")` at the start of your script, or pass `--log-level DEBUG` to the CLI if available.

**Q: My script is slow—how can I speed it up?**
A: Use `metainformant.core.execution.parallel` for parallel processing, enable caching with `metainformant.core.io.cache`, and profile with `pytest-benchmark`. See [docs/tasks/performance_tuning.md](../docs/tasks/performance_tuning.md).

**Q: Where do I find example code?**
A: Check module READMEs in `src/metainformant/<module>/README.md`, tutorials in [docs/TUTORIALS.md](TUTORIALS.md), and the `examples/` directory if present.

---

## Need Help?

- **Documentation issues?** Edit any `.md` file and open a PR — documentation is collaborative.
- **Code questions?** Open a [GitHub Discussion](https://github.com/docxology/metainformant/discussions).
- **Chat:** Join `#metainformant:matrix.org` for real-time help.
- **Bug reports:** Use [GitHub Issues](https://github.com/docxology/metainformant/issues) with the bug template.

Happy exploring!
