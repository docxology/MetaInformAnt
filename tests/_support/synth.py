"""Shared synthetic-data factories for METAINFORMANT tests.

These builders replace the copy-pasted ``_make_*`` helpers that were
duplicated across ``tests/*/`` test files. Every factory is deterministic
(seed-controlled), produces small in-memory real data (no test doubles), and is
typed. Import them as::

    from tests._support.synth import (
        make_expression_frame,
        make_ortholog_frame,
        make_phenotype_map,
        make_genotype_map,
        make_sample_metadata_tsv,
        make_species_tree,
    )

REAL-IMPLEMENTATION policy: factories build real numpy/pandas objects and
real files under ``tmp_path``; they never fake the unit under test.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Expression data
# ---------------------------------------------------------------------------


def make_expression_frame(
    genes: int = 5,
    samples: int = 6,
    seed: int = 42,
    mean: float = 10.0,
    sd: float = 3.0,
    gene_prefix: str = "gene_",
    sample_prefix: str = "sample_",
) -> pd.DataFrame:
    """Genes x samples expression DataFrame (genes as index)."""
    rng = np.random.RandomState(seed)
    data = rng.normal(mean, sd, size=(genes, samples))
    index = [f"{gene_prefix}{i}" for i in range(genes)]
    columns = [f"{sample_prefix}{j}" for j in range(samples)]
    return pd.DataFrame(data, index=index, columns=columns)


def make_expression_dict(
    samples: int = 20,
    genes: int = 5,
    seed: int = 42,
    mean: float = 10.0,
    sd: float = 3.0,
) -> Dict[str, Dict[str, float]]:
    """Nested {gene_id: {sample_id: value}} mapping (random module based)."""
    rnd = random.Random(seed)
    sample_ids = [f"sample_{i}" for i in range(samples)]
    return {f"gene_{g}": {s: rnd.gauss(mean, sd) for s in sample_ids} for g in range(genes)}


def make_de_matrix(
    n_cells: int = 80,
    n_genes: int = 100,
    n_de_genes: int = 10,
    seed: int = 42,
) -> Tuple[List[List[float]], List[int], List[str]]:
    """Expression matrix with known differentially expressed genes.

    Group 0 (first half of cells) is baseline; group 1 (second half) has
    elevated expression in the first ``n_de_genes`` genes.
    """
    rng = np.random.RandomState(seed)
    half = n_cells // 2
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    matrix = rng.exponential(1.0, size=(n_cells, n_genes))
    matrix[half:, :n_de_genes] += rng.exponential(5.0, size=(half, n_de_genes))
    groups = [0] * half + [1] * half
    return matrix.tolist(), groups, gene_names


# ---------------------------------------------------------------------------
# Orthology / cross-species
# ---------------------------------------------------------------------------


def make_ortholog_frame(
    pairs: Sequence[Tuple[str, str]],
    source_col: str = "human",
    target_col: str = "mouse",
) -> pd.DataFrame:
    """Ortholog DataFrame from a list of (source, target) gene pairs."""
    return pd.DataFrame({source_col: [p[0] for p in pairs], target_col: [p[1] for p in pairs]})


def make_species_tree(leaves: Sequence[str] = ("human", "mouse"), clade_name: str = "clade_1") -> Dict[str, Any]:
    """Nested-dict phylogenetic tree.

    Two leaves produce a flat root; three or more produce one internal
    clade (named ``clade_name``) containing the first two leaves plus
    remaining leaf children.
    """
    if len(leaves) < 2:
        raise ValueError("make_species_tree needs at least two leaves")
    first = {"name": leaves[0], "distance": 0.1}
    second = {"name": leaves[1], "distance": 0.15}
    if len(leaves) == 2:
        children = [first, second]
    else:
        children = [
            {"name": clade_name, "distance": 0.3, "children": [first, second]},
        ] + [{"name": leaf, "distance": 0.5} for leaf in leaves[2:]]
    return {"name": "root", "distance": 0.0, "children": children}


# ---------------------------------------------------------------------------
# Phenotype / genotype
# ---------------------------------------------------------------------------


def make_phenotype_map(
    n: int = 20,
    seed: Optional[int] = None,
    slope: float = 2.0,
    intercept: float = 5.0,
) -> Dict[str, float]:
    """Deterministic {sample_id: quantitative phenotype} mapping."""
    if seed is not None:
        rnd = random.Random(seed)
        return {f"sample_{i}": rnd.gauss(intercept + slope * i, 1.0) for i in range(n)}
    return {f"sample_{i}": float(i * slope + intercept) for i in range(n)}


def make_genotype_map(
    n: int = 20,
    n_variants: int = 3,
    seed: int = 42,
    prefix: str = "rs",
) -> Dict[str, List[int]]:
    """{variant_id: [dosage per sample]} with values in {0, 1, 2}."""
    rnd = random.Random(seed)
    return {f"{prefix}{v + 1000}": [rnd.choice([0, 1, 2]) for _ in range(n)] for v in range(n_variants)}


# ---------------------------------------------------------------------------
# File-producing factories (real files under caller-provided directories)
# ---------------------------------------------------------------------------


def make_sample_metadata_tsv(
    directory: Path,
    rows: int = 3,
    filename: str = "metadata.tsv",
    extra_columns: Optional[Mapping[str, Sequence[Any]]] = None,
) -> Path:
    """Write a small AMALGKIT-style sample metadata TSV; returns its path."""
    header = [
        "sample_id",
        "species_name",
        "data_type",
        "tissue",
        "filename",
        "SRR_accession",
    ]
    lines = ["\t".join(header)]
    for i in range(rows):
        lines.append(
            "\t".join(
                [
                    f"sample_{i}",
                    f"species_{i}",
                    "rna_seq",
                    "whole_body",
                    f"sample_{i}.fastq.gz",
                    f"SRR{1000000 + i}",
                ]
            )
        )
    if extra_columns:
        raise ValueError("extra_columns not supported in minimal factory; extend as needed")
    path = Path(directory) / filename
    path.write_text("\n".join(lines) + "\n")
    return path


def make_fasta(
    directory: Path,
    sequences: Sequence[Tuple[str, str]],
    filename: str = "sequences.fasta",
    line_width: int = 60,
) -> Path:
    """Write a real FASTA file from (header, sequence) pairs; returns path."""
    chunks = []
    for header, seq in sequences:
        chunks.append(f">{header}")
        for start in range(0, len(seq), line_width):
            chunks.append(seq[start : start + line_width])
    path = Path(directory) / filename
    path.write_text("\n".join(chunks) + "\n")
    return path


def make_assoc_results(
    n: int = 100,
    seed: int = 42,
    n_chroms: int = 22,
) -> List[Dict[str, Any]]:
    """Synthetic GWAS association results across several chromosomes."""
    rng = np.random.default_rng(seed)
    results: List[Dict[str, Any]] = []
    for i in range(n):
        results.append(
            {
                "variant_id": f"rs{100000 + i}",
                "chrom": str((i % n_chroms) + 1),
                "pos": int(rng.integers(1_000_000, 250_000_000)),
                "p_value": float(rng.uniform(1e-10, 1.0)),
                "beta": float(rng.normal(0, 0.5)),
                "se": float(abs(rng.normal(0.1, 0.05))),
            }
        )
    return results


def make_spatial_cells(
    n: int = 50,
    n_types: int = 3,
    seed: int = 42,
    side: float = 100.0,
) -> Tuple[List[str], np.ndarray]:
    """(cell-type labels, (n, 2) coordinates) for spatial analyses."""
    rng = np.random.RandomState(seed)
    coords = rng.rand(n, 2) * side
    type_names = [f"type{chr(ord('A') + t)}" for t in range(max(n_types, 1))]
    types = rng.choice(type_names, size=n).tolist()
    return types, coords


def make_config_yaml(
    directory: Path,
    config: Mapping[str, Any],
    filename: str = "config.yaml",
) -> Path:
    """Write a config YAML mapping; returns path. Uses yaml if available."""
    import yaml  # deferred: repo dependency

    path = Path(directory) / filename
    path.write_text(yaml.safe_dump(dict(config), sort_keys=False))
    return path
