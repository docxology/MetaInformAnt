from __future__ import annotations

import random
from typing import Dict, Literal, Optional

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def neighbor_joining_tree(id_to_seq: Dict[str, str]):
    records = [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in id_to_seq.items()]
    alignment = MultipleSeqAlignment(records)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def upgma_tree(id_to_seq: Dict[str, str]):
    records = [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in id_to_seq.items()]
    alignment = MultipleSeqAlignment(records)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return tree


def to_newick(tree) -> str:
    from io import StringIO

    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue()


def upgma_tree(id_to_seq: Dict[str, str]):
    records = [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in id_to_seq.items()]
    alignment = MultipleSeqAlignment(records)
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    return constructor.upgma(dm)


def bootstrap_support(
    id_to_seq: Dict[str, str], n_replicates: int = 100, method: str = "nj", random_state: int | None = None
) -> Dict[str, float]:
    """Simple bootstrap with optional deterministic seed.

    Returns mapping of clade name (sorted, semicolon-joined) to support fraction.
    """
    seq_ids = list(id_to_seq.keys())
    if not seq_ids:
        return {}
    L = min(len(s) for s in id_to_seq.values())
    if L == 0:
        return {}

    rng = random.Random(random_state)

    def build_tree(resampled: Dict[str, str]):
        if method == "nj":
            return neighbor_joining_tree(resampled)
        return upgma_tree(resampled)

    def splits(tree):
        # Represent clades by sorted tuple of terminal names
        return [tuple(sorted([t.name for t in clade.get_terminals()])) for clade in tree.get_nonterminals()]

    counts: Dict[tuple[str, ...], int] = {}
    for _ in range(n_replicates):
        idxs = [rng.randrange(L) for _ in range(L)]
        resampled = {k: "".join(v[i] for i in idxs) for k, v in id_to_seq.items()}
        tree = build_tree(resampled)
        for s in splits(tree):
            counts[s] = counts.get(s, 0) + 1
    return {";".join(s): c / n_replicates for s, c in counts.items()}


def to_ascii(tree) -> str:
    """Return an ASCII art representation of the tree."""
    from io import StringIO

    handle = StringIO()
    Phylo.draw_ascii(tree, file=handle)
    return handle.getvalue()


def basic_tree_stats(tree) -> Dict[str, int]:
    """Return simple stats for a Phylo tree: number of terminals and clades."""
    return {
        "num_terminals": len(tree.get_terminals()),
        "num_clades": sum(1 for _ in tree.find_clades()),
    }


def bootstrap_support(
    id_to_seq: Dict[str, str],
    *,
    n_replicates: int = 100,
    method: Literal["nj", "upgma"] = "nj",
    random_state: Optional[int] = None,
) -> Dict[frozenset[str], float]:
    """Bootstrap clade support using simple column resampling.

    Returns a mapping from clade (as frozenset of leaf names) to support in [0,1].
    """
    if not id_to_seq:
        return {}
    ids = list(id_to_seq.keys())
    L = min(len(seq) for seq in id_to_seq.values())
    if L == 0 or n_replicates <= 0:
        return {}

    def build_tree(data: Dict[str, str]):
        return neighbor_joining_tree(data) if method == "nj" else upgma_tree(data)

    def clade_splits(tree) -> set[frozenset[str]]:
        leaves = {t.name for t in tree.get_terminals()}
        splits: set[frozenset[str]] = set()
        for clade in tree.get_nonterminals():
            terminals = {t.name for t in clade.get_terminals()}
            if 0 < len(terminals) < len(leaves):
                splits.add(frozenset(terminals))
        return splits

    # Original tree splits
    base_tree = build_tree(id_to_seq)
    base_splits = clade_splits(base_tree)
    counts: dict[frozenset[str], int] = {s: 0 for s in base_splits}

    import random

    rng = random.Random(random_state)

    for _ in range(n_replicates):
        positions = [rng.randrange(L) for _ in range(L)]
        boot_data: Dict[str, str] = {i: "".join(id_to_seq[i][pos] for pos in positions) for i in ids}
        tree = build_tree(boot_data)
        for s in clade_splits(tree):
            if s in counts:
                counts[s] += 1

    return {s: counts[s] / n_replicates for s in counts}


def nj_tree_from_kmer(id_to_seq: Dict[str, str], *, k: int = 3, metric: str = "cosine"):
    """Build NJ tree from k-mer distances via `distances.kmer_distance_matrix`."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    from . import distances as distances_mod

    ids = list(id_to_seq.keys())
    dm_values = distances_mod.kmer_distance_matrix(id_to_seq, k=k, metric=metric)
    # Biopython DistanceMatrix expects lower-triangular distances without diagonal
    matrix: list[list[float]] = []
    for i in range(len(ids)):
        row = []
        for j in range(i):
            row.append(dm_values[i][j])
        row.append(0.0)  # diagonal element
        matrix.append(row)
    dm = DistanceMatrix(names=ids, matrix=matrix)
    constructor = DistanceTreeConstructor()
    return constructor.nj(dm)
