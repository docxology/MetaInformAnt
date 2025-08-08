from __future__ import annotations

from typing import Dict

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
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


def to_newick(tree) -> str:
    from io import StringIO

    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue()


