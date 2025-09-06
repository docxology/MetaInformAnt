from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

from Bio.Seq import Seq


def translate_dna(seq: str, *, to_stop: bool = False, table: int | str = 1) -> str:
    return str(Seq(seq).translate(table=table, to_stop=to_stop))


@dataclass
class ORF:
    start: int
    end: int
    frame: int
    strand: int
    protein: str
    start_codon: str
    stop_codon: str | None


def find_orfs(seq: str, *, min_aa: int = 50, include_reverse: bool = True) -> List[ORF]:
    strands = [(seq, +1)]
    if include_reverse:
        from metainformant.dna.sequences import reverse_complement

        strands.append((reverse_complement(seq), -1))

    orfs: list[ORF] = []
    for s, strand in strands:
        n = len(s)
        for frame in range(3):
            i = frame
            while i + 3 <= n:
                codon = s[i : i + 3]
                if codon == "ATG":
                    # scan to stop
                    j = i + 3
                    while j + 3 <= n:
                        stop = s[j : j + 3]
                        if stop in {"TAA", "TAG", "TGA"}:
                            prot = translate_dna(s[i:j])
                            if len(prot) >= min_aa:
                                start_pos = i if strand == 1 else n - (i + 3)
                                end_pos = j + 3 if strand == 1 else n - j
                                orfs.append(
                                    ORF(
                                        start=start_pos,
                                        end=end_pos,
                                        frame=frame,
                                        strand=strand,
                                        protein=prot,
                                        start_codon=codon,
                                        stop_codon=stop,
                                    )
                                )
                            break
                        j += 3
                    i = j
                else:
                    i += 3
    return orfs
