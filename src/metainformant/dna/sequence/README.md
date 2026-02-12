# Sequence

Core DNA sequence operations: FASTA I/O, composition analysis, consensus calling, k-mer profiling, motif discovery, and restriction enzyme mapping.

## Contents

| File | Purpose |
|------|---------|
| `core.py` | FASTA read/write, reverse complement, GC content, ORFs, sequence validation |
| `composition.py` | GC/AT skew, melting temperature, nucleotide and dinucleotide frequencies |
| `consensus.py` | Consensus sequence generation with IUPAC ambiguity and bootstrap support |
| `kmer.py` | K-mer counting, spectrum analysis, microsatellite and homopolymer detection |
| `motifs.py` | Position weight matrices, motif discovery, conservation scoring |
| `restriction.py` | Restriction site mapping, virtual digests, and enzyme compatibility |

## Key Functions

| Function | Description |
|----------|-------------|
| `read_fasta()` | Parse a FASTA file into a dict of sequences |
| `write_fasta()` | Write sequences to FASTA format |
| `reverse_complement()` | Return the reverse complement of a DNA sequence |
| `gc_content()` | Fraction of G+C bases in a sequence |
| `gc_skew()` | (G-C)/(G+C) strand asymmetry metric |
| `melting_temperature()` | Estimate Tm using nearest-neighbor or Wallace method |
| `generate_consensus()` | Majority-rule consensus from aligned sequences |
| `consensus_with_ambiguity()` | Consensus using full IUPAC ambiguity codes |
| `bootstrap_consensus()` | Consensus with bootstrap confidence values |
| `count_kmers()` | Count all k-mers of a given length in a sequence |
| `find_microsatellites()` | Detect short tandem repeats (microsatellites) |
| `create_pwm()` | Build a position weight matrix from aligned motif instances |
| `discover_motifs()` | De novo motif discovery from a set of sequences |
| `find_restriction_sites()` | Locate enzyme recognition sites in a sequence |
| `virtual_digest()` | Simulate restriction enzyme digestion into fragments |

## Usage

```python
from metainformant.dna.sequence.core import read_fasta, reverse_complement, gc_content
from metainformant.dna.sequence.kmer import count_kmers, find_microsatellites
from metainformant.dna.sequence.restriction import find_restriction_sites, virtual_digest

seqs = read_fasta("sequences.fasta")
rc = reverse_complement("ATCGATCG")
kmers = count_kmers("ATCGATCGATCG", k=3)
sites = find_restriction_sites("ATCGAATTCGATCG", enzymes=["EcoRI"])
```
