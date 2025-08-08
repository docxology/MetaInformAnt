### DNA: Overview

Capabilities
- FASTA I/O (`sequences.read_fasta`)
- Pairwise alignment (`alignment.global_align`, `alignment.local_align`)
- Lightweight MSA (`msa.align_msa`, optional external tools)
- Phylogeny (NJ/UPGMA, Newick export)
- Population genetics (diversity, Tajima's D, Fst)
 - Accession checks and Entrez fetch (NCBI)
- Consensus sequence from alignment
- FASTQ quality summaries

```mermaid
flowchart TD
  A[FASTA] --> B[MSA]
  A --> C[Pairwise Alignments]
  B --> D[Distance Matrix]
  D --> E[NJ/UPGMA Tree]
  A --> F[Population Stats]
```

See: [Sequences](./sequences.md), [Alignment](./alignment.md), [MSA](./msa.md), [Phylogeny](./phylogeny.md), [Population](./population.md).

- Extras
  - [Entrez/NCBI](./ncbi.md)
  - [Accession Validation](./accessions.md)
  - [Consensus](./consensus.md)
  - [FASTQ](./fastq.md)


