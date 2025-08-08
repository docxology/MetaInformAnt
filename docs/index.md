### METAINFORMANT Documentation

Welcome to the METAINFORMANT docs. This site is organized by domain and core utilities. Use the navigation below to explore each analysis area.

- [Architecture](./architecture.md)
- [CLI](./cli.md)
- [Core Utilities](./core.md)

- DNA
  - [Overview](./dna/index.md)
  - [Sequences](./dna/sequences.md)
  - [Pairwise Alignment](./dna/alignment.md)
  - [Multiple Sequence Alignment (MSA)](./dna/msa.md)
  - [Phylogeny](./dna/phylogeny.md)
  - [Population Genetics](./dna/population.md)

- RNA
  - [Overview](./rna/index.md)
  - [amalgkit Wrapper](./rna/amalgkit.md)
  - [Workflow](./rna/workflow.md)
  - [Configs](./rna/configs.md)

- Simulation
  - [Overview](./simulation/index.md)
  - [Sequence Generators](./simulation/sequences.md)
  - [RNA Counts](./simulation/rna_counts.md)
  - [Agents & GridWorld](./simulation/agents.md)

- Math
  - [Overview](./math/index.md)
  - [Price Equation](./math/price.md)
  - [Selection Models](./math/selection.md)
  - [Drift–Diffusion Model](./math/ddm.md)

- Visualization
  - [Overview](./visualization/index.md)
  - [Phylogenetic Trees](./visualization/trees.md)
  - [Plots](./visualization/plots.md)
  - [Animations](./visualization/animations.md)

- Other Domains
  - [Protein](./protein/index.md)
  - [Ontology](./ontology/index.md)
  - [Phenotype](./phenotype/index.md)
  - [Epigenome](./epigenome/index.md)
  - [Ecology](./ecology/index.md)

```mermaid
graph TD
  A[CLI] --> B[DNA]
  A --> C[RNA]
  A --> D[Simulation]
  A --> E[Math]
  A --> F[Visualization]
  A --> G[Protein]
  A --> H[Ontology]
  A --> I[Phenotype]
  A --> J[Epigenome]
  A --> K[Ecology]
  subgraph Core
    L[config]
    M[io]
    N[logging]
    O[text]
    P[parallel]
    Q[hash]
  end
  L -.-> B
  M -.-> B
  L -.-> C
  M -.-> C
  L -.-> D
  M -.-> D
  L -.-> E
  M -.-> F
```


