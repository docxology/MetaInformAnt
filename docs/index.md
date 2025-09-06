### METAINFORMANT Documentation

Welcome to the METAINFORMANT docs. This site is organized by domain and core utilities. Use the navigation below to explore each analysis area.

- [Architecture](./architecture.md)
- [CLI](./cli.md)
- [Core Utilities](./core.md)
- [Setup](./setup.md)
- [Testing](./testing.md)

- DNA
  - [Overview](./dna/index.md)
  - [Sequences](./dna/sequences.md)
  - [Pairwise Alignment](./dna/alignment.md)
  - [Multiple Sequence Alignment (MSA)](./dna/msa.md)
  - [Phylogeny](./dna/phylogeny.md)
  - [Population Genetics](./dna/population.md)

- RNA
  - [Overview](./rna/index.md)
  - [amalgkit Wrapper](./rna/amalgkit/amalgkit.md)
  - [Workflow](./rna/workflow.md)
  - [Configs](./rna/configs.md)
  - [Steps](./rna/steps.md)

- Single-Cell Genomics
  - [Overview](./singlecell/index.md)
  - [Preprocessing](./singlecell/preprocessing.md)
  - [Dimensionality Reduction](./singlecell/dimensionality.md)
  - [Clustering](./singlecell/clustering.md)
  - [Trajectory Analysis](./singlecell/trajectory.md)
  - [Visualization](./singlecell/visualization.md)
  - [Integration](./singlecell/integration.md)

- Quality Control
  - [Overview](./quality/index.md)
  - [FASTQ Analysis](./quality/fastq.md)

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

See also: the top-level project README for quickstarts.

```mermaid
graph TD
  A[CLI] --> B[DNA]
  A --> C[RNA]
  A --> SC[Single-Cell]
  A --> QC[Quality Control]
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
    R[paths]
    S[cache]
    T[db]
  end
  L -.-> B
  M -.-> B
  L -.-> C
  M -.-> C
  L -.-> SC
  M -.-> SC
  L -.-> QC
  M -.-> QC
  L -.-> D
  M -.-> D
  L -.-> E
  M -.-> F
  P -.-> SC
  S -.-> SC
```
