# Architecture

High-level architecture of METAINFORMANT.

```mermaid
flowchart LR
  subgraph CLI
    A[metainformant.__main__\nargparse CLI]
  end
  subgraph Core[Core Utilities]
    C1[config]
    C2[io]
    C3[logging]
    C4[text]
    C5[parallel]
    C6[hash]
    C7[paths]
    C8[cache]
    C9[db]
  end
  subgraph DNA
    D1[sequences]
    D2[alignment]
    D3[msa]
    D4[phylogeny]
    D5[population]
    D6[entrez/ncbi]
    D7[variants]
  end
  subgraph RNA
    R1[amalgkit]
    R2[workflow]
    R3[configs]
  end
  subgraph GWAS[Genome-Wide Association]
    G1[quality]
    G2[structure]
    G3[association]
    G4[visualization]
  end
  subgraph SingleCell[Single-Cell Genomics]
    SC1[preprocessing]
    SC2[dimensionality]
    SC3[clustering]
    SC4[trajectory]
    SC5[visualization]
    SC6[integration]
  end
  subgraph Quality[Quality Control]
    Q1[fastq]
  end
  subgraph Simulation
    S1[sequences]
    S2[rna]
    S3[agents]
  end
  subgraph Math
    M1[price]
    M2[selection]
    M3[ddm]
  end
  subgraph Viz
    V1[trees]
    V2[plots]
    V3[animations]
  end
  A --> D1 & D2 & D3 & D4 & D5
  A --> R1 & R2
  A --> G1 & G2 & G3 & G4
  A --> SC1 & SC2 & SC3 & SC4 & SC5 & SC6
  A --> Q1
  A --> S1 & S2 & S3
  A --> M1 & M2 & M3
  A --> V1 & V2 & V3
  C1 -.-> A
  C2 -.-> A
  C1 -.-> D1
  C2 -.-> D1
  C1 -.-> R2
  C2 -.-> R1
  C2 -.-> SC1
  C5 -.-> SC2
  C8 -.-> SC6
  C2 -.-> Q1
  C5 -.-> R2
  C8 -.-> R2
```

## Project directories and conventions

- **`config/`**: Declarative configuration and options for runs. Read by `metainformant.core.config` and consumed across domains. Environment variables may override values.
- **`data/`**: Canonical datasets and local databases. Treated as read-mostly inputs and long-lived artifacts under versioned subfolders.
- **`output/`**: All run and test outputs. Ephemeral, reproducible, safe to delete. Modules must default to writing here unless a user-specified path is provided.

Guidelines:

- `core` owns config loading and path resolution; domain modules do not hardcode absolute paths.
- Prefer parameters and env variables to override locations, but default to `config/`, `data/`, and `output/` at repo root.
- Tests and CLI runs should not write outside `output/`.

See also: [CLI](./cli.md), [Core](./core.md).
