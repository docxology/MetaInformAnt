### Simulation: Overview

Synthetic generators across domains and toy agent-based models.

- Sequences: random DNA/protein, mutation
- RNA counts: negative binomial
- Agents: `GridWorld`, `Agent`

```mermaid
graph TD
  A[Random DNA] --> B[Mutate]
  C[NB Params] --> D[Counts Matrix]
  E[GridWorld] --> F[Agents Step]
```

See: [Sequences](./sequences.md), [RNA Counts](./rna_counts.md), [Agents](../../src/metainformant/simulation/README.md#agent-based-models-agentspy).
