# Information Decomposition

Partial information decomposition (PID), co-information, unique information, redundant information, synergistic information, dual total correlation, and O-information for multivariate biological data analysis.

## Key Concepts

**Co-information** (interaction information) generalises mutual information to N variables via inclusion-exclusion on entropies. For three variables: CI(X;Y;Z) = I(X;Y) - I(X;Y|Z). Positive co-information indicates redundancy; negative indicates synergy.

**Partial information decomposition** (PID) decomposes the total information I(X,Y;T) that two source variables carry about a target into four non-negative components: redundancy (shared by both sources), unique_X (only from X), unique_Y (only from Y), and synergy (only from joint observation). Uses the minimum mutual information (MMI) approach from Barrett (2015).

**Unique information** UI(X->T\\Y) = I(X;T) - Redundancy quantifies the information that source X exclusively provides about target T, beyond what source Y also provides.

**Redundant information** equals the minimum individual mutual information: Red = min_i I(S_i;T). This is the common information floor that every source individually provides.

**Synergistic information** is information about the target available only from jointly observing all sources. For two sources: Syn = I(X,Y;T) - I(X;T) - I(Y;T) + Redundancy. For N sources: Syn = I(S_joint;T) - sum(I(S_i;T)) + (n-1)*Redundancy.

**Dual total correlation** (binding information) DTC = H(X_joint) - sum(H(X_i|X_{not_i})) measures shared multivariate dependence. Non-negative, zero iff at least one variable is independent of the rest.

**O-information** (Rosas et al. 2019) quantifies the balance between redundancy and synergy in a multivariate system: O = TC - DTC. Positive values indicate redundancy-dominated systems; negative values indicate synergy-dominated systems.

## Function Reference

### `co_information(variables, base=2.0) -> float`

Compute co-information via inclusion-exclusion on entropies.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `variables` | `List[Sequence[Any]]` | At least 2 discrete variable sequences (same length) |
| `base` | `float` | Logarithm base |

**Returns** co-information value. Positive = redundancy, negative = synergy.

### `partial_information_decomposition(x, y, target, base=2.0) -> Dict`

Decompose I(X,Y;T) into redundant, unique, and synergistic components using MMI.

**Returns** dict with:
- `redundancy`: Shared information both sources carry
- `unique_x`: Information only X carries
- `unique_y`: Information only Y carries
- `synergy`: Information from joint observation only
- `total`: Total mutual information I(X,Y;T)

### `unique_information(source, target, other_source, base=2.0) -> float`

Compute unique information that source provides about target beyond other_source.

### `redundant_information(sources, target, base=2.0) -> float`

Compute redundant information shared by all sources about target. Equals min_i I(S_i;T) under MMI.

### `synergistic_information(sources, target, base=2.0) -> float`

Compute synergistic information available only from joint observation of all sources. For N sources: Syn = I(S_joint;T) - sum(I(S_i;T)) + (n-1)*Redundancy.

### `dual_total_correlation(variables, base=2.0) -> float`

Compute dual total correlation (binding information). DTC = (1-n)*H(X_joint) + sum(H(X_{not_i})).

### `o_information(variables, base=2.0) -> float`

Compute O-information: O = (n-2)*H(X_joint) + sum(H(X_i)) - sum(H(X_{not_i})). Positive = redundancy-dominated, negative = synergy-dominated.

## Usage Examples

```python
from metainformant.information.metrics.advanced.decomposition import (
    co_information, partial_information_decomposition,
    unique_information, redundant_information,
    synergistic_information, dual_total_correlation,
    o_information,
)

# Gene expression variables and a phenotype target
gene_a   = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]
gene_b   = [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0]
phenotype = [0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0]

# Partial information decomposition
pid = partial_information_decomposition(gene_a, gene_b, phenotype)
print(f"Redundancy: {pid['redundancy']:.4f} bits")
print(f"Unique (A): {pid['unique_x']:.4f} bits")
print(f"Unique (B): {pid['unique_y']:.4f} bits")
print(f"Synergy:    {pid['synergy']:.4f} bits")
print(f"Total:      {pid['total']:.4f} bits")

# Unique information: what does gene_a uniquely provide?
ui_a = unique_information(gene_a, phenotype, gene_b)
print(f"Unique info (A beyond B): {ui_a:.4f} bits")

# Redundant information across three regulatory genes
gene_c = [1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0]
red = redundant_information([gene_a, gene_b, gene_c], phenotype)
print(f"Redundancy (3 sources): {red:.4f} bits")

# Synergistic information
syn = synergistic_information([gene_a, gene_b], phenotype)
print(f"Synergy: {syn:.4f} bits")

# Co-information for three variables
ci = co_information([gene_a, gene_b, gene_c])
if ci > 0:
    print(f"Redundancy-dominated: CI = {ci:.4f} bits")
else:
    print(f"Synergy-dominated: CI = {ci:.4f} bits")

# Dual total correlation
dtc = dual_total_correlation([gene_a, gene_b, gene_c])
print(f"Dual total correlation: {dtc:.4f} bits")

# O-information: redundancy vs synergy balance
o_val = o_information([gene_a, gene_b, gene_c])
print(f"O-information: {o_val:.4f} ({'redundancy' if o_val > 0 else 'synergy'})")
```

## Biological Applications

**Gene regulatory networks**: PID reveals whether transcription factors provide redundant or synergistic information about target gene expression. Synergy-dominated interactions suggest combinatorial regulation.

**Biomarker selection**: Unique information identifies features that contribute non-overlapping predictive power, guiding panel design.

**Neural coding**: O-information quantifies whether neural populations encode information redundantly (robust coding) or synergistically (efficient coding).

**Multi-omics integration**: Decomposing information across genomic, transcriptomic, and epigenomic layers reveals which layers contribute unique versus shared predictive power for clinical outcomes.

## Configuration

Environment variable prefix: `INFO_`

## Related Modules

- `metainformant.information.metrics.core.syntactic` -- entropy, MI, conditional MI
- `metainformant.information.metrics.core.estimation` -- bias-corrected estimation
- `metainformant.information.metrics.advanced.channel` -- channel capacity
- `metainformant.information.metrics.advanced.geometry` -- divergence measures, Fisher-Rao
