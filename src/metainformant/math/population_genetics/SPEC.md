# Specification: population_genetics

## 🎯 Scope
Population genetics submodule.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 9 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `coalescent.py`
- `core.py`
- `demography.py`
- `effective_size.py`
- `fst.py`
- `ld.py`
- `selection.py`
- `statistics.py`

### Behavior Contracts

#### `fst.py`
- `fst_from_allele_freqs(pop1_freqs, pop2_freqs=None)` — moment F_ST pooled over
  loci: `sum_j var_p_j / sum_j (var_p_j + mean_i p_ij (1 - p_ij))`. For two
  populations this is algebraically identical to `(Ht - Hs) / Ht`, and the
  single-list `[p1, p2]` call equals the two-list call for one locus.
- `fst_from_allele_freq_matrix(pop_freqs)` — same estimator for k >= 2
  populations; `var_p` uses the unbiased `k/(k-1)` sample variance for `k > 2`
  and the population variance for `k == 2`.
- `weirs_fst(population_counts)` — Weir & Cockerham (1984) single-locus
  variance components (`a`, `b`, `c`) from per-population haplotype (gene-copy)
  counts; `c = 0` because haplotype counts carry no diploid heterozygosity, and
  `theta = sum_u a_u / sum_u (a_u + b_u + c_u)` clamped to [0, 1]. The legacy
  `(counts, labels)` CRC32-jitter signature was removed.

#### `ld.py`
- `ld_coefficients(pA, pa, pB, pb, pAB)` — frequency mode returns signed `(D, D')`.
- `ld_coefficients(genotypes)` — each row is one PHASED haplotype `[g1, g2]`;
  returns `{"D", "D_prime", "r_squared"}`.
- `ld_coefficients(genotypes, phased=False)` — each row is an unphased diploid
  `[[a1, a2], [b1, b2]]`; haplotype frequencies are recovered by a deterministic
  EM routine before computing D, D' and r^2.
