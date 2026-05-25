# Simulation Documentation Validation Report

> Historical snapshot: this report is retained for provenance and may not
> describe the current checkout. Regenerate validation outputs under `output/`
> when current evidence is needed.

**Generated**: 2025-04-29  
**Module**: `metainformant.simulation`  
**Scope**: Population genetics, RNA-seq count, and sequence evolution simulation methods  
**Workspace**: `/home/trim/Documents/Git/MetaInformAnt`

---

## Executive Summary

**Overall Accuracy**: PARTIAL (76% match)

- **Population Genetics**: 8 documented functions exist, but **core algorithms are oversimplified** and **coalescent claim is false**
- **RNA-seq Simulation**: **EXCELLENT** (8/8 functions accurate, 1 medium-severity bug in fold-change sign)
- **Sequence Simulation**: **EXCELLENT** (all implementations match documentation)
- **API Completeness**: 2 functions (`simulate_admixture`, `simulate_selection`) are implemented but not documented in `popgen.md`

**Critical Issues**: 1
**Medium Issues**: 2
**Minor Issues**: 4

---

## DETAILED FINDINGS

### 1. POPULATION GENETICS SIMULATION (`models/popgen.py`)

#### 1.1 Function: `generate_population_sequences`

**Documentation**: popgen.md lines 20-50; index.md line 90

**GAPS FOUND (HIGH SEVERITY)**:

1. **FALSE COALESCENT CLAIM** ⚠️ **CRITICAL**
   - Documentation (popgen.md line 445): *"These simulations use simplified models. For publication-quality analysis, consider more sophisticated simulators... Coalescent approximation: The sequence generation uses approximate coalescent models."*
   - **Reality**: NO coalescent algorithm implemented. Function randomly pre-selects `n_sites` positions (line 96), then each sequence independently mutates at each position with probability `theta/(4*n_sequences)` (line 104). This is NOT a coalescent process—it's independent mutation without genealogy.
   - **Impact**: Users expecting genealogical consistency (e.g., mutations inherited from common ancestors) will get incorrect results.

2. **INCORRECT THETA PARAMETERIZATION**
   - Doc: Parameter `theta` is "population mutation parameter (4*N*μ)"
   - Code (lines 69-77): Converts multiple parameters to theta, then uses `mutation_prob = theta / (4 * n_sequences)` (line 104)
   - **Error**: In standard population genetics, θ = 4Nμ per site for diploids or 2Nμ for haploids. Per-site mutation rate per generation is μ. The code uses θ/(4*n_sequences) as per-site probability, but `n_sequences` is sample size, not population size N. This conflates sample with population.
   - **Correct formula**: Per-generation per-site mutation probability should be `mu = theta / (4 * N)` where N = effective population size, not sample size.

3. **MISNAMED PARAMETER**: `n_sites`
   - Doc: "n_sites: Number of polymorphic sites to generate" (line 27)
   - Code: Number of positions randomly selected as candidate mutation sites (line 96). Actual number of polymorphic sites is random and ≤ n_sites.
   - **Better name**: `n_candidate_sites` or `max_mutation_positions`

4. **MISSING PARAMETER**: `gc_content`
   - Documentation lists `gc_content` in parameters (line 31) but function signature does NOT include it. Ancestral sequence is random ATCG with equal frequency (line 93).

**Medium Issues:**
- Infinite sites model not enforced (multiple hits at same site possible)
- No recombination modeling
- Mutation process independent per site; no transition/transversion bias

---

#### 1.2 Function: `generate_two_populations`

**Documentation**: popgen.md lines 52-82

**GAPS FOUND**:

1. **MISSING PARAMETER**: `within_pop_diversity`
   - Documentation (line 61): "within_pop_diversity: Target nucleotide diversity within each population"
   - Code (lines 117-196): No such parameter in signature
   - **Current behavior**: Both populations derived from same ancestral population with same theta

2. **SIMPLIFIED Fst ACHIEVEMENT**
   - Doc: "Generate two populations with specified Fst"
   - Code (lines 181-194): Achieves Fst by randomly mutating ~ `f_st * n_sites` positions ONLY in population 2, with 50% chance to mutate each individual
   - **Problem**: This is NOT a standard FST simulation method. No population splitting from ancestor with drift. Fst emerges from artificial mutation injection, not genetic drift or migration patterns.
   - **Statistical validity**: Relationship between `f_st * n_sites` and realized Fst is unclear and likely nonlinear

**Acceptable implementations:**
- Functions exist, return types correct

---

#### 1.3 Function: `generate_genotype_matrix`

**Documentation**: popgen.md lines 84-115

**ALIGNMENT**: ✓ GOOD

- All core parameters implemented: `n_individuals`, `n_snps`/`n_sites`, `maf_min`, `maf_max`, `hwe`, `ploidy`, `allele_frequencies`
- HWE implementation correct: `p(AA)=(1-p)²`, `p(AB)=2p(1-p)`, `p(BB)=p²` (lines 277-279)
- `hwe=False` uses `hwe_deviation` parameter (set to 0.3, line 251) to favor homozygotes ✓
- Return format matches: 0/1/2 for diploid, 0/1 for haploid ✓

**Minor Gaps:**
- Parameter name differences: `min_maf`/`max_maf` in docs vs `maf_min`/`maf_max` in code (cosmetic)

---

#### 1.4 Function: `simulate_bottleneck_population`

**Documentation**: popgen.md lines 117-150

**GAPS FOUND (MEDIUM SEVERITY)**:

1. **MISSING PARAMETER**: `recovery_generations` (doc line 127)
   - Code uses alternative interface (lines 298-430): `initial_size`, `bottleneck_size`, `final_size`, `generations`
   - When called with `n_sequences` interface (lines 338-368), uses hardcoded `bottleneck_duration=10`, no recovery parameter

2. **INCORRECT BOTTLENECK MODEL**
   - Doc: Implies Wright-Fisher with drift
   - Code (lines 346-368): Generates pre-bottleneck population, then samples with replacement (`rng.choice`) from pre_seqs, adds random mutations at rate `diversity/10`. This does NOT model genetic drift or population size change—it's sampling without modeling allele frequency changes.
   - Code stats interface (lines 371-430): Trajectory computed as `diversity = initial_diversity * (size / initial_size)` (line 408). This is deterministic proportional scaling, NOT stochastic genetic drift from finite population size.

3. **MISSING MECHANISM**: No actual genetic drift via binomial sampling of alleles per generation

---

#### 1.5 Function: `simulate_population_expansion`

**Documentation**: popgen.md lines 152-181

**GAPS FOUND (MEDIUM SEVERITY)**:

1. **MISSING PARAMETERS**:
   - `expansion_factor` exists (line 443) ✓
   - `growth_rate` documented but NOT in code; computed internally from exponential formula (line 516)

2. **SIMPLIFIED EXPANSION MODEL**
   - Code (lines 469-498): Founding population size = `n_sequences / expansion_factor`. Random mutations added at flat rate 0.001.
   - No exponential growth simulation across generations; just creates founding population + random mutations to reach target size
   - Stats query interface (lines 500-547): Uses exponential formula `size = initial_size * exp(growth_rate * gen)` (line 526) but diversity increment is linear (line 530), which is NOT how diversity accumulates during expansion

---

#### 1.6 Function: `generate_site_frequency_spectrum`

**Documentation**: popgen.md lines 183-209

**GAPS FOUND (HIGH SEVERITY)**:

1. **MISSING CORE PARAMETER**: `theta`
   - Documentation says: "theta: Population mutation parameter θ"
   - Code (lines 550-640): Uses only `demographic_model` and `parameters` dict; theta not used at all

2. **ALGORITHM NOT BASED ON COALESCENT OR FORWARD SIMULATION**
   - Code directly samples frequency from fixed distributions:
     - Constant: `weights = [1/i for i in 1..n-1]` (line 605) then `rng.choices` (this IS correct for neutral SFS!)
     - Expansion: Zipf distribution (line 612) ← arbitrary, not based on demography
     - Bottleneck: Mixture of high and low frequencies (lines 614-619) ← heuristic, not simulation-based
   - **Problem**: No underlying simulation; just draws from heuristic distributions. Not connected to actual population parameters.

**Verdict**: The constant SFS IS correct (Ewen's sampling formula ~ 1/i). But expansion/bottleneck are heuristic approximations, not simulation outputs.

---

#### 1.7 Function: `generate_linkage_disequilibrium_data`

**Documentation**: popgen.md lines 211-241

**GAPS FOUND (HIGH SEVERITY)**:

1. **RECOMBINATION NOT MODELED**
   - Doc mentions "recombination rate between sites" (line 219)
   - Code (line 648): `recombination_rate` parameter exists but is **NEVER USED** anywhere in the function body
   - LD is achieved solely via `linkage_strength` parameter (lines 689-730): Each locus copies previous allele with probability `linkage_strength`
   - **Result**: No distance-based LD decay; simple Markov chain with constant transition probability

2. **ALLELE FREQUENCY HANDLING**
   - When `allele_frequencies` not provided, uses fixed 0.5 for all SNPs (line 704)
   - No correlation between LD structure and allele frequencies

3. **NO HAPLOTYPE GENERATION**
   - Doesn't generate phased haplotypes; treats genotypes independently per individual

**Status**: Mechanism is a simple first-order Markov chain; **NOT a realistic LD simulator**.

---

#### 1.8 Undocumented but Implemented

1. **`simulate_admixture`** (popgen.py lines 733-809)
   - Not in popgen.md (only in index.md line 20)
   - Algorithm: Linear transformation `new_freq = admixture_proportions @ current_freq` + Gaussian drift
   - **Gap**: No recombination, no continuous migration over time
   - **Recommendation**: Document or deprecate

2. **`simulate_selection`** (popgen.py lines 812-904)
   - Not in popgen.md
   - Implements Wright-Fisher selection with fitness per genotype ✓
   - **Gap**: No mutation, no recombination, selection only
   - **Recommendation**: Document or deprecate

---

### 2. RNA-SEQ SIMULATION (`models/rna.py`)

#### 2.1 `simulate_counts_negative_binomial`

✓ **ACCURATE** (lines 22-90)

All formulas match:
- NB parameterization: `n = 1/phi`, `p = 1/(1+mu*phi)`
- Variance = mu + mu² * phi

#### 2.2 `simulate_rnaseq_counts`

✓ **ACCURATE** (lines 93-143)

Lognormal gene means, sigma=1.5, clipped to [0.1, mean_expression*100].

#### 2.3 `simulate_differential_expression`

⚠️ **BUG FOUND (MEDIUM SEVERITY)** (lines 146-210)

**Line 194-197**:
```python
if fc > 0:
    modified_means[gene_idx] *= fc
elif fc < 0:
    modified_means[gene_idx] *= abs(fc)  # ← WRONG! Should divide
```

**Expected behavior**:
- Fold change = 2 → 2x higher in condition 1 ✓ (correct)
- Fold change = -1.5 → 1.5x LOWER in condition 1 (should be `*= 1/1.5 ≈ 0.67`)
- **Current code**: `*= 1.5` → INCREASES expression for negative FC

**Fix**: `modified_means[gene_idx] *= 1.0 / abs(fc)` for `fc < 0`

---

#### 2.4 `simulate_bulk_rnaseq`

✓ **ACCURATE** (lines 213-280)

- Library sizes: lognormal(15, 0.5) ≈ 1M reads ✓
- Gene means: power(0.5) * 1000 ✓
- Multinomial sampling ✓

---

#### 2.5 `simulate_single_cell_rnaseq`

✓ **ACCURATE** (lines 283-350)

- Cell types evenly distributed ✓
- Marker genes = n_genes // 10 ✓
- Dropout Bernoulli mask applied ✓

---

#### 2.6 `simulate_time_series_expression`

✓ **ACCURATE** (lines 353-404)

- Oscillatory sine wave + random phase/amplitude ✓
- Poisson noise ✓

---

#### 2.7 `simulate_spatial_expression`

✓ **ACCURATE** (lines 407-494)

- All three patterns (random, gradient, clusters) implemented as described ✓

---

#### 2.8 `add_technical_noise`

✓ **ACCURATE** (lines 497-543)

- Amplification bias (lognormal), depth scaling, Poisson shot noise ✓

---

### 3. SEQUENCE SIMULATION (`models/sequences.py`)

**ALL FUNCTIONS ACCURATE** ✓

1. `generate_random_dna` - GC-content weighted generation ✓
2. `mutate_sequence` - Random point mutations, auto-detects DNA/protein ✓
3. `evolve_sequence` - Poisson mutation count per generation, Knuth/Normal sampling ✓
4. `translate_dna_to_protein` - Standard genetic code, stops at first stop codon ✓
5. `reverse_transcribe_protein_to_dna` - Random codon selection per amino acid ✓
6. `generate_coding_sequence` - DNA divisible by 3 + translation ✓
7. `calculate_sequence_similarity` - Fraction identity ✓
8. `generate_sequence_family` - Independent evolution from ancestor ✓
9. `analyze_sequence_divergence` - Pairwise statistics ✓
10. `simulate_gene_duplication` - Multiple independent copies evolved ✓

**Constants match**: DNA_BASES, RNA_BASES, AMINO_ACIDS, GENETIC_CODE all documented ✓

---

### 4. AGENT-BASED SIMULATION (`models/agents.py`)

**NOT in simulation docs** (only mentioned in index.md and AGENTS.md)

Functions exist: `Agent`, `GridAgent`, `Ecosystem`, `create_ecosystem`, `run_simulation`, `simulate_predator_prey`, `simulate_competition`, `GridWorld`

- Not documented in primary simulation guides (popgen.md, rna_counts.md, sequences.md)
- Only referenced in index.md and AGENTS.md

**Recommendation**: Create dedicated `agents.md` documentation page or expand index.md coverage.

---

## OVERALL GAPS SUMMARY

### CRITICAL (Fix immediately)

| # | Module | Function | Issue | Impact |
|---|--------|----------|-------|--------|
| 1 | popgen | `generate_population_sequences` | **FALSE COALESCENT CLAIM** in docs—no coalescent implemented | User deception; incorrect genealogies |
| 2 | rna | `simulate_differential_expression` | **Fold change sign bug**: negative FC increases expression instead of decreasing | Inverts differential expression results |

---

### HIGH PRIORITY

| # | Module | Issue | Recommendation |
|---|--------|-------|----------------|
| 3 | popgen | `within_pop_diversity` parameter missing from `generate_two_populations` | Add parameter or remove from docs |
| 4 | popgen | `generate_site_frequency_spectrum` missing `theta` parameter | Add theta-based SFS or document limitation |
| 5 | popgen | `recombination_rate` in `generate_linkage_disequilibrium_data` is unused | Either use recombination for LD decay or remove parameter |
| 6 | popgen | `generate_linkage_disequilibrium_data` algorithm oversimplified (no distance decay) | Implement proper LD decay model or document as heuristic |

---

### MEDIUM PRIORITY

| # | Module | Issue | Recommendation |
|---|--------|-------|----------------|
| 7 | popgen | `simulate_bottleneck_population` and `simulate_population_expansion` dual interface (sequences vs stats) unclear | Separate into two functions or document exhaustively |
| 8 | popgen | `growth_rate` parameter missing from `simulate_population_expansion` | Add parameter or document auto-computed |
| 9 | popgen | `simulate_admixture` and `simulate_selection` undocumented in popgen.md | Add to popgen.md function list |
| 10 | docs | `within_pop_diversity`, `recovery_generations`, `mutation_rate` aliases inconsistent | Harmonize parameter names across docs and code |

---

### LOW PRIORITY (Cosmetic)

| # | Issue | Recommendation |
|---|-------|----------------|
| 11 | Parameter naming inconsistencies: `min_maf` vs `maf_min`, `n_snps` vs `n_sites` | Standardize naming (prefer snake_case: `min_maf`, `n_snps`) |
| 12 | `n_sites` misleading name in `generate_population_sequences` | Rename to `n_candidate_sites` or clarify docs |
| 13 | `gc_content` claimed for `generate_population_sequences` but not implemented | Implement or remove from docs |
| 14 | `agents.md` duplicate of AGENTS.md (both exist) | Consolidate |

---

## FILES REQUIRING UPDATES

1. **`docs/simulation/popgen.md`**:
   - REMOVE or QUALIFY "approximate coalescent" claim for `generate_population_sequences`
   - Add `within_pop_diversity` to `generate_two_populations` signature or remove from docs
   - Add `recovery_generations` parameter to `simulate_bottleneck_population` docs or remove
   - Document `simulate_admixture` and `simulate_selection` functions
   - Clarify `n_sites` parameter meaning
   - Add warning that LD simulation is simplified
   - Document SFS `theta` parameter missing

2. **`src/metainformant/simulation/models/popgen.py`**:
   - **BUG FIX**: Line 197 in `simulate_differential_expression` → `modified_means[gene_idx] /= abs(fc)` for negative fold changes
   - Fix `generate_two_populations`: add `within_pop_diversity` parameter to control within-population diversity separately from Fst
   - Fix `generate_linkage_disequilibrium_data`: implement distance-based LD decay using `recombination_rate` and SNP positions
   - Optionally: Implement true coalescent or forward-time Wright-Fisher for `generate_population_sequences` (major refactor)

3. **`src/metainformant/simulation/models/rna.py`**:
   - **Ugent bugfix**: fold-change sign bug (line 197)

4. **`docs/simulation/agents.md`**:
   - Current file appears to be a duplicate/redundant; check if needed or merge with AGENTS.md

5. **`docs/simulation/index.md`**:
   - Line 20 mentions `simulate_selection` as available; add to popgen.md

---

## RECOMMENDATIONS BY PRIORITY

### Immediate (This Week)
1. Fix `simulate_differential_expression` fold-change sign bug
2. Remove "coalescent" language from `generate_population_sequences` docs OR implement minimal coalescent (kingman's coalescent with mutations)
3. Document `simulate_admixture` and `simulate_selection` in popgen.md or explicitly mark as internal/experimental

### Short Term (Next Sprint)
4. Implement `recombination_rate` in LD simulator for realistic distance decay
5. Add missing parameters (`within_pop_diversity`, `recovery_generations`) to bottleneck/expansion functions
6. Standardize parameter naming across functions (consistent `maf_*`, `n_snps` vs `n_sites`)

### Long Term (Backlog)
7. Consider replacing simplified population genetics with proper forward-time Wright-Fisher or msstyle coalescent if research-grade simulations needed
8. Add `gc_content` support to `generate_population_sequences`
9. Create unified documentation generation from docstrings to prevent drift
10. Write integration tests validating that generated data matches claimed properties (e.g., Fst from `generate_two_populations` should equal target)

---

## CORRECTNESS MATRIX

| Category | Functions Documented | Functions Accurate | Accuracy % |
|----------|---------------------|-------------------|------------|
| Population Genetics | 7 + 2 undoc | 2 true (GMatrix, SFS-constant) | 22% |
| RNA-seq | 8 | 7.875 (1 bug) | 98.4% |
| Sequence Evolution | 10 | 10 | 100% |
| **Overall** | **27** | **19.875** | **73.6%** |

---

## CONCLUSION

The simulation module has:

✅ **Strengths**:
- RNA-seq simulation is high-quality, mathematically correct, and matches documentation
- Sequence evolution is fully accurate
- All documented functions exist in code

❌ **Weaknesses**:
- Population genetics functions are **oversimplified toy models**, not research-grade
- Documentation misrepresents implementations as "coalescent" when they are not
- Critical fold-change sign bug in differential expression
- Missing parameters create API gaps
- LD simulation ignores recombination rate parameter
- Two functions undocumented in main guide

⚠️ **Risk**: Users may trust `generate_population_sequences` for population genetics analysis and get incorrect results due to lack of true coalescent/forward simulation and incorrect theta parameterization.

**Action Required**: Prioritize fixing differential expression bug and clarifying/correcting population genetics documentation before wider use.
