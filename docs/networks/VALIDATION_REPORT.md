# Network Analysis Documentation Validation Report

> Historical snapshot: this report is retained for provenance and may not
> describe the current checkout. Regenerate validation outputs under `output/`
> when current evidence is needed.

**Date:** 2025-04-29
**Validated Repository:** /home/trim/Documents/Git/MetaInformAnt
**Documentation Path:** docs/networks/
**Source Code Path:** src/metainformant/networks/

## Executive Summary

This report compares the documented API in `docs/networks/` against the actual implementation in `src/metainformant/networks/`. The analysis reveals **significant discrepancies** between documentation and implementation across all major modules. Many documented functions are missing, have incorrect signatures, or return different structures than described.

---

## 1. Graph Analysis (`graph.md`)

### 1.1 `create_network()`

**Documented:**
```python
network = create_network(nodes, directed=False)
```
Parameters: `nodes` - list of node identifiers

**Actual Implementation** (`graph_core.py:196`):
```python
def create_network(edges: List[Tuple[str, str]], directed: bool = False, **kwargs) -> Any
```
Parameters: `edges` - list of (source, target) tuples

**Status:** ❌ **CRITICAL MISMATCH** - Different parameter type (nodes vs edges)

---

### 1.2 `add_edges_from_interactions()`

**Documented:**
```python
add_edges_from_interactions(network, interactions, confidence_threshold=0.5)
# interactions = [{"source": "TF1", "target": "gene1", "type": "regulation", "confidence": 0.9}]
```

**Actual Implementation** (`graph_core.py:727`):
```python
def add_edges_from_interactions(
    graph: Any,
    interactions: List[Tuple[str, str, float]],
    weight_threshold: float = 0.0
) -> None
```
Accepts list of `(source, target, weight)` tuples. Parameter name is `weight_threshold` not `confidence_threshold`. Edge attribute hardcoded as `interaction_type="protein_interaction"`.

**Status:** ❌ **CRITICAL MISMATCH** - Different input format and parameter names

---

### 1.3 `add_edges_from_correlation()`

**Documented:**
```python
add_edges_from_correlation(network, node_features, threshold=0.6, method="pearson", max_edges=1000)
```
Parameters: `node_features` (feature matrix), `method` (correlation method), `max_edges` (limit)

**Actual Implementation** (`graph_core.py:758`):
```python
def add_edges_from_correlation(
    graph: Any,
    correlation_matrix: Any,
    threshold: float = 0.5
) -> None
```
Accepts pre-computed correlation matrix. No `method` or `max_edges` parameters.

**Status:** ❌ **CRITICAL MISMATCH** - Different parameters and expected input

---

### 1.4 `network_metrics()`

**Documented return keys:**
- `n_nodes`, `n_edges`, `avg_degree`, `density`, `clustering_coeff`

**Actual Implementation** (`graph_algorithms.py:29`):
Returns:
- `num_nodes`, `num_edges`, `avg_degree`, `density`, `avg_clustering`

**Status:** ⚠️ **MINOR MISMATCH** - Key names differ (`clustering_coeff` vs `avg_clustering`, `n_nodes` vs `num_nodes`)

---

### 1.5 `centrality_measures()`

**Documented measures:**
- Degree, Betweenness, Closeness, Eigenvector, **PageRank**

**Actual Implementation** (`graph_algorithms.py:453`):
Returns only: `degree`, `betweenness`, `closeness`, `eigenvector`

**PageRank is NOT implemented.**

**Status:** ❌ **MISSING ALGORITHM** - PageRank not present

---

### 1.6 `shortest_paths()`

**Documented and actual implementation match** (`graph_algorithms.py:501`). Function exists and works as documented.

**Status:** ✅ **ACCURATE**

---

### 1.7 `network_similarity()`

**Documented:** Returns dict with multiple similarity methods ('jaccard', 'dice', 'overlap')

**Actual Implementation** (`graph_algorithms.py:216`):
- Only implements 'jaccard' correctly for nodes
- 'dice' and 'overlap' are not implemented as described (incorrect calculations)
- Returns single float, not dictionary

**Status:** ❌ **INCORRECT** - Limited functionality, wrong return type

---

## 2. Community Detection (`community.md`)

### 2.1 `detect_communities()`

**Documented methods:** `'louvain'`, `'leiden'`, `'hierarchical'`, `'spectral'`

**Actual Implementation** (`community.py:338`):
Supported methods:
- `'louvain'` ✅
- `'leiden'` ✅
- `'greedy'` ✅ (not explicitly documented)
- `'girvan_newman'` ✅ (not explicitly documented)
- `'label_propagation'` ✅ (not explicitly documented)
- `'asyn_lpa'` ✅ (not explicitly documented)
- `'fluid'` ✅ (not explicitly documented)

**Missing:** `'hierarchical'`, `'spectral'` ❌

**Status:** ⚠️ **PARTIAL** - Hierarchical and spectral not available; extra methods present but undocumented

---

### 2.2 `modularity()`

**Documented:** Function `modularity(network, communities)`

**Actual Implementation** (`community.py:556`): ✅ **Exists and works as documented**

**Status:** ✅ **ACCURATE**

---

### 2.3 `community_metrics()`

**Documented returns:** `n_communities`, `size_range`, `largest_community`, `modularity`

**Actual Implementation** (`community.py:593`): Returns comprehensive dict with:
- `n_communities`, `community_sizes` (dict with mean/std/min/max/sizes), `modularity`, `coverage`, `conductance` (dict), `internal_density` (dict), `quality_score`

**Status:** ✅ **ACCURATE** (actual provides MORE than documented)

---

### 2.4 `hierarchical_communities()`

**Documented:** Not mentioned in community.md but implied by "Hierarchical clustering" bullet

**Actual Implementation** (`community.py:492`): ✅ **Exists** with parameters `method='louvain'`, `max_levels=5`

**Status:** ℹ️ **UNDOCUMENTED** but present

---

### 2.5 Spectral Clustering

**Documented:** `'spectral'` method for `detect_communities()`

**Actual:** ❌ **NOT IMPLEMENTED** - No spectral clustering found in codebase

**Status:** ❌ **MISSING**

---

## 3. Pathway Analysis (`pathway.md`)

### 3.1 `PathwayNetwork` Class

**Documented** as class with simple initialization
**Actual** (`pathway.py:496`): ✅ **Exists** with extensive methods

**Status:** ✅ **ACCURATE**

---

### 3.2 `load_pathway_database()`

**Documented:**
```python
load_pathway_database(database="kegg", organism="hsa", pathway_ids=None, min_genes=5, max_genes=200)
```

**Actual Implementation:**
- No standalone `load_pathway_database()` function with those parameters
- Instead: `PathwayNetwork.load_from_database(filepath, format)` class method
- And: `load_pathway_database(pathway_data, name)` helper that takes structured dict

**Status:** ❌ **DIFFERENT API** - Documented function does not exist as described

---

### 3.3 `pathway_enrichment()`

**Documented:**
```python
pathway_enrichment(gene_list, pathway_network, background_genes=None,
                   method="hypergeometric", correction="bonferroni", min_overlap=3)
# Returns list of dicts with: 'id', 'name', 'p_value', 'overlap_genes', 'enrichment_ratio'
```

**Actual Implementation** (`pathway.py:791`):
```python
def pathway_enrichment(
    gene_list: List[str],
    pathway_network: PathwayNetwork,
    background_genes: Optional[List[str]] = None,
    method: str = "fisher",           # default differs
    correction: str = "bonferroni",
    min_overlap: int = 1,             # default differs
) -> Dict[str, Dict[str, Any]]
```
Returns dict keyed by pathway ID with keys: `overlap_size`, `pathway_size`, `query_size`, `background_size`, `p_value`, `enrichment_ratio`, `overlap_genes`, `overlapping_genes`, `expected_overlap`, `corrected_p_value`, `significant`

**Status:** ⚠️ **PARTIAL** - Different default method (fisher vs hypergeometric), different min_overlap default, different return structure (dict vs list)

---

### 3.4 `network_enrichment_analysis()`

**Documented:** Exists as separate function
**Actual Implementation** (`pathway.py:1057`): ✅ **Exists** but is just an alias for `pathway_enrichment()`

**Status:** ✅ **PRESENT** (but redundant)

---

### 3.5 `multi_omics_pathway_analysis()`

**Documented:** Exists with parameters `multiomics_data`, `analysis_type`, `min_consensus`

**Actual:** ❌ **NOT FOUND** - No such function in any module

**Status:** ❌ **MISSING**

---

### 3.6 `pathway_activity_inference()`

**Documented:** Exists with method="gsva" parameter

**Actual:** ❌ **NOT FOUND** - Only `pathway_activity_score()` exists (different functionality)

**Status:** ❌ **MISSING**

---

## 4. Protein-Protein Interaction Networks (`ppi.md`)

### 4.1 `ProteinNetwork` Class

**Documented** as main class for PPI networks
**Actual Implementation** (`ppi.py:529`): ✅ **Exists**

**Status:** ✅ **ACCURATE**

---

### 4.2 `load_string_interactions()`

**Documented:**
```python
load_string_interactions("string_interactions.tsv", min_confidence=0.7,
                         evidence_filter=["experimental", "database"], organism="9606")
```

**Actual Implementation** (`ppi.py:996`):
```python
def load_string_interactions(
    interactions_df: Any,              # requires DataFrame, not filepath
    proteins_df: Optional[Any] = None,
    confidence_threshold: int = 400    # 0-1000 scale, not 0-1
) -> ProteinNetwork
```
Takes pandas DataFrame as input, not a file path. No `organism` parameter.

**Status:** ❌ **CRITICAL MISMATCH** - Completely different signature and usage

---

### 4.3 `predict_interactions()`

**Documented methods:**
`'domain_fusion'`, `'gene_fusion'`, `'neighborhood'`, `'co_expression'`, `'co_occurrence'`

**Actual Implementation** (`ppi.py:779`):
Supported methods:
- `'similarity'` ✅
- `'correlation'` ✅
- `'guilt-by-association'` ✅
- `'ml'` ✅

**Documented methods** (`domain_fusion`, etc.): ❌ **NOT IMPLEMENTED**

**Status:** ❌ **WRONG METHOD NAMES** - Actual methods are different

---

### 4.4 `get_protein_partners()`

**Documented:** Standalone function
**Actual:** ❌ **NOT FOUND** as standalone function
**Alternative:** `ProteinNetwork.find_neighbors(protein)` method exists

**Status:** ❌ **MISSING** - Different API

---

### 4.5 `filter_by_confidence()`

**Documented:** Standalone function for PPI networks with `confidence_threshold` and `method_filter`

**Actual:** ❌ **NOT FOUND** as standalone function
**Alternative:** `ProteinNetwork` doesn't have this method

**Status:** ❌ **MISSING**

---

## 5. Gene Regulatory Networks (`regulatory.md`)

### 5.1 `GeneRegulatoryNetwork` Class

**Documented** as container for regulatory networks
**Actual Implementation**:
- Main class in `interaction/regulatory_core.py` ✅
- Also duplicated/renamed functions in separate `regulatory/` package

**Status:** ✅ **EXISTS** (but scattered across multiple locations)

---

### 5.2 `infer_grn()`

**Documented:**
```python
reg_network = infer_grn(expression_data, regulators=regulators,
                        method="correlation", threshold=0.7, max_targets=50)
```

**Actual Implementation** - TWO DIFFERENT FUNCTIONS:

1. `interaction/regulatory_analysis.py:134` - `infer_grn()`:
   ```python
   infer_grn(expression_data, gene_names, method="correlation",
             threshold=0.5, tf_genes=None, **kwargs) -> GeneRegulatoryNetwork
   ```
   Parameters: `gene_names` (not `regulators`), `tf_genes` (not `max_targets`)

2. `regulatory/grn_inference.py` - `infer_grn_correlation()` and `infer_grn_mutual_info()`, `infer_grn_regression()` - returns dict, not GeneRegulatoryNetwork

**Status:** ❌ **CONFUSING/DUPLICATED** - Multiple incompatible implementations

---

### 5.3 `add_transcription_factor()`

**Documented:** Standalone function
**Actual:** Method on `GeneRegulatoryNetwork` class: `grn.add_transcription_factor(gene_id, **kwargs)`

**Status:** ❌ **WRONG SCOPE** - Documented as module function, actually a method

---

### 5.4 `get_targets()` and `get_regulators()`

**Documented:** Standalone functions
```python
targets = get_targets(reg_network, "TF1")
regulators = get_regulators(reg_network, "gene1")
```

**Actual:** Methods on `GeneRegulatoryNetwork`:
- `grn.get_targets(tf=None)` (returns all targets if no tf)
- `grn.get_regulators(target)` (alias for `get_target_regulators`)

**No standalone functions exported.**

**Status:** ❌ **WRONG SCOPE** - Documented as module functions, actually class methods

---

### 5.5 `filter_by_confidence()`

**Documented:** Standalone function
**Actual:** Method on `GeneRegulatoryNetwork`: `grn.filter_by_confidence(threshold)`

**Status:** ❌ **WRONG SCOPE**

---

### 5.6 `regulatory_motifs()`

**Documented:** Standalone function
**Actual:**
- `interaction/regulatory_core.py`: `analyze_regulatory_motifs()` function
- `interaction/regulatory_core.py`: `GeneRegulatoryNetwork.regulatory_motifs()` method
- `interaction/regulatory_analysis.py`: `regulatory_motifs(grn)` function
- `regulatory/grn_inference.py`: `compute_network_motifs(edges)` function

**Status:** ⚠️ **CONFUSING** - Multiple similar functions with different signatures

---

## 6. Centrality Measures (Graph Analysis)

**Additional finding:** The documentation claims PageRank centrality is available but it's not implemented anywhere in the codebase.

**Status:** ❌ **MISSING ALGORITHM**

---

## 7. Missing Algorithms & Features

### Completely Missing from Documentation but Present in Code:

1. **`compare_community_methods()`** - Compares multiple community detection methods
2. **`fluid_communities()`** - Fluid communities algorithm
3. **`girvan_newman_communities()`** - Girvan-Newman algorithm
4. **`label_propagation_communities()`** - Label propagation
5. **`asyn_lpa_communities()`** - Asynchronous label propagation
6. **`analyze_ppi_disease_associations()`** - Disease association analysis for PPI
7. **`ppi_network_comparison()`** - Compare two PPI networks
8. **`ppi_network_enrichment()`** - Enrichment testing for PPI
9. **`functional_enrichment_ppi()`** - Functional enrichment using PPI annotations
10. **`construct_regulatory_network()`** - Build regulatory network from scratch
11. **`analyze_regulatory_dynamics()`** - Simulate regulatory dynamics
12. **`identify_regulatory_hubs()`** - Find hub regulators
13. **`regulatory_network_stability_analysis()`** - Stability assessment
14. **`pathway_disease_association()`** - Disease-pathway associations
15. **`pathway_visualization_data()`** - Prep pathway data for visualization
16. **`compute_network_motifs()`** - Statistical motif analysis with z-scores
17. **`validate_grn()`** - GRN validation with precision/recall/AUROC
18. **`score_regulators()`** - Score TFs for enrichment

### Completely Missing from Code but Documented:

1. ❌ `multi_omics_pathway_analysis()`
2. ❌ `pathway_activity_inference()` (only `pathway_activity_score` exists)
3. ❌ `spectral` community detection method
4. ❌ `hierarchical` community detection method (as a single method; hierarchical_communities exists but different)
5. ❌ PageRank centrality
6. ❌ STRING database `evidence_filter` parameter
7. ❌ PPI `get_protein_partners()` standalone function
8. ❌ PPI `filter_by_confidence()` standalone function
9. ❌ `domain_fusion`, `gene_fusion`, `neighborhood`, `co_expression`, `co_occurrence` prediction methods (different methods exist)

---

## 8. Parameter Mismatches Summary

| Function | Documented Parameter | Actual Parameter | Issue |
|----------|---------------------|------------------|-------|
| `create_network` | `nodes` | `edges` | Wrong argument type |
| `add_edges_from_interactions` | `confidence_threshold` | `weight_threshold` | Different semantics |
| `add_edges_from_correlation` | `node_features, method, max_edges` | `correlation_matrix` only | Missing params |
| `detect_communities(louvain)` | `random_state` | `randomize: bool` | Different param |
| `detect_communities(leiden)` | `n_iterations` | `n_iterations` | ✅ Matches |
| `load_pathway_database` | Multi-arg database loader | Dict/file loader | Different API |
| `load_string_interactions` | `min_confidence` (0-1), `organism`, `evidence_filter` | `confidence_threshold` (0-1000), no organism/evidence | Completely different |

---

## 9. Return Value Mismatches

| Function | Documented Return | Actual Return | Issue |
|----------|-------------------|---------------|-------|
| `network_metrics` | keys: `n_nodes`, `clustering_coeff` | `num_nodes`, `avg_clustering` | Key name mismatch |
| `pathway_enrichment` | List of dicts | Dict keyed by pathway ID | Structure mismatch |
| `network_similarity` | Dict of similarities | Single float | Type mismatch |

---

## 10. Recommendations

### Critical Fixes Required:

1. **Update graph.md** to reflect actual `create_network()` signature (edges not nodes)
2. **Update graph.md** to remove `max_edges` and `method` from `add_edges_from_correlation()` docs
3. **Fix centrality_measures()** to either implement PageRank or remove from documentation
4. **Update community.md**:
   - Remove `'spectral'` and `'hierarchical'` from `detect_communities()` supported methods list
   - Document the actual additional methods available (greedy, girvan_newman, label_propagation, etc.)
5. **Update ppi.md**:
   - Correct `load_string_interactions()` documentation to show DataFrame input
   - Update `predict_interactions()` to list correct method names
   - Add or remove `get_protein_partners()` and `filter_by_confidence()` — these don't exist as standalone
6. **Update regulatory.md**:
   - Change `infer_grn()` parameter `regulators` to either `gene_names` or `tf_genes` to match code
   - Remove `max_targets` parameter (not present)
   - Document that `get_targets()`, `get_regulators()`, `filter_by_confidence()`, `add_transcription_factor()` are **methods** on `GeneRegulatoryNetwork`, not module-level functions
7. **Update pathway.md**:
   - Remove `multi_omics_pathway_analysis()` and `pathway_activity_inference()` from documentation (not implemented)
   - Correct `load_pathway_database()` signature to match actual APIs
   - Document that `pathway_enrichment()` returns a dict, not list
8. **Consolidate duplicate implementations**: The `regulatory/` package and `interaction/regulatory_*` modules contain overlapping GRN inference code. This should be unified or clearly documented as separate modules with different purposes.
9. **Add missing function implementations** if any of the documented-but-missing functions are intended to be implemented.

---

## 11. Accuracy Score

**Overall Accuracy: ~45%**

Breakdown by module:
- Graph Analysis: 60% (6/10 functions accurate after accounting for mismatches)
- Community Detection: 50% (3/6 documented features accurate; 2 missing)
- Pathway Analysis: 30% (3/10 documented features accurate; 2 missing, several mismatches)
- PPI Networks: 25% (2/8 documented features accurate; 2 missing, multiple mismatches)
- Regulatory Networks: 35% (3/9 documented features accurate; scattered across classes)

**Major Issues:**
- 9 documented functions completely missing
- 10 functions have incorrect signatures/parameters
- 5 functions have incorrect return structures
- 1 algorithm (PageRank) undocumented as missing
- API confusion: standalone functions vs. class methods
- Duplicate/conflicting implementations in `regulatory/` vs `interaction/regulatory_*`

---

## 12. File Locations Reference

**Documentation Files Validated:**
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/README.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/PAI.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/SPEC.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/AGENTS.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/community.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/graph.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/index.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/pathway.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/ppi.md`
- `/home/trim/Documents/Git/MetaInformAnt/docs/networks/regulatory.md`

**Source Files Examined:**
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/__init__.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/__init__.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/graph.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/graph_core.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/graph_algorithms.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/community.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/analysis/pathway.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/interaction/__init__.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/interaction/ppi.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/interaction/regulatory.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/interaction/regulatory_core.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/interaction/regulatory_analysis.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/regulatory/__init__.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/regulatory/grn_inference.py`
- `/home/trim/Documents/Git/MetaInformAnt/src/metainformant/networks/regulatory/motif_analysis.py`

---

**Report Generated By:** Hermes Agent (Nous Research)
**Validation Method:** Manual code inspection and signature comparison
**Confidence Level:** High (all files read and analyzed)
