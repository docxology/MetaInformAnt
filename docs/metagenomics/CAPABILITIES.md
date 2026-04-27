# Metagenomics Module Capabilities

Complete function and class reference for the metagenomics module, organized by subpackage.

---

## Table of Contents

1. [Amplicon (`amplicon/`)](#amplicon)
2. [Shotgun (`shotgun/`)](#shotgun)
3. [Diversity (`diversity/`)](#diversity)
4. [Functional (`functional/`)](#functional)
5. [Comparative (`comparative/`)](#comparative)
6. [Visualization (`visualization/`)](#visualization)
7. [Algorithm Comparison](#algorithm-comparison)
8. [Function Matrix](#function-matrix)

---

## Amplicon

### Core Data Structures

#### `OTU`

**Dataclass**: Operational Taxonomic Unit from centroid-based clustering.

**Fields**:
- `centroid_id` (str): Identifier of the representative sequence
- `centroid_sequence` (str): The nucleotide sequence
- `member_ids` (list[str]): All sequence IDs assigned to this OTU
- `member_sequences` (list[str]): All sequences in OTU
- `size` (int): Number of members

**Example**:
```python
otu = OTU(
    centroid_id="OTU_1",
    centroid_sequence="ACGTACGT...",
    member_ids=["seq1", "seq2", "seq3"],
    size=3,
)
```

---

#### `ClusteringResult`

**Dataclass**: Result of `cluster_otus()`.

**Fields**:
- `otus` (list[OTU]): List of OTU objects
- `threshold` (float): Identity threshold used (e.g., 0.97)
- `total_sequences` (int): Input sequences count
- `num_otus` (int): Number of OTUs formed
- `otu_table` (dict[str, int]): `{centroid_id: count}` — abundance per OTU

**Usage**:
```python
result = cluster_otus(sequences, threshold=0.97)
print(f"{result.num_otus} OTUs from {result.total_sequences} sequences")
for otu in result.otus[:5]:
    print(f"  {otu.centroid_id}: size={otu.size}")
```

---

#### `TaxonomyAssignment`

**Dataclass**: Taxonomic classification of a sequence.

**Fields**:
- `sequence_id` (str): Identifier of classified sequence
- `lineage` (list[tuple[str, str]]): Ordered list of `(rank, name)` tuples, e.g., `[("domain","Bacteria"), ("phylum","Firmicutes"), ("class","Bacilli")]`
- `confidence` (dict[str, float]): Per-rank confidence scores (0–1)
- `method` (str): Classification method used (`"naive_bayes"`, `"blast"`)
- `best_hit_id` (str): Best-matching reference sequence ID
- `best_hit_identity` (float): Identity to best hit (0–1)

---

### Functions

#### `cluster_otus(sequences, threshold=0.97, abundance=None, sort_by="length", prefilter=True) -> ClusteringResult`

**Purpose**: Greedy centroid-based OTU clustering (VSEARCH/UCLUST style).

**Algorithm**:
1. Sort sequences by abundance (descending) or length (descending)
2. First sequence → centroid of OTU #1
3. For each subsequent sequence:
   - Compare to all existing centroids using `calculate_identity()`
   - If identity ≥ threshold for any centroid, assign to that OTU (best match)
   - Else, create new OTU with this sequence as centroid

**Pre-filter**: Uses `_quick_identity()` (8-mer overlap) to skip expensive alignment if k-mer similarity < threshold - 0.10.

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: nucleotide_sequence}`
- `threshold` (float): Identity threshold 0.0–1.0 (default 0.97 ≈ 97% identity)
- `abundance` (dict[str, int] | None): Precomputed sequence abundances for sorting; if None and `sort_by="abundance"`, falls back to length sorting
- `sort_by` (str): `"length"` (default) or `"abundance"` centroid seeding order
- `prefilter` (bool): Enable k-mer pre-filter for speed (default True)

**Returns**: `ClusteringResult`

**Complexity**: O(n × c) where c = average candidates examined per sequence after pre-filter (~1–10 vs. n without filter).

**Example**:
```python
seqs = {
    "s1": "ACGTACGTACGT",
    "s2": "ACGTACGTACGA",  # 1 mismatch
    "s3": "TTTTAAAACCCC",
}
result = cluster_otus(seqs, threshold=0.90)
print(result.num_otus)  # Likely 2: {s1,s2} cluster, s3 separate
```

**When to use**: Legacy 16S analysis where 97% OTUs are desired; large datasets where ASV denoising is too slow.

---

#### `calculate_identity(seq1, seq2) -> float`

**Purpose**: Compute pairwise global sequence identity using Needleman-Wunsch with affine gap penalties.

**Scoring**:
- Match: +2
- Mismatch: -1
- Gap open: -5
- Gap extend: -1

**Identity definition**: `(# identical positions in alignment) / (alignment length)`.

**Parameters**:
- `seq1`, `seq2` (str): Uppercase nucleotide strings (ACGT)

**Returns**: `float` in [0.0, 1.0]

**Example**:
```python
identity = calculate_identity("ACGTACGT", "ACGTACGA")  # 7/8 = 0.875
```

**Notes**:
- Strips any existing gap characters (`-`, `.`) before alignment
- Raises `ValueError` if either sequence empty
- Full traceback computed; identity counted during traceback

**When to use**: Called internally by `cluster_otus()`; exported for direct pairwise comparisons.

---

#### `filter_chimeras(sequences, reference_db=None, abundance=None, min_divergence=1.0, min_score=0.28) -> dict[str, bool]`

**Purpose**: Detect chimeric sequences using UCHIME de novo algorithm.

**Algorithm**:
1. Sort sequences by abundance (most abundant first; assumed non-chimeric)
2. For each query sequence (starting from 3rd most abundant):
   - Find two best-matching parent candidates among more abundant sequences using `_quick_identity()`
   - For each possible breakpoint (left/right split), build chimeric model from parents:
     - Model A: left from parent A, right from parent B
     - Model B: left from parent B, right from parent A
   - Compute divergence between query and optimal chimeric model
   - Chimera score = divergence_ratio; flag if score ≥ `min_score`

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: sequence}`
- `reference_db` (dict[str, str] | None): If provided, use these as parent pool instead of abundance-sorted input (reference-based chimera check)
- `abundance` (dict[str, int] | None): Abundance counts for sorting; if None, uses length
- `min_divergence` (float): Minimum divergence ratio between chimeric model and parents (default 1.0)
- `min_score` (float): Threshold to flag as chimeric (0.0–1.0, default 0.28)

**Returns**: `dict[str, bool]` — `True` = chimeric, `False` = non-chimeric

**Example**:
```python
seqs = {
    "parent1": "AAAAAAAAAA",
    "parent2": "TTTTTTTTTT",
    "chimera": "AAAAATTTTT",  # hybrid of parent1 + parent2
}
result = filter_chimeras(seqs)
print(result["chimera"])  # True (likely flagged)
print(result["parent1"])  # False (most abundant → assumed good)
```

**Note**: Chimera detection is approximate; known false positives in homopolymer regions. Always inspect flagged sequences manually in production.

---

#### `denoise_sequences(sequences, error_rates=None, quality_scores=None, omega_a=1e-40, min_abundance=1) -> DenoisingResult`

**Purpose**: DADA2-style denoising to resolve exact Amplicon Sequence Variants (ASVs) from amplicon reads, modeling PCR and sequencing errors.

**Algorithm summary** (not full DADA2, simplified version):
1. **Dereplication**: Count unique sequences and abundances
2. **Error rate estimation** (if `quality_scores` provided): Learn position-specific substitution rates from Phred quality scores
3. **Merging error variants**: For each sequence, compute probability that it arose as an error from a more abundant sequence using Poisson model; merge if probability exceeds threshold
4. **Output**: List of denoised ASVs with final abundances

**Parameters**:
- `sequences` (dict[str, str]): Read sequences (keys arbitrary, values DNA strings)
- `error_rates` (dict | None): Pre-computed error model; if None, estimated from data or default rates used
- `quality_scores` (dict[str, str] | None): `{read_id: Phred string}` for accurate error rates
- `omega_a` (float): Parameter controlling merging aggressiveness (lower = more aggressive merging)
- `min_abundance` (int): Minimum final abundance for an ASV to be retained

**Returns**: `DenoisingResult` with `asvs`, `abundances`, `error_model`, `num_asvs`, `reads_removed`

**Example**:
```python
reads = {"r1": "ACGTACGT", "r2": "ACGTACGA", "r3": "ACGTACGA", ...}
result = denoise_sequences(reads)
print(f"Input: {len(reads)} reads → {result.num_asvs} ASVs")
for asv, count in zip(result.asvs, result.abundances):
    print(f"  {asv}: {count} reads")
```

**Recommendation**: Use ASV denoising for high-resolution microbiome analysis (e.g., strain-level tracking, single-nucleotide variants). Requires higher sequencing depth than OTU clustering.

---

#### `merge_paired_reads(forward, reverse, min_overlap=20, max_mismatch_ratio=0.2) -> dict[str, str]`

**Purpose**: Merge overlapping paired-end reads into a single consensus sequence.

**Algorithm**:
1. For each read pair, find largest overlap ≥ `min_overlap`
2. Align overlapping region with allowed mismatch ratio
3. In conflict positions (mismatch or gap), choose base from read with higher Phred quality (if available)
4. Concatenate non-overlapping parts

**Parameters**:
- `forward` (dict[str, str]): `{id: sequence}` (R1 reads)
- `reverse` (dict[str, str]): `{id: reverse-complemented sequence}` — caller must RC reverse reads first
- `min_overlap` (int): Minimum overlapping bases required (default 20)
- `max_mismatch_ratio` (float): Maximum allowed mismatches in overlap (default 0.2 = 20%)

**Returns**: `dict[str, str]` merged sequences keyed by read ID

**Example**:
```python
fwd = {"pair1": "ACGTACGTACGTACGT"}
rev_rc = {"pair1": "ACGTACGTACGTACGT"[::-1].translate(complement_table)}  # RC
merged = merge_paired_reads(fwd, rev_rc)
print(merged["pair1"])  # "ACGTACGTACGTACGT" (perfect merge)
```

---

#### `classify_taxonomy(sequences, reference_db, reference_taxonomy=None, method="naive_bayes", confidence_threshold=0.8, k=8, bootstrap_n=100) -> list[TaxonomyAssignment]`

**Purpose**: Assign taxonomic lineage to sequences using reference database.

**Methods**:

1. **`naive_bayes`** (RDP-style):
   - Extract k-mers (default k=8) from query sequence
   - For each k-mer, compute conditional probability `P(taxon | kmer)` from training database
   - Multiply probabilities across all k-mers (naive assumption of independence)
   - Assign at each rank (domain → species) with confidence = bootstrap proportion

2. **`blast`** (k-mer overlap):
   - Find top N most similar reference sequences by k-mer Jaccard
   - Take LCA of top hits as assignment

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: sequence}`
- `reference_db` (dict[str, str]): training database `{ref_id: sequence}`
- `reference_taxonomy` (dict[str, list]): `{ref_id: [(rank, name), ...]}` required for `naive_bayes`; optional for `blast` (can use LCA from hits)
- `method` (str): `"naive_bayes"` (default) or `"blast"`
- `confidence_threshold` (float): Minimum per-rank confidence for assignment (default 0.8). Ranks below threshold get `"unclassified"`
- `k` (int): k-mer size (default 8)
- `bootstrap_n` (int): Number of bootstrap replicates for confidence estimation (default 100, slower)

**Returns**: `list[TaxonomyAssignment]` one per input sequence

**Example**:
```python
ref_db = {"ref1": "ACGT...", "ref2": "TTAA..."}
ref_tax = {
    "ref1": [("domain","Bacteria"), ("phylum","Actinobacteria")],
    "ref2": [("domain","Bacteria"), ("phylum","Bacteroidetes")],
}
assignments = classify_taxonomy(
    {"query1": "ACGT..."},
    ref_db, ref_tax,
    method="naive_bayes",
    confidence_threshold=0.7,
)
assignments[0].lineage
# [("domain","Bacteria"), ("phylum","Actinobacteria")]
```

**Performance**: `naive_bayes` builds k-mer → taxon frequency table from reference (slow one-time cost); classification is fast O(k-mers). `blast` is O(n × m) over references.

---

## Shotgun

### Assembly Data Structures

#### `Contig`

**Fields**:
- `contig_id` (str): Unique identifier (e.g., `"contig_001"`)
- `sequence` (str): Nucleotide sequence (capital ACGT)
- `length` (int): Sequence length in bp
- `coverage` (float): Mean read depth across contig
- `gc_content` (float): Fraction G+C (0–1)

---

#### `AssemblyStats`

**Fields**:
- `total_contigs` (int)
- `total_length` (int) — sum of all contig lengths
- `n50` (int) — shortest contig such that ≥50% total assembly is in contigs ≥ that length
- `l50` (int) — fewest contigs accounting for ≥50% total length
- `n90` (int) — 90% summary metric
- `l90` (int) — complementary to N90
- `gc_content` (float) — weighted average GC%
- `mean_coverage` (float)
- `largest_contig` (int)
- `smallest_contig` (int)

**Calculation**:
Sort contigs by length descending; accumulate lengths until reaching 50% of total; N50 = length of last contig added; L50 = count of contigs used.

---

### Assembly Functions

#### `assemble_contigs(reads, k_range=None, min_contig_length=200, min_kmer_coverage=2) -> list[Contig]`

**Purpose**: De Bruijn graph-based de novo assembly.

**Algorithm** (conceptual, simplified implementation):
1. For each k in `k_range`:
   - Extract all k-mers from reads with counts
   - Build de Bruijn graph: nodes = (k-1)-mers, edges = k-mers
   - Tip removal: Delete low-coverage dead-ends (< `min_kmer_coverage`)
   - Bubble popping: Merge alternative paths with similar coverage
   - Extract contigs as non-branching paths
2. Merge contigs across k values (remove duplicates, keep longest)

**Parameters**:
- `reads` (dict[str, str] | list[str]): Read sequences
- `k_range` (list[int] | None): k-mer sizes; default `[21, 33, 55]` adjusted by read length
- `min_contig_length` (int): Discard contigs shorter than this (default 200)
- `min_kmer_coverage` (int): Minimum k-mer count to keep node (default 2)

**Returns**: `list[Contig]`

**Complexity**: O(number of k-mers × average node degree). k-mer counting is O(total bases). Memory ~ #unique k-mers × (k bytes + count). For 1 Gb reads, 31-mers: ~10⁸ unique kmers → several GB.

**Limitations**:
- No repeat resolution (unitigs only, not full unitig-to-contig scaffolding)
- No paired-end / mate-pair linking (scaffolding stub exists but not implemented)
- No use of long reads for gap bridging

**Future**: Integration with external assemblers via wrappers.

---

#### `scaffold_contigs(contigs, paired_reads=None, insert_size=500, min_links=3) -> list[Scaffold]`

**Purpose**: Order and orient contigs using paired-end or mate-pair linkages.

**Algorithm** (stub; not fully implemented):
1. Map paired reads to contigs via k-mer anchoring
2. For each contig pair, count number of read pairs with one in each contig in correct orientation/insert size
3. Build graph: nodes = contigs, edges = linkage counts
4. Find linear chains (paths) with greedy alignment
5. Join contigs with `N` gaps estimated from insert size

**Parameters**:
- `contigs` (list[Contig]): Contig objects from `assemble_contigs()`
- `paired_reads` (dict[str, tuple[str,str]] | None): Paired read dictionaries (R1,R2); if None, uses previously loaded reads
- `insert_size` (int): mean insert length (default 500 bp)
- `min_links` (int): Minimum supporting read pairs to join contigs (default 3)

**Returns**: `list[Scaffold]`

**Status**: Placeholder scaffold implementation; for production use MEGAHIT/metaSPAdes.

---

#### `calculate_assembly_stats(contigs) -> AssemblyStats`

**Purpose**: Compute standard assembly quality metrics.

**Metrics computed**:
- `total_contigs`, `total_length`
- `n50`, `l50`, `n90`, `l90`
- `gc_content` (weighted by contig length)
- `mean_coverage` (weighted)
- `largest_contig`, `smallest_contig`

**Parameters**: `contigs` (list[Contig])

**Returns**: `AssemblyStats`

**Example**:
```python
stats = calculate_assembly_stats(contigs)
print(f"N50 = {stats.n50:,} bp")
```

---

### Binning Functions

#### `calculate_tetranucleotide_freq(sequence, normalize=True) -> list[float]`

**Purpose**: Compute tetranucleotide frequency (TNF) profile for a DNA sequence.

**Algorithm**:
1. Slide 4-nt window across sequence; count occurrences of each 4-mer (256 possibilities: AAAA, AAAC, ..., TTTT)
2. Divide by (len(sequence) - 3) to get frequencies
3. Optionally normalize to zero mean, unit variance

**Parameters**:
- `sequence` (str): DNA string
- `normalize` (bool): If True, subtract mean and divide by std across 256 bins

**Returns**: `list[float]` length 256

**Usage**: TNF is a compositional signature of genomic DNA; different taxa have distinct 4-mer frequencies due to codon bias, GC% pattern.

---

#### `bin_contigs(contigs, coverage=None, method="composition", n_bins=None, min_contig_length=1000) -> BinningResult`

**Purpose**: Group contigs into metagenome-assembled genomes (MAGs) using unsupervised clustering.

**Methods**:

1. **`composition`**: Cluster contigs based solely on TNF vectors (256-D). Uses k-means.
2. **`coverage`**: Cluster based on coverage profiles across multiple samples (if multi-sample co-assembly). Coverage vector per contig = (cov_sample1, cov_sample2, ...).
3. **`combined`**: Concatenate normalized TNF + normalized coverage (recommended). Normalization: z-score per feature (mean=0, std=1).

**Parameters**:
- `contigs` (dict[str, str]): `{contig_id: sequence}`
- `coverage` (dict[str, list[float]] | None): `{contig_id: [cov_sample1, cov_sample2, ...]}`; required for `coverage` or `combined` methods
- `method` (str): `"composition"`, `"coverage"`, or `"combined"`
- `n_bins` (int | None): Number of bins (clusters). If None, estimated from total assembly size / typical bacterial genome size (≈3 Mbp) → guess of strain count
- `min_contig_length` (int): Skip contigs shorter than this for clustering (default 1000) to avoid spurious short contigs

**Returns**: `BinningResult` with `bins` (list[GenomeBin]), `unbinned_contigs` (int), total counts

**Algorithm details**:
1. **Feature extraction**: For each contig ≥ `min_contig_length`, compute feature vector:
   - Composition: 256 TNF bins (normalized)
   - Coverage: per-sample depth values (normalized to median=1)
   - Combined: 256 + n_samples dimensions
2. **Dimensionality reduction** (optional future): PCA to 10–50 PCs to denoise
3. **Clustering**: k-means with k-means++ initialization
4. **Assignment**: Each contig assigned to nearest cluster centroid

**Quality assessment** (CheckM-style, done via separate `refine_bins()`):
- **Completeness**: Fraction of single-copy marker genes found (universal genes like rpoB, gyrB, etc.)
- **Contamination**: Fraction of markers found in >1 copy (mixed strain) or presence of non-universal markers

**Example**:
```python
contigs = {"c1": "ACGT...", "c2": "GGCC...", ...}
coverage = {"c1": [50.0, 55.0], "c2": [52.0, 48.0], ...}  # 2 samples
result = bin_contigs(contigs, coverage, method="combined", n_bins=10)
print(f"Assigned {result.binned_contigs} contigs to {len(result.bins)} bins")

# Refine for quality
from metainformant.metagenomics.shotgun.binning import refine_bins
hq, mq, lq = refine_bins(result.bins, contigs)
print(f"High quality bins: {len(hq)}, Medium: {len(mq)}, Low: {len(lq)}")
```

---

#### `refine_bins(bins, contigs, completeness_threshold=0.5, contamination_threshold=0.1) -> tuple[list[GenomeBin], list[GenomeBin], list[GenomeBin]]`

**Purpose**: Split/merge/filter bins based on marker gene completeness and contamination.

**Classification**:
- **High-quality (HQ)**: completeness ≥ 90%, contamination ≤ 5%
- **Medium-quality (MQ)**: completeness ≥ 50%, contamination ≤ 10%
- **Low-quality (LQ)**: completeness < 50% or contamination > 10%

**Parameters**:
- `bins` (list[GenomeBin]): Bins from `bin_contigs()`
- `contigs` (dict[str, str]): Full contig sequences for marker search
- `completeness_threshold` (float): Minimum completeness to keep bin
- `contamination_threshold` (float): Maximum contamination allowed

**Returns**: `(hq_bins, mq_bins, lq_bins)` — three lists

**Algorithm**:
1. For each bin, search contigs for universal single-copy marker genes (HMMs or consensus patterns)
2. Compute completeness = (# markers found) / (# universal markers)
3. Compute contamination = (# markers found in ≥2 copies) / (# markers found)
4. Filter bins below completeness threshold or above contamination threshold
5. (Future) Split contaminated bins by clustering contigs based on TNF divergence

---

### Profiling Functions

#### `build_kmer_index(reference_sequences, taxonomy=None, k=31) -> KmerIndex`

**Purpose**: Build a k-mer lookup index from reference genomes for taxonomic classification.

**Data model**:
```python
@dataclass
class KmerIndex:
    kmer_to_taxon: dict[str, str]          # k-mer → taxon ID (LCA)
    taxon_lineages: dict[str, list]        # taxon_id → [(rank, name), ...]
    k: int
    total_kmers: int
    num_taxa: int
```

**Algorithm**:
1. For each reference genome:
   - Extract all k-mers (sliding window)
   - Insert into `kmer_to_taxon`:
     - If k-mer not present: map to this taxon
     - If present with another taxon: update to LCA of current and stored taxon
2. Store taxonomy for later reporting

**Parameters**:
- `reference_sequences` (dict[str, str]): `{genome_id: sequence}` (could be species genomes, marker genes)
- `taxonomy` (dict[str, list] | None): `{genome_id: [(rank, name), ...]}`; if None, lineage not stored
- `k` (int): k-mer size, default 31 (good balance specificity/sensitivity)

**Returns**: `KmerIndex`

**Example**:
```python
refs = {"genusA_sp1": "ACGT...", "genusA_sp2": "ACGT..."}
tax = {"genusA_sp1": [("genus","GenusA"), ("species","sp1")], ...}
index = build_kmer_index(refs, tax, k=31)
print(f"Indexed {index.total_kmers} k-mers from {index.num_taxa} taxa")
```

---

#### `profile_community(reads, database=None, reference_sequences=None, k=31, min_kmer_hits=2, confidence_threshold=0.5) -> CommunityProfile`

**Purpose**: Classify reads taxonomically using k-mer LCA approach (Kraken-style).

**Algorithm**:
For each read:
1. Extract all k-mers
2. Look up each k-mer in database; collect matching taxa
3. Compute LCA of all taxonomic hits (most specific common ancestor)
4. Confidence = fraction of k-mers agreeing with LCA (vs. up-assigned to parent)
5. If confidence < `confidence_threshold`, classify as `"unclassified"`

**Parameters**:
- `reads` (list[str] | dict[str, str]): Unclassified sequences
- `database` (KmerIndex | None): Pre-built index from `build_kmer_index()`
- `reference_sequences` (dict[str, str] | None): Alternative: build index on-the-fly from these references
- `k` (int): k-mer size (must match database)
- `min_kmer_hits` (int): Minimum k-mers matching to consider classification (default 2)
- `confidence_threshold` (float): Minimum fraction agreeing to assign (default 0.5)

**Returns**: `CommunityProfile`

**Example**:
```python
reads = ["ACGTACGT...", "TTTTAAAC...", ...]
profile = profile_community(reads, database=index)
for taxon in profile.taxa[:10]:
    print(f"{taxon.name}: {taxon.relative_abundance:.2%} "
          f"({taxon.read_count} reads)")
```

**Output interpretation**:
- Reads with no hits → `unclassified_reads` count
- Relative abundances sum to 1.0 (excluding unclassified) or can include unclassified in total

---

## Diversity

Diversity functions are pure statistical calculations; they do not parse data. Input is abundance vectors (counts or relative abundances).

### Alpha Diversity

#### `alpha_diversity(abundances, metric="shannon") -> dict`

**Purpose**: Compute within-sample diversity.

**Supported metrics**:

| Metric | Input type | Range | Interpretation |
|--------|------------|-------|---------------|
| `shannon` | counts or prop | [0, log(S)] | Entropy; higher = more diverse, accounts for evenness |
| `simpson` | counts or prop | [0, 1] | Probability two random individuals are different species |
| `invsimpson` | counts or prop | [1, S] | Reciprocal Simpson; effective number of species |
| `chao1` | integer counts | [S_obs, ∞) | Richness estimator; S_est = S_obs + (singleton^2 / doubleton) |
| `ace` | integer counts | [S_obs, ∞) | Abundance-based Coverage Estimator; uses ≥10-abundance class |
| `observed` | counts | [1, S_max] | Simple species count (richness) |
| `fisher_alpha` | counts | (0, ∞) | Parameter of Fisher's log-series; higher = more diverse |
| `pielou_evenness` | counts or prop | [0, 1] | J = H / ln(S); 1 = perfectly even |

**Parameters**:
- `abundances` (list[int] | list[float]): Per-species counts (integers) or frequencies
- `metric` (str): One of the above

**Returns**: `dict` with keys:
- `value` (float): Computed diversity
- `metric` (str): metric name
- `n_species` (int): Number of species with non-zero abundance
- `total_count` (float): Sum of abundances

**Example**:
```python
counts = [100, 50, 25, 10, 5]  # 5 species
result = alpha_diversity(counts, metric="shannon")
print(f"Shannon diversity = {result['value']:.3f}")

# Compare metrics
for m in ["shannon", "simpson", "chao1", "observed"]:
    r = alpha_diversity(counts, metric=m)
    print(f"{m}: {r['value']:.3f}")
```

---

### Beta Diversity

#### `beta_diversity(samples, metric="bray_curtis") -> dict`

**Purpose**: Compute pairwise dissimilarity between samples.

**Supported metrics**:

| Metric | Data type | Range | Notes |
|--------|-----------|-------|-------|
| `bray_curtis` | counts/abundances | [0, 1] | Quantitative; 1 = no shared abundance |
| `jaccard` | presence/absence | [0, 1] | Qualitative; shared species / total species |
| `aitchison` | compositional data | [0, √2] | Euclidean on CLR-transformed abundances |

**Parameters**:
- `samples` (list[list[float]]): Outer list = samples; inner = per-species abundances (same length across samples)
- `metric` (str): `"bray_curtis"`, `"jaccard"`, or `"aitchison"`

**Returns**: `dict` with:
- `distance_matrix` (list[list[float]]): Symmetric n×n matrix
- `metric` (str)
- `n_samples` (int)

**Example**:
```python
samples = [
    [100, 50, 0, 10],    # sample A
    [10, 5, 0, 50],      # sample B
    [0, 0, 100, 0],      # sample C
]
result = beta_diversity(samples, metric="bray_curtis")
print(result["distance_matrix"][0][1])  # A vs B dissimilarity
```

**Note**: For `aitchison`, uses CLR (centered log-ratio) transform prior to Euclidean; handles zeros via small pseudocount (1e-6).

---

### Ordination

#### `ordination(distance_matrix, method="pcoa", n_components=2) -> dict`

**Purpose**: Dimensionality reduction of distance matrix for visualization.

**Methods**:
- **PCoA** (Principal Coordinates Analysis): Eigen-decomposition of distance matrix (classical MDS). Preserves distances as much as possible in low-D.
- **NMDS** (Non-metric Multidimensional Scaling): Optimizes rank ordering; focuses on ordinal relationships, not metric distances.

**Parameters**:
- `distance_matrix` (list[list[float]]): Square symmetric matrix
- `method` (str): `"pcoa"` (default) or `"nmds"`
- `n_components` (int): Target dimensions (2 or 3 for plotting)

**Returns**: `dict` with:
- `coordinates` (list[list[float]]): Sample coordinates (n × n_components)
- `eigenvalues` (list[float]): Eigenvalues (PCoA) or stress values (NMDS) per component
- `variance_explained` (list[float]): Fraction of total variance per component (PCoA only)

**Example**:
```python
pcoa = ordination(dist_matrix, method="pcoa", n_components=2)
xs = [coord[0] for coord in pcoa["coordinates"]]
ys = [coord[1] for coord in pcoa["coordinates"]]
# Plot: plt.scatter(xs, ys)
```

**Note**: NMDS not fully stress-optimized in current implementation; uses simple scikit-learn `MDS` with non-metric=True.

---

### Rarefaction

#### `rarefaction_curve(abundances, depths=None, n_iterations=10, seed=None) -> dict`

**Purpose**: Generate rarefaction curve to assess sampling depth sufficiency.

**Algorithm**:
1. Build pool of individual organisms: `[taxon_idx repeated count times]`
2. For each sampling depth:
   - Randomly subsample `depth` individuals without replacement (repeat `n_iterations` times)
   - Record mean observed species across iterations ± std
3. Check saturation: Last 3 points within 5% → curve considered saturated

**Parameters**:
- `abundances` (list[int]): Counts per species
- `depths` (list[int] | None): Specific depths to evaluate; if None, auto-generate 20 evenly spaced points from 1 to total
- `n_iterations` (int): Random subsampling replicates (default 10; increase to 100 for smooth curves)
- `seed` (int | None): Random seed for reproducibility

**Returns**: `dict`:
- `depths` (list[int]): Depths evaluated
- `mean_species` (list[float]): Mean observed species at each depth
- `std_species` (list[float]): Standard deviation
- `is_saturated` (bool): Has plateaued?

**Example**:
```python
abundances = [1000, 500, 200, 50, 10, 5, 1, 1]  # community with many rare taxa
curve = rarefaction_curve(abundances, n_iterations=50)
plt.errorbar(curve["depths"], curve["mean_species"], yerr=curve["std_species"])
plt.xlabel("Sequencing depth (reads)")
plt.ylabel("Observed ASVs")
```

**Interpretation**: If curve is still rising steeply at max depth, undersampled; consider deeper sequencing.

---

#### `rarefy(abundances, depth, seed=None) -> list[int]`

**Purpose**: Subsample community to even depth (rarefaction).

**Parameters**:
- `abundances` (list[int]): Counts per taxon
- `depth` (int): Target total reads after subsampling
- `seed` (int | None): For reproducibility

**Returns**: `list[int]` rarefied counts (same length as input, sum = `depth`)

**Example**:
```python
counts = [1000, 500, 200]
rarefied = rarefy(counts, depth=1000)
print(f"Original: {sum(counts)}, Rarefied: {sum(rarefied)}")
```

**Use case**: Standardize to common depth before beta diversity to correct for library size differences.

---

### PERMANOVA

#### `permanova(distance_matrix, groups, n_permutations=999, seed=None) -> dict`

**Purpose**: Permutational multivariate analysis of variance — test whether group dissimilarities are greater than expected by chance.

**Hypothesis**:
- H₀: No difference in community composition between groups
- H₁: Groups differ in multivariate space

**Test statistic**: Pseudo-F = (SS_between / df_between) / (SS_within / df_within)

**Permutation test**: Shuffle group labels `n_permutations` times, compute pseudo-F each time; p-value = (1 + #permuted ≥ observed) / (n_permutations + 1)

**Parameters**:
- `distance_matrix` (list[list[float]]): Square matrix
- `groups` (list[str]): Group label per sample (e.g., `["control","control","treat","treat"]`)
- `n_permutations` (int): Number of permutations (default 999; use 9999 for publication)
- `seed` (int | None): Random seed

**Returns**: `dict`:
- `pseudo_f` (float)
- `p_value` (float)
- `r_squared` (float) — fraction of variance explained by groups
- `n_permutations` (int)

**Example**:
```python
beta = beta_diversity(samples, metric="bray_curtis")
permanova_result = permanova(
    beta["distance_matrix"],
    groups=["A","A","B","B","C","C"],
    n_permutations=9999,
)
print(f"Pseudo-F = {permanova_result['pseudo_f']:.3f}, "
      f"p = {permanova_result['p_value']:.4g}, R² = {permanova_result['r_squared']:.3f}")
```

**Interpretation**: p < 0.05 → reject H₀; groups have significantly different community composition.

---

## Functional

### Gene Annotation

**Status**: Minimal implementations. Details TBD.

#### `predict_genes(contigs) -> list[Gene]`

Predict open reading frames (ORFs) from contigs using Glimmer-style interpolated Markov models or simple six-frame translation with long ORF heuristic.

**Parameters**: `contigs` (list[Contig])

**Returns**: `list[Gene]` with coordinates, strand, amino acid sequence

---

#### `annotate_genes(genes, database="eggNOG", ...) -> list[AnnotatedGene]`

Assign functional annotations (COG, KEGG, PFAM) to predicted genes via homology search (DIAMOND/BLAST) or HMM (HMMER).

**Parameters**: `genes`, `database`

**Returns**: `list[AnnotatedGene]` with gene ontology, enzyme commission numbers, pathway IDs

---

### Pathway Analysis

#### `reconstruct_pathways(annotated_genes, method="minimal") -> dict[str, float]`

Compute pathway abundances by counting pathway-specific enzymes or reactions.

**Parameters**:
- `annotated_genes` (list[AnnotatedGene])
- `method` (str): `"minimal"` (presence/absence), `"abundance"` (sum gene abundances), `"coverage"` (pathway completeness %)

**Returns**: `dict[str, float]` pathway_name → abundance/completeness score

---

#### `calculate_completeness(pathway_abundances, reference="KEGG") -> list[PathwayCompleteness]`

Assess how complete each pathway is relative to reference genome (0–100%).

---

## Comparative

#### `differential_abundance(feature_table, group_a, group_b, method="aldex2") -> dict`

**Purpose**: Identify features (OTUs, ASVs, taxa, genes, pathways) that differ in abundance between two groups.

**Supported methods**:

| method | Statistics | Compositional-aware? | Recommended |
|--------|------------|---------------------|-------------|
| `"aldex2"` | CLR-transformed Welch's t-test + Monte Carlo | Yes | 16S default |
| `"ancom"` | ANCOM test (bias-corrected) | Yes | 16S compositional |
| `"wilcoxon"` | Rank-sum test | No | Nonparametric, small n |
| `"t_test"` | Standard t-test | No | Absolute abundances only |

**Parameters**:
- `feature_table`: 2D array-like (features × samples) or (samples × features) depending on convention
- `group_a`, `group_b`: Sample indices or boolean masks
- `method`: statistical method

**Returns**: `dict[int | str, dict]` mapping feature index/ID to statistics:
- `log_fc` (float): Log2 fold change (A vs B)
- `p_value` (float)
- `q_value` (float): BH-FDR corrected
- `significant` (bool): q < 0.05

**Example**:
```python
diff = differential_abundance(
    feature_table=otu_counts.T,  # samples × OTUs
    group_a=[0,1,2],
    group_b=[3,4,5],
    method="aldex2",
)
sig_features = [i for i, r in diff.items() if r["significant"]]
```

---

## Visualization

**Current**: Basic plot stubs. For production plots, use `metainformant.visualization` instead.

### Exported Functions (in `visualization/plots.py`)

| Function | Description |
|----------|-------------|
| `stacked_bar(abundances, top_n=20, ...)` | Stacked bar chart of top taxa |
| `krona_chart(profile)` | Interactive Krona radiant (circular hierarchy) |
| `ordination_plot(coords, groups)` | PCoA/NMDS scatter plot with ellipses |
| `heatmap(data, row_labels, col_labels)` | Heatmap with dendrograms |
| `rarefaction_curve(curve_data)` | Plot rarefaction curves |

---

## Algorithm Comparison

### Clustering: OTU vs. ASV

| | OTU Clustering | ASV Denoising |
|---|---------------|---------------|
| **Resolution** | 97% identity clusters (OTU) | Exact sequence variants (100%) |
| **Speed** | Fast (O(n) with prefilter) | Moderate (error modeling) |
| **Sensitivity** | Lumps true variants | Resolves single-nucleotide differences |
| **Error handling** | Noise filtered by clustering threshold | Explicit error model (Poisson) |
| **Reproducibility** | Depends on order of input, threshold | Deterministic given same data |
| **Best for** | Large studies, legacy comparison | Strain tracking, rare biosphere |


### Taxonomic Classification: Naive Bayes vs. BLAST

| | Naive Bayes (RDP) | BLAST (k-mer LCA) |
|---|-------------------|-------------------|
| **Speed** | Fast after training | Slower (linear scan) |
| **Training** | Requires labeled reference DB | No training (direct lookup) |
| **Confidence** | Bootstrap support per rank | LCA depth, no per-rank confidence |
| **Memory** | Small (k-mer probs) | Large (store all ref kmers) |
| **Best for** | Large reference sets, fast batch | Small references, one-off queries |

---

### Differential Abundance: Compositional vs. Standard

Compositional data (relative abundances) have sum constraint → spurious correlations. Use compositional-aware tests:

- **ALDEx2**: CLR transform + Monte-Carlo Dirichlet sampling
- **ANCOM**: Tests each feature against the geometric mean of all others (log-ratio)
- **Standard t-test**: Only valid if using absolute counts from spike-in normalization

---

## Function Matrix

| Submodule | Function | Class | Output |
|-----------|----------|-------|--------|
| amplicon | `cluster_otus` | `ClusteringResult` | OTU table |
| amplicon | `calculate_identity` | function | float |
| amplicon | `filter_chimeras` | dict[str, bool] | chimera flags |
| amplicon | `denoise_sequences` | `DenoisingResult` | ASV table |
| amplicon | `merge_paired_reads` | dict[str, str] | merged reads |
| amplicon | `classify_taxonomy` | list[TaxonomyAssignment] | classified sequences |
| amplicon | `build_taxonomy_tree` | TaxonomyNode | tree |
| shotgun | `assemble_contigs` | list[Contig] | assembly |
| shotgun | `scaffold_contigs` | list[Scaffold] | scaffolds |
| shotgun | `calculate_assembly_stats` | AssemblyStats | metrics |
| shotgun | `bin_contigs` | BinningResult | bins |
| shotgun | `refine_bins` | (hq,mq,lq) | quality bins |
| shotgun | `build_kmer_index` | KmerIndex | index |
| shotgun | `profile_community` | CommunityProfile | taxon abundances |
| diversity | `alpha_diversity` | dict | alpha value |
| diversity | `beta_diversity` | dict | distance matrix |
| diversity | `ordination` | dict | coords |
| diversity | `rarefaction_curve` | dict | curve data |
| diversity | `rarefy` | list[int] | subsampled counts |
| diversity | `permanova` | dict | test stats |
| functional | `predict_genes` | list[Gene] | (stub) |
| functional | `annotate_genes` | list[AnnotatedGene] | (stub) |
| functional | `reconstruct_pathways` | dict[str, float] | (stub) |
| functional | `calculate_completeness` | list[PathwayCompleteness] | (stub) |
| comparative | `differential_abundance` | dict[feature → {logfc,p,q}] | stats |

---

## Configuration

| Variable | Purpose | Default |
|----------|---------|--------|
| `META_OTU_THRESHOLD` | Default OTU identity threshold for `cluster_otus()` | `0.97` |
| `META_` prefix (future) | General module overrides via `core.config` | N/A |

All other parameters are function arguments.

---

## Performance Tips

- **OTU clustering**: Use `prefilter=True` (default) for large datasets; reduces comparisons to near-linear.
- **Taxonomic classification**: Build KmerIndex once, reuse across many samples.
- **Diversity metrics**: All O(n_features × n_samples) or better; trivial for typical sizes (<10k features, <1000 samples).
- **Assembly**: Most expensive operation; consider running external assemblers (MEGAHIT) on HPC and only loading contigs into METAINFORMANT for downstream.

---

## Limitations

- **Assembly**: Not production-grade; use MEGAHIT/MetaSPAdes for real data
- **Binning**: Basic k-means; CheckM not integrated yet
- **Annotation**: Stubs only; no ORF prediction or annotation implemented
- **Profiling**: k-mer LCA simple; Kraken2/Bracken faster and more accurate

The module provides a solid foundation and algorithmic reference; for production, integrate external tools where indicated.
