# Metagenomics Module Capabilities

Complete function and class reference for the metagenomics module, organized by subpackage. The `metainformant.metagenomics` package exposes its six subpackages only, so import each function from the module that defines it (for example `metainformant.metagenomics.amplicon.otu_clustering`).

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
from metainformant.metagenomics.amplicon.otu_clustering import OTU

otu = OTU(
    centroid_id="OTU_1",
    centroid_sequence="ACGTACGTACGT",
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
from metainformant.metagenomics.amplicon.otu_clustering import cluster_otus

sequences = {"s1": "ACGTACGTACGT", "s2": "ACGTACGTACGA", "s3": "TTTTAAAACCCC"}
result = cluster_otus(sequences, threshold=0.90)
print(f"{result.num_otus} OTUs from {result.total_sequences} sequences")
for otu in result.otus[:5]:
    print(f"  {otu.centroid_id}: size={otu.size}")
```

---

#### `ErrorModel`, `ASV` and `DenoisingResult`

**Dataclasses** from `asv_denoising`.

- `ErrorModel`: `transition_rates` (4x4 list, rows and columns A, C, G, T; entry `[i][j]` is P(observe j | true base i); one matrix for all positions), `quality_error_rates` (dict, Phred 0–41 → `10 ** (-Q / 10)`), `positions_modeled` (int, longest read length)
- `ASV`: `sequence`, `abundance` (reads assigned, including merged variants), `sample_counts` (dict, left empty by `denoise_sequences`), `quality` (mean Phred score over all positions of the reads assigned to the ASV; 0.0 when no quality scores were supplied)
- `DenoisingResult`: `asvs` (list[ASV]), `error_model`, `total_reads`, `num_asvs`, `reads_removed` (reads merged into other ASVs)

---

#### `TaxonomyAssignment`

**Dataclass**: Taxonomic classification of a sequence.

**Fields**:
- `sequence_id` (str): Identifier of classified sequence
- `lineage` (list[tuple[str, str]]): One `(rank, name)` tuple per rank in `TAXONOMIC_RANKS` (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`), e.g., `[("domain","Bacteria"), ("phylum","Firmicutes"), ("class","unclassified"), ...]`
- `confidence` (dict[str, float]): Per-rank confidence scores (0–1)
- `method` (str): Classification method used (`"naive_bayes"`, `"blast"`)
- `best_hit_id` (str | None): Best-matching reference sequence ID (`"blast"` only)
- `best_hit_identity` (float): 11-mer Jaccard score of the best hit, 0–1 (`"blast"` only; not a sequence identity)

---

#### `TaxonomyNode`

**Dataclass**: Node of the tree built by `build_taxonomy_tree()`.

**Fields**: `name`, `rank`, `count` (sequences ending at this node), `children` (dict of name → node), `parent`. **Methods**: `add_child(name, rank)`, `get_lineage()`, `total_descendants()`, `to_dict()`.

---

### Functions

#### `cluster_otus(sequences, threshold=None, abundance=None, sort_by="length", prefilter=True) -> ClusteringResult`

**Purpose**: Greedy centroid-based OTU clustering (VSEARCH/UCLUST style).

**Algorithm**:
1. Sort sequences by abundance (descending) or length (descending)
2. First sequence → centroid of OTU #1
3. For each subsequent sequence:
   - Compare to all existing centroids using `calculate_identity()`
   - If identity ≥ threshold for any centroid, assign to that OTU (best match)
   - Else, create new OTU with this sequence as centroid

**Pre-filter**: A conservative k-mer count test (the q-gram lemma on 8-mers) skips an alignment only when it proves that the pair cannot reach the threshold: sequences that differ in at most `e` alignment columns share at least `longer length − 8 + 1 − 8e` 8-mers (repeated 8-mers counted), and `e` is bounded by the threshold and the shorter length. A qualifying pair is never skipped, so the clustering is identical with and without the pre-filter.

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: nucleotide_sequence}`
- `threshold` (float | None): Identity threshold in (0, 1]. An explicit value always wins. If `None`, the `META_OTU_THRESHOLD` environment variable is used when set (it is validated: `ValueError` for anything but a number in (0, 1]; a blank value counts as unset), otherwise 0.97
- `abundance` (dict[str, int] | None): Precomputed sequence abundances for sorting; if None and `sort_by="abundance"`, falls back to length sorting
- `sort_by` (str): `"length"` (default) or `"abundance"` centroid seeding order
- `prefilter` (bool): Skip alignments that provably cannot reach the threshold (default True; speed only, same result)

**Returns**: `ClusteringResult`

**Raises**: `ValueError` for empty input or an invalid threshold.

**Complexity**: About one O(L₁×L₂) pure-Python alignment per sequence–centroid pair that passes the pre-filter.

**Example**:
```python
from metainformant.metagenomics.amplicon.otu_clustering import cluster_otus

seqs = {
    "s1": "ACGTACGTACGT",
    "s2": "ACGTACGTACGA",  # 1 mismatch
    "s3": "TTTTAAAACCCC",
}
result = cluster_otus(seqs, threshold=0.90)
print(result.num_otus)  # 2: {s1, s2} cluster, s3 separate
```

**When to use**: Legacy 16S analysis where 97% OTUs are desired; hundreds to a few thousand sequences (the alignment step is pure Python).

---

#### `calculate_identity(seq1, seq2) -> float`

**Purpose**: Compute pairwise global sequence identity using Needleman-Wunsch with affine gap penalties.

**Scoring**:
- Match: +2
- Mismatch: -1
- Gap open: -5
- Gap extend: -1

**Identity definition**: `(# identical positions in alignment) / (alignment length, gap columns included)`.

**Parameters**:
- `seq1`, `seq2` (str): Nucleotide strings (upper-cased internally)

**Returns**: `float` in [0.0, 1.0]

**Example**:
```python
from metainformant.metagenomics.amplicon.otu_clustering import calculate_identity

identity = calculate_identity("ACGTACGT", "ACGTACGA")  # 7/8 = 0.875
```

**Notes**:
- Removes `-` characters before alignment
- Raises `ValueError` if either sequence is empty
- Full dynamic-programming matrices (no banding); identity is counted during the traceback

**When to use**: Called internally by `cluster_otus()`; exported for direct pairwise comparisons.

---

#### `filter_chimeras(sequences, reference_db=None, abundance=None, min_divergence=1.0, min_score=0.28) -> dict[str, bool]`

**Purpose**: Detect chimeric sequences with a simplified UCHIME-style de novo test.

**Algorithm**:
1. Sort sequences by abundance (most abundant first; descending length when `abundance` is None)
2. Without `reference_db`, the first two sequences are assumed non-chimeric. For each later query:
   - Find the two best-matching parent candidates among the non-chimeric sequences accepted so far (or among `reference_db`) using the shared 8-mer fraction (`_quick_identity()`)
   - For breakpoints at 10% intervals, build chimeric models from the two parents:
     - Model A: left from parent A, right from parent B
     - Model B: left from parent B, right from parent A
   - A model counts only if its position-wise identity to the query exceeds that of the closer parent by at least `min_divergence` percentage points
   - Chimera score = (model identity − closer parent identity) / closer parent identity; the query is flagged if the best score is ≥ `min_score`

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: sequence}`
- `reference_db` (dict[str, str] | None): If provided, use these as the parent pool instead of the abundance-sorted input (reference-based chimera check)
- `abundance` (dict[str, int] | None): Abundance counts for sorting; if None, uses length
- `min_divergence` (float): Minimum margin of the chimeric model over the closer parent, in percentage points of identity (the UCHIME `mindiv` convention, so 1.0 means 1%; default 1.0). A negative or NaN value raises `ValueError`
- `min_score` (float): Threshold to flag as chimeric (0.0–1.0, default 0.28)

**Returns**: `dict[str, bool]` — `True` = chimeric, `False` = non-chimeric

**Example**:
```python
from metainformant.metagenomics.amplicon.otu_clustering import filter_chimeras

seqs = {
    "parent1": "AAAAAAAAAA",
    "parent2": "TTTTTTTTTT",
    "chimera": "AAAAATTTTT",  # hybrid of parent1 + parent2
}
result = filter_chimeras(seqs)
print(result["chimera"])  # True
print(result["parent1"])  # False (one of the first two sequences → assumed good)
```

**Note**: Chimera detection is approximate; sequences shorter than 10 nt are never flagged. Inspect flagged sequences manually in production.

---

#### `estimate_error_rates(quality_scores, sequences=None) -> ErrorModel`

**Purpose**: Learn a base-substitution error model from Phred quality scores.

**Parameters**:
- `quality_scores` (list[list[float]]): One list of Phred scores per read
- `sequences` (list[str] | None): One sequence per read; when given, the relative frequencies of the substitutions observed against the per-position consensus decide which base replaces the true one (otherwise the three alternatives are equally likely)

**Method**: The error probability of a base call is `10 ** (-Q / 10)`, bounded to [1e-12, 0.75], and it is averaged as a probability (which exceeds the probability at the mean Q whenever the scores vary). Each row of `transition_rates` is `[1 - error, error * share, ...]` and sums to 1.

**Returns**: `ErrorModel`. **Raises**: `ValueError` for empty `quality_scores`.

**Example**:
```python
from metainformant.metagenomics.amplicon.asv_denoising import estimate_error_rates

model = estimate_error_rates([[30, 30, 25, 20], [28, 30, 22, 18]])
print(1 - model.transition_rates[0][0])  # about 0.0050 (the probability at the mean Q would be 0.0029)
```

---

#### `denoise_sequences(sequences, error_rates=None, quality_scores=None, abundance=None, omega_a=1e-40, min_abundance=1, band_size=16) -> DenoisingResult`

**Purpose**: Simplified DADA2-style denoising to resolve Amplicon Sequence Variants (ASVs) from amplicon reads, modeling sequencing errors.

**Algorithm summary** (not full DADA2):
1. **Dereplication**: Count unique sequences and abundances (or take counts from `abundance`); drop those below `min_abundance`
2. **Error model**: use `error_rates`, else estimate one from `quality_scores`, else a default Phred-30 model (error probability 0.001)
3. **Merging error variants**: visit sequences by descending abundance. For each, compute `lambda = (ASV abundance) × P(variant | ASV)` against every earlier ASV whose length differs by at most `band_size` (no alignment) and the abundance p-value `P(X >= abundance)` for `X ~ Poisson(lambda)`, evaluated as `poisson.sf(abundance - 1, lambda)`. Merge into the first ASV with p > `omega_a`; otherwise start a new ASV
4. **Output**: ASVs with final abundances and mean Phred quality

**Parameters**:
- `sequences` (dict[str, str]): Read sequences (keys arbitrary, values DNA strings)
- `error_rates` (ErrorModel | None): Pre-computed error model
- `quality_scores` (dict[str, list[int]] | None): `{read_id: [Phred score per position]}`. They set the error probability of every mismatch (the mean of `10 ** (-Q / 10)` over the reads that share a sequence), feed the error-model estimate, and give each ASV its `quality`
- `abundance` (dict[str, int] | None): Per-sequence counts; if None, duplicates in the input are counted
- `omega_a` (float): P-value threshold for merging (default 1e-40; a smaller value merges more variants into their parents)
- `min_abundance` (int): Minimum abundance for a sequence to be considered (default 1)
- `band_size` (int): Maximum length difference between a variant and a parent for them to be compared (default 16)

**Returns**: `DenoisingResult` with `asvs`, `error_model`, `total_reads`, `num_asvs`, `reads_removed`

**Example**:
```python
from metainformant.metagenomics.amplicon.asv_denoising import denoise_sequences

reads = {"r1": "ACGTACGT", "r2": "ACGTACGA", "r3": "ACGTACGA"}
result = denoise_sequences(reads)
print(f"Input: {len(reads)} reads → {result.num_asvs} ASVs")
for asv in result.asvs:
    print(f"  {asv.sequence}: {asv.abundance} reads")
```

**Recommendation**: Use ASV denoising for high-resolution microbiome analysis (e.g., strain-level tracking, single-nucleotide variants). Requires higher sequencing depth than OTU clustering.

---

#### `merge_paired_reads(forward, reverse, min_overlap=20, max_mismatch_ratio=0.2, quality_forward=None, quality_reverse=None) -> dict[str, str]`

**Purpose**: Merge overlapping paired-end reads into a single consensus sequence.

**Algorithm**:
1. Reverse-complement each reverse read internally
2. Find the largest overlap ≥ `min_overlap` whose mismatch ratio is ≤ `max_mismatch_ratio`
3. At mismatching overlap positions, keep the base with the higher Phred score when both quality dicts are given; otherwise keep the forward base
4. Concatenate the non-overlapping parts

**Parameters**:
- `forward` (dict[str, str]): `{id: sequence}` (R1 reads)
- `reverse` (dict[str, str]): `{id: sequence}` (R2 reads as sequenced; they are reverse-complemented internally)
- `min_overlap` (int): Minimum overlapping bases required (default 20)
- `max_mismatch_ratio` (float): Maximum allowed mismatches in overlap (default 0.2 = 20%)
- `quality_forward`, `quality_reverse` (dict[str, list[int]] | None): Phred scores per read

**Returns**: `dict[str, str]` merged sequences keyed by read ID (pairs that do not merge are omitted). Raises `ValueError` if the two dicts share no read ID.

**Example**:
```python
from metainformant.metagenomics.amplicon.asv_denoising import merge_paired_reads

# A 62 bp amplicon read as a 40 bp forward read and a 42 bp reverse read (20 bp overlap)
fwd = {"pair1": "ATCGGCTAAGCTTGCAGGATCCTAGGCTAACGTTAGCCGA"}
rev = {"pair1": "TAGCTAGCCGATCGATCCGTAATCGGCTAACGTTAGCCTAGG"}
merged = merge_paired_reads(fwd, rev, min_overlap=20)
print(len(merged["pair1"]))  # 62
```

---

#### `classify_taxonomy(sequences, reference_db, reference_taxonomy=None, method="naive_bayes", confidence_threshold=0.8, k=8, bootstrap_n=100, seed=42) -> list[TaxonomyAssignment]`

**Purpose**: Assign taxonomic lineage to sequences using a reference database.

**Methods**:

1. **`naive_bayes`** (RDP-style):
   - Extract k-mers (default k=8) from the query sequence
   - Score each taxon at each rank with a prior from the reference k-mer totals and Laplace-smoothed k-mer likelihoods (naive independence assumption); k-mers absent from every reference carry no information and are ignored
   - Confidence = fraction of `bootstrap_n` trials, each resampling `ceil(n_kmers / k)` k-mers with replacement, that reproduce the best taxon at that rank (a trial without a known k-mer, or with a tie, counts as disagreeing)
   - Abstain: a rank whose confidence is below `confidence_threshold` is `"unclassified"`, and so is every rank below it. A query with no known k-mer, or whose best two taxa tie, cannot be classified

2. **`blast`** (k-mer similarity vote, not BLAST):
   - Score every reference by 11-mer Jaccard similarity
   - Take a score-weighted vote among the top 5 hits at each rank; confidence is the vote share

**Parameters**:
- `sequences` (dict[str, str]): `{seq_id: sequence}`
- `reference_db` (dict[str, str]): reference database `{ref_id: sequence}`
- `reference_taxonomy` (dict[str, list] | None): `{ref_id: [(rank, name), ...]}`; if None, each lineage is parsed from the reference ID split on `;`
- `method` (str): `"naive_bayes"` (default) or `"blast"`
- `confidence_threshold` (float): Minimum per-rank confidence for assignment (default 0.8). Ranks below the threshold get `"unclassified"`
- `k` (int): k-mer size for `naive_bayes` (default 8)
- `bootstrap_n` (int): Number of bootstrap trials for `naive_bayes` (default 100)
- `seed` (int | None): Seed of the `naive_bayes` bootstrap (default 42). Every query uses its own generator started from this seed, so results are reproducible and independent of the other queries; `None` draws fresh entropy

**Returns**: `list[TaxonomyAssignment]` one per input sequence

**Example**:
```python
from metainformant.metagenomics.amplicon.taxonomy import classify_taxonomy

ref_db = {"r1": "ATCGATCGATCG", "r2": "GGGGCCCCTTTT"}
ref_tax = {
    "r1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
    "r2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")],
}
assignments = classify_taxonomy(
    {"query1": "ATCGATCGATCG"},
    ref_db, ref_tax,
    method="naive_bayes",
    confidence_threshold=0.7,
)
print(assignments[0].lineage[:2])
# [('domain', 'Bacteria'), ('phylum', 'Firmicutes')]
```

**Performance**: Neither method caches reference data between queries: `naive_bayes` rebuilds the per-rank k-mer tables from the whole reference for every query sequence, and `blast` computes an 11-mer set for every reference for every query. The cost therefore grows with queries × reference size; keep reference sets small.

---

#### `build_taxonomy_tree(classifications) -> TaxonomyNode`

**Purpose**: Build a hierarchical tree from classification results. The root is a `TaxonomyNode` named `root` (rank `root`); a lineage stops at its first `unclassified` rank, which becomes an `unclassified_<rank>` child. Each sequence is counted once, at the deepest node it reaches, so `root.total_descendants()` equals the number of classifications.

---

#### `calculate_confidence(assignments, min_confidence=0.0) -> dict[str, dict[str, float]]`

**Purpose**: Summarize per-rank confidence. Returns `{rank: {"mean", "median", "fraction_confident", "n"}}`, where `fraction_confident` is the fraction of assignments with confidence ≥ `min_confidence`.

---

## Shotgun

### Assembly Data Structures

#### `Contig`

**Fields**:
- `contig_id` (str): Unique identifier (e.g., `"contig_000000"`)
- `sequence` (str): Nucleotide sequence (capital ACGT)
- `length` (int): Sequence length in bp
- `coverage` (float): Mean k-mer count along the contig
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
- `gc_content` (float) — GC fraction over all contigs
- `mean_coverage` (float) — mean over the contigs with non-zero coverage
- `largest_contig` (int)
- `smallest_contig` (int)
- `mean_length` (float), `median_length` (float)

**Calculation**:
Sort contigs by length descending; accumulate lengths until reaching 50% of total; N50 = length of last contig added; L50 = count of contigs used.

---

#### `Scaffold`

**Fields**: `scaffold_id`, `contigs` (list[Contig], longest first), `gaps` (list[int], N-run lengths between consecutive contigs), `sequence` (with N gaps), `total_length`.

---

### Assembly Functions

#### `assemble_contigs(reads, k_range=None, min_contig_length=200, min_kmer_coverage=2) -> list[Contig]`

**Purpose**: Simplified de Bruijn graph de novo assembly.

**Algorithm**:
1. For each k in `k_range`:
   - Count the k-mers of all reads (forward strand as given) and keep those seen at least `min_kmer_coverage` times
   - Build the de Bruijn graph: nodes = (k-1)-mers, edges = k-mers
   - Extract contigs as non-branching paths (no tip removal or bubble popping)
2. Merge contigs across k values: contigs are processed longest first and any contig that is a substring of one already kept is dropped

**Parameters**:
- `reads` (dict[str, str] | list[str]): Read sequences
- `k_range` (list[int] | None): k-mer sizes (each at least 11; even values log a warning). Default depends on the mean read length: `[21]` below 100 bp, `[21, 33]` below 200 bp, `[21, 33, 55]` otherwise
- `min_contig_length` (int): Discard contigs shorter than this (default 200)
- `min_kmer_coverage` (int): Minimum k-mer count to keep an edge (default 2)

**Returns**: `list[Contig]`, longest first, with IDs `contig_000000`, `contig_000001`, ... The output is deterministic: graph start nodes are visited in sorted order, so the same reads and parameters give the same contigs and IDs in every process.

**Complexity**: K-mer counting is O(total bases). All distinct k-mers are held in an in-memory Python `Counter` before filtering, so memory grows with the number of distinct k-mers and sequencing errors add many singletons.

**Limitations**:
- No tip removal, bubble popping or repeat resolution (unitigs only)
- Paired-end information is not used during assembly (see `scaffold_contigs`)
- No use of long reads for gap bridging

---

#### `scaffold_contigs(contigs, paired_reads=None, insert_size=500, insert_std=100, min_links=3) -> list[Scaffold]`

**Purpose**: Group contigs with paired-end linkages.

**Algorithm**:
1. Index the 31-mers of both strands of every contig
2. Map each read of a pair to the contigs sharing a 31-mer with it (a mate sequenced from the opposite strand still maps)
3. For each pair of contigs, count the read pairs with one read in each contig
4. Join contigs whose count is ≥ `min_links` (union-find, strongest links first)
5. List the joined contigs longest first and join them with runs of `N`

**Parameters**:
- `contigs` (list[Contig]): Contig objects from `assemble_contigs()`
- `paired_reads` (dict[str, tuple[str, str]] | None): `{pair_id: (forward, reverse)}`; if None, each contig becomes its own scaffold
- `insert_size` (int): mean insert length (default 500 bp); a gap is `max(10, insert_size - min(length of the preceding contig, insert_size))` N's
- `insert_std` (int): accepted but not used
- `min_links` (int): Minimum supporting read pairs to join contigs (default 3)

**Returns**: `list[Scaffold]`, longest first

**Limitations**: Contigs are joined in the orientation they already have; relative orientation is not inferred. For production use MEGAHIT/metaSPAdes.

---

#### `calculate_assembly_stats(contigs) -> AssemblyStats`

**Purpose**: Compute standard assembly quality metrics.

**Metrics computed**:
- `total_contigs`, `total_length`
- `n50`, `l50`, `n90`, `l90`
- `gc_content` (over all contigs)
- `mean_coverage` (over contigs with non-zero coverage)
- `largest_contig`, `smallest_contig`, `mean_length`, `median_length`

**Parameters**: `contigs` (list[Contig])

**Returns**: `AssemblyStats`. **Raises**: `ValueError` for an empty list.

**Example**:
```python
from metainformant.metagenomics.shotgun.assembly import Contig, calculate_assembly_stats

stats = calculate_assembly_stats([Contig("c1", "A" * 1000), Contig("c2", "T" * 500)])
print(f"N50 = {stats.n50:,} bp")  # 1,000 bp
```

---

### Binning Data Structures

#### `GenomeBin`

**Fields**:
- `bin_id` (str), `contig_ids` (list[str]), `total_length` (int), `num_contigs` (int), `gc_content` (float), `mean_coverage` (float)
- `completeness` (float | None): expected single-copy markers present ÷ expected markers; `None` until assessed
- `contamination` (float | None): extra marker copies ÷ expected markers (CheckM definition; can exceed 1.0); `None` until assessed
- `quality_score` (float | None): completeness − 5 × contamination; `None` until assessed
- `marker_counts` (dict[str, int]): copy number of each expected marker found in the bin

#### `BinningResult`

**Fields**: `bins` (list[GenomeBin]), `unbinned_contigs` (list[str]), `total_contigs`, `binned_contigs`, `high_quality_bins` (int | None: bins with completeness ≥ 0.9 and contamination ≤ 0.05; `None` when `bin_contigs` was not given `marker_genes`, 0 when there are no bins).

#### `BACTERIAL_SINGLE_COPY_MARKERS`

The default expected marker set: a tuple of 36 bacterial single-copy gene symbols (`rpoB`, `rpoC`, `rpsB`, ..., `ychF`).

---

### Binning Functions

#### `calculate_tetranucleotide_freq(sequence, normalize=True) -> list[float]`

**Purpose**: Compute tetranucleotide frequency (TNF) profile for a DNA sequence.

**Algorithm**:
1. Slide a 4-nt window across the sequence and its reverse complement; count each of the 256 possible 4-mers (AAAA, AAAC, ..., TTTT)
2. If `normalize`, divide the counts by their sum so the frequencies add up to 1.0

**Parameters**:
- `sequence` (str): DNA string (at least 4 bases, otherwise `ValueError`)
- `normalize` (bool): If True, scale the counts to sum to 1.0; otherwise return raw counts

**Returns**: `list[float]` of length 256

**Usage**: TNF is a compositional signature of genomic DNA; different taxa have distinct 4-mer frequencies due to codon bias, GC% pattern.

---

#### `bin_contigs(contigs, coverage=None, method="composition", n_bins=None, min_contig_length=1000, seed=42, marker_genes=None, marker_set=None) -> BinningResult`

**Purpose**: Group contigs into putative genomes (MAGs) using unsupervised clustering.

**Methods**:

1. **`composition`**: Cluster contigs on their 256-D TNF vectors (default).
2. **`coverage`**: Cluster on log-transformed coverage (one value per contig, or a vector across samples).
3. **`combined`**: Concatenate TNF and log-transformed coverage. Without `coverage` it reduces to composition-only.

**Parameters**:
- `contigs` (dict[str, str]): `{contig_id: sequence}`
- `coverage` (dict[str, float] | dict[str, list[float]] | None): `{contig_id: depth}` or `{contig_id: [depth in sample 1, depth in sample 2, ...]}`; values enter the features as `log(max(depth, 0.01))`
- `method` (str): `"composition"`, `"coverage"`, or `"combined"`
- `n_bins` (int | None): Number of bins (clusters). If None, the smaller of one bin per 3 Mbp of binned sequence and half the number of binnable contigs, raised to at least 2 and capped at the number of binnable contigs
- `min_contig_length` (int): Contigs shorter than this stay unbinned (default 1000)
- `seed` (int): Seed of the k-means++ initialization (default 42)
- `marker_genes` (Mapping[str, Sequence[str]] | None): Single-copy marker gene annotations per contig (contig ID → marker names found on it). When given, every bin is assessed (see `assess_bin_quality()`) and `high_quality_bins` is filled
- `marker_set` (Collection[str] | None): Expected markers for `marker_genes` (default `BACTERIAL_SINGLE_COPY_MARKERS`)

**Returns**: `BinningResult`. Bin `completeness`, `contamination`, `quality_score` and `BinningResult.high_quality_bins` are `None` unless `marker_genes` is given.

**Raises**: `ValueError` for empty `contigs` or an unknown method.

**Algorithm details**: pure-Python k-means with k-means++ initialization; each contig goes to its nearest centroid.

**Example**:
```python
from metainformant.metagenomics.shotgun.binning import bin_contigs

contigs = {"c1": "ATCG" * 500, "c2": "GCTA" * 500, "c3": "AAAA" * 500}
coverage = {"c1": 10.0, "c2": 10.5, "c3": 50.0}
result = bin_contigs(contigs, coverage, method="combined", n_bins=2)
print(f"Assigned {result.binned_contigs} contigs to {len(result.bins)} bins")
print(result.high_quality_bins)  # None: no marker_genes were given
```

---

#### `assess_bin_quality(bins, contigs=None, marker_genes=None, marker_set=None) -> list[GenomeBin]`

**Purpose**: Compute CheckM-style bin quality from single-copy marker gene annotations. Updates the bins in place and returns them.

For `M` expected markers (the size of `marker_set`, by default the 36 `BACTERIAL_SINGLE_COPY_MARKERS`):
- **Completeness** = (expected markers present at least once in the bin) / `M`
- **Contamination** = (extra copies: the sum of `copies - 1` over the markers present) / `M`
- **Quality score** = completeness − 5 × contamination (Parks et al.), negative for a heavily contaminated bin

**Parameters**:
- `bins` (list[GenomeBin])
- `contigs`: accepted but unused (no marker search is performed on sequences)
- `marker_genes` (Mapping[str, Sequence[str]] | None): contig ID → marker names found on that contig; a name listed twice counts as two copies, and contigs missing from the mapping carry no markers. The annotations must come from a real marker search (for example HMMER with TIGRFAM/Pfam marker profiles); annotation names outside the expected set are ignored
- `marker_set` (Collection[str] | None): Expected markers, matching the annotation names (for example a lineage-specific CheckM or BUSCO set). An empty set raises `ValueError`

**Without `marker_genes`**, completeness, contamination and quality score are set to `None` for every bin (and `marker_counts` is emptied): they are never estimated from genome size or sequence composition.

**Example**:
```python
from metainformant.metagenomics.shotgun.binning import (
    BACTERIAL_SINGLE_COPY_MARKERS, GenomeBin, assess_bin_quality,
)

bin1 = GenomeBin("b1", ["c1"], total_length=3_000_000)
print(assess_bin_quality([bin1])[0].completeness)  # None
markers = {"c1": list(BACTERIAL_SINGLE_COPY_MARKERS)}
print(assess_bin_quality([bin1], marker_genes=markers)[0].completeness)  # 1.0
```

---

#### `count_high_quality_bins(bins, min_completeness=0.9, max_contamination=0.05) -> int | None`

**Purpose**: Count the assessed bins with completeness ≥ `min_completeness` and contamination ≤ `max_contamination`. Returns `None` when any bin has no assessed quality, and 0 for an empty list.

---

#### `refine_bins(bins, contigs, completeness_threshold=0.5, contamination_threshold=0.1, marker_genes=None, marker_set=None) -> list[GenomeBin]`

**Purpose**: Filter and split bins on their marker gene completeness and contamination.

**Parameters**:
- `bins` (list[GenomeBin]): Bins from `bin_contigs()`
- `contigs` (dict[str, str]): Contig sequences (used to split contaminated bins by TNF)
- `completeness_threshold` (float): Minimum completeness to keep a bin (default 0.5)
- `contamination_threshold` (float): Maximum contamination accepted without splitting (default 0.1)
- `marker_genes`, `marker_set`: As in `assess_bin_quality()`. `marker_genes` is required whenever `bins` is not empty (`ValueError` otherwise; an empty `bins` returns `[]`)

**Returns**: `list[GenomeBin]` with updated quality, renamed `refined_bin_0000`, `refined_bin_0001`, ...

**Algorithm**:
1. Assess every bin with `marker_genes`
2. Drop bins below `completeness_threshold`
3. A bin whose contamination exceeds `contamination_threshold` is dropped if it has one contig, kept unchanged if it has two, and split in two by TNF k-means if it has more; the sub-bins are assessed with the same annotations and kept when their completeness reaches `completeness_threshold`

Merging of small bins is not implemented.

---

### Profiling Functions

#### `build_kmer_index(reference_sequences, taxonomy=None, k=31) -> KmerIndex`

**Purpose**: Build a k-mer lookup index from reference genomes for taxonomic classification.

**Data model**:
```python
@dataclass
class KmerIndex:
    kmer_to_taxon: dict[str, str]          # k-mer → taxon ID (LCA when shared)
    taxon_lineages: dict[str, list]        # taxon_id → [(rank, name), ...]
    k: int
    total_kmers: int
    num_taxa: int
```

**Algorithm**:
1. Extract every forward-strand k-mer of each reference genome
2. A k-mer found in one reference maps to that reference; a k-mer shared by several references maps to the LCA of their lineages (an empty lineage, indexed as `"unclassified"`, when they share no rank)
3. Store the lineage of every taxon

**Parameters**:
- `reference_sequences` (dict[str, str]): `{genome_id: sequence}`
- `taxonomy` (dict[str, list] | None): `{genome_id: [(rank, name), ...]}`; if None, each reference is its own taxon
- `k` (int): k-mer size, default 31 (at least 11, otherwise `ValueError`)

**Returns**: `KmerIndex`

**Example**:
```python
from metainformant.metagenomics.shotgun.profiling import build_kmer_index

refs = {"genusA_sp1": "ACGTTGCAACGGTTAACCGGTTAACGTACGTTGCA", "genusA_sp2": "TTGGCCAATTGGCCAAGGCCTTAAGGCCTTAAGGC"}
tax = {
    "genusA_sp1": [("genus", "GenusA"), ("species", "sp1")],
    "genusA_sp2": [("genus", "GenusA"), ("species", "sp2")],
}
index = build_kmer_index(refs, tax, k=31)
print(f"Indexed {index.total_kmers} k-mers from {index.num_taxa} taxa")
```

---

#### `profile_community(reads, database=None, reference_sequences=None, reference_taxonomy=None, k=31, min_kmer_hits=2, confidence_threshold=0.5) -> CommunityProfile`

**Purpose**: Classify reads taxonomically with a k-mer LCA approach (Kraken-style).

**Algorithm**:
For each read:
1. Extract all k-mers and look each up in the index; k-mers indexed under an empty lineage (shared by taxa with no common ancestor) are ignored
2. If fewer than `min_kmer_hits` informative hits remain, the read is unclassified
3. Take the taxon with the most hits; if it holds at least `confidence_threshold` of the hits, the read is assigned to it
4. Otherwise assign the read to the LCA of all taxa it hit; if those taxa share no rank, the read stays unclassified

**Parameters**:
- `reads` (list[str] | dict[str, str]): Reads
- `database` (KmerIndex | None): Pre-built index from `build_kmer_index()`
- `reference_sequences` (dict[str, str] | None), `reference_taxonomy` (dict[str, list] | None): Alternative: build the index on the fly from these references
- `k` (int): k-mer size when building on the fly (ignored when `database` is given)
- `min_kmer_hits` (int): Minimum informative k-mer hits to classify a read (default 2)
- `confidence_threshold` (float): Minimum fraction of hits agreeing on the top taxon (default 0.5)

**Returns**: `CommunityProfile` (`taxa`, `total_reads`, `classified_reads`, `unclassified_reads`, `classification_rate`, `rank_summaries`)

**Example**:
```python
from metainformant.metagenomics.shotgun.profiling import profile_community

# index: the KmerIndex built above
reads = ["ACGTTGCAACGGTTAACCGGTTAACGTACGTTGCA"]
profile = profile_community(reads, database=index, min_kmer_hits=1)
for taxon in profile.taxa[:10]:
    print(f"{taxon.name}: {taxon.relative_abundance:.2%} "
          f"({taxon.read_count} reads)")
```

**Output interpretation**:
- Reads with too few informative hits, or whose hit taxa share no ancestor → `unclassified_reads`
- A taxon's `relative_abundance` is its fraction of the classified reads, so the abundances sum to 1.0 over the classified reads

---

#### `calculate_relative_abundance(profile, rank=None, min_abundance=0.0) -> dict[str, float]`

**Purpose**: Normalize a profile's read counts to relative abundances summing to 1.0, optionally aggregated at a rank (taxa without that rank count as `"unclassified"`); taxa below `min_abundance` are grouped as `"Other"`.

---

## Diversity

Diversity functions are pure statistical calculations; they do not parse data. Input is abundance vectors (counts or relative abundances, as noted per metric).

### Alpha Diversity

#### `alpha_diversity(abundances, metric="shannon") -> dict`

**Purpose**: Compute within-sample diversity.

**Supported metrics**:

| Metric | Input type | Range | Interpretation |
|--------|------------|-------|---------------|
| `shannon` | counts or prop | [0, ln S] | Entropy; higher = more diverse, accounts for evenness |
| `simpson` | counts or prop | (0, 1] | Σ pᵢ² (dominance form: probability that two random individuals are the same species); 1 − D is the Gini–Simpson diversity |
| `invsimpson` | counts or prop | [1, S] | 1 / Σ pᵢ²; effective number of species |
| `chao1` | integer counts | [S_obs, ∞) | Richness estimator; S_est = S_obs + f1² / (2·f2), or S_obs + f1·(f1 − 1)/2 when f2 = 0 and f1 > 1 |
| `ace` | integer counts | [S_obs, ∞) | Abundance-based Coverage Estimator; species with abundance ≤ 10 are the rare class |
| `observed` | counts or prop | [1, S_max] | Simple species count (richness) |
| `fisher_alpha` | integer counts | (0, ∞] | Parameter of Fisher's log-series, the root of S = α·ln(1 + N/α), found by bracketed bisection; `inf` when every species is a singleton |
| `pielou_evenness` | counts or prop | [0, 1] | J = H / ln(S); 1 = perfectly even (also 1 when S ≤ 1) |

**Parameters**:
- `abundances` (list | NumPy array): Per-species abundances. `chao1`, `ace` and `fisher_alpha` need integer-valued counts (floats such as `3.0` are fine) and raise `ValueError` for proportions or any other non-integer value
- `metric` (str): One of the above

**Returns**: `dict` with keys:
- `value` (float): Computed diversity
- `metric` (str): metric name
- `n_species` (int): Number of species with non-zero abundance
- `total_count` (float): Sum of abundances

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import alpha_diversity

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
| `jaccard` | presence/absence | [0, 1] | Qualitative; 1 − shared species / total species |
| `aitchison` | compositional data | [0, ∞) | Euclidean distance on CLR-transformed abundances (pseudocount 0.5) |

**Parameters**:
- `samples` (list[list[float]]): Outer list = samples (at least two); inner = per-species abundances (same length across samples)
- `metric` (str): `"bray_curtis"`, `"jaccard"`, or `"aitchison"`

**Returns**: `dict` with:
- `distance_matrix` (list[list[float]]): Symmetric n×n matrix
- `metric` (str)
- `n_samples` (int)

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import beta_diversity

samples = [
    [100, 50, 0, 10],    # sample A
    [10, 5, 0, 50],      # sample B
    [0, 0, 100, 0],      # sample C
]
result = beta_diversity(samples, metric="bray_curtis")
print(result["distance_matrix"][0][1])  # A vs B dissimilarity
```

---

### Ordination

#### `ordination(distance_matrix, method="pcoa", n_components=2, seed=42) -> dict`

**Purpose**: Dimensionality reduction of distance matrix for visualization.

**Methods**:
- **PCoA** (Principal Coordinates Analysis): Eigen-decomposition of the double-centered squared distance matrix (classical MDS); NumPy `eigh` when NumPy is importable, pure-Python power iteration otherwise. Deterministic. At most `n_samples − 1` components are returned.
- **NMDS** (Non-metric Multidimensional Scaling): Gradient-descent stress minimization against the raw distances from a random start. It has no rank-order (isotonic) step, so it behaves like metric MDS.

**Parameters**:
- `distance_matrix` (list[list[float]]): Square symmetric matrix
- `method` (str): `"pcoa"` (default) or `"nmds"`
- `n_components` (int): Target dimensions (2 or 3 for plotting)
- `seed` (int | None): Seed of the NMDS starting configuration (default 42, so repeated calls give identical results; `None` draws fresh entropy). PCoA ignores it

**Returns**: `dict` with:
- `coordinates` (list[list[float]]): Sample coordinates (n × n_components)
- `eigenvalues` (list[float]): Eigenvalues (PCoA) or the last stress values (NMDS)
- `variance_explained` (list[float]): For PCoA, each returned component's share of the summed positive eigenvalues of the returned components; for NMDS, `1 − final stress` repeated for each component

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import ordination

dist_matrix = [
    [0.0, 0.2, 0.7, 0.8],
    [0.2, 0.0, 0.6, 0.7],
    [0.7, 0.6, 0.0, 0.3],
    [0.8, 0.7, 0.3, 0.0],
]
pcoa = ordination(dist_matrix, method="pcoa", n_components=2)
xs = [coord[0] for coord in pcoa["coordinates"]]
ys = [coord[1] for coord in pcoa["coordinates"]]
# Plot: plt.scatter(xs, ys)
```

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
- `seed` (int | None): Random seed. The default `None` draws fresh entropy, so results vary between runs; pass an integer for reproducibility

**Returns**: `dict`:
- `depths` (list[int]): Depths evaluated
- `mean_species` (list[float]): Mean observed species at each depth
- `std_species` (list[float]): Standard deviation
- `is_saturated` (bool): Has plateaued?

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import rarefaction_curve

abundances = [1000, 500, 200, 50, 10, 5, 1, 1]  # community with many rare taxa
curve = rarefaction_curve(abundances, n_iterations=50, seed=42)
# plt.errorbar(curve["depths"], curve["mean_species"], yerr=curve["std_species"])
```

**Interpretation**: If curve is still rising steeply at max depth, undersampled; consider deeper sequencing.

---

#### `rarefy(abundances, depth, seed=None) -> list[int]`

**Purpose**: Subsample community to even depth (rarefaction).

**Parameters**:
- `abundances` (list[int]): Counts per taxon
- `depth` (int): Target total reads after subsampling
- `seed` (int | None): Random seed (default `None`: unseeded; pass an integer for reproducibility)

**Returns**: `list[int]` rarefied counts (same length as input, sum = `depth`)

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import rarefy

counts = [1000, 500, 200]
rarefied = rarefy(counts, depth=1000, seed=42)
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
- `seed` (int | None): Random seed (default `None`: unseeded, so the p-value varies between runs; pass an integer for reproducibility)

**Returns**: `dict`:
- `pseudo_f` (float)
- `p_value` (float)
- `r_squared` (float) — fraction of variance explained by groups
- `n_permutations` (int)

**Example**:
```python
from metainformant.metagenomics.diversity.metrics import beta_diversity, permanova

samples = [[100, 50, 0, 10], [90, 55, 2, 12], [10, 5, 0, 50], [12, 8, 1, 45], [0, 0, 100, 0], [1, 0, 90, 3]]
beta = beta_diversity(samples, metric="bray_curtis")
permanova_result = permanova(
    beta["distance_matrix"],
    groups=["A", "A", "B", "B", "C", "C"],
    n_permutations=9999,
    seed=42,
)
print(f"Pseudo-F = {permanova_result['pseudo_f']:.3f}, "
      f"p = {permanova_result['p_value']:.4g}, R² = {permanova_result['r_squared']:.3f}")
```

**Interpretation**: p < 0.05 → reject H₀; groups have significantly different community composition.

---

## Functional

### Data Structures

- `ORF`: `sequence_id`, `start`, `end`, `strand` (`"+"` or `"-"`), `nucleotide_seq`, `protein_seq`, `frame`, `partial`; properties `length_nt`, `length_aa`
- `HMMHit`: `query_id`, `target_id`, `target_name`, `score`, `e_value`, `bias`, `domain_score`, `domain_e_value`, `query_start`, `query_end`, `target_start`, `target_end`, `description`. `e_value`, `bias` and `domain_e_value` are NaN (not computed), and the start/end fields are placeholders (the whole query, no target coordinates)
- `FunctionalAnnotation`: `orf_id`, `gene_families`, `kegg_orthologs`, `cog_categories`, `pfam_domains`, `ec_numbers`, `best_hit`, `description`, `confidence`
- `PathwayDefinition`: `pathway_id`, `name`, `description`, `database`, `required_kos`, `required_ecs`, `optional_kos`, `reaction_steps`, `category`, `subcategory`; property `total_required`
- `PathwayResult`: `pathway_id`, `pathway_name`, `completeness`, `matched_kos`, `matched_ecs`, `missing_kos`, `missing_ecs`, `total_required`, `total_matched`, `redundancy`, `confidence`

### Gene Annotation

#### `predict_orfs(sequence, sequence_id="seq", min_length=100, allow_partial=True, translation_table=11) -> list[ORF]`

Scan all six reading frames of one nucleotide sequence for ORFs delimited by start and stop codons (a six-frame scan, not Glimmer or Prodigal). `min_length` is in nucleotides. `translation_table` is an NCBI genetic code id (default 11, bacterial); it sets the amino-acid assignments, the stop codons (TGA is tryptophan in table 4) and, through the table's start codons, which of ATG, GTG and TTG can open an ORF. An id that is not an NCBI genetic code raises `ValueError`. Every in-frame start upstream of a stop yields its own, nested ORF; each (strand, start, end) appears once. Partial ORFs at the sequence edges are included when `allow_partial` is True.

**Returns**: `list[ORF]` sorted by position; reverse-strand ORFs report forward-strand coordinates.

---

#### `annotate_genes(sequences, hmm_db=None, e_value_threshold=None, min_score=25.0) -> list[FunctionalAnnotation]`

Score protein sequences (`{gene_id: protein_sequence}`) against a user-supplied profile database `{profile_id: {position: {amino_acid: score}}}` (position-specific scoring matrices; no database is bundled). A profile hits when its best ungapped alignment score, averaged over matched positions, reaches `min_score`. Profile IDs starting with `KO:` or `COG:` fill `kegg_orthologs` or `cog_categories` with the prefix removed, IDs starting with `PF` fill `pfam_domains` unchanged, and any other ID fills `gene_families`; `EC:<a.b.c.d>` in a gene family ID fills `ec_numbers`. Without `hmm_db` every annotation is empty.

`e_value_threshold` must be `None`: the scoring has no calibrated null model, so E-values cannot be computed, hit E-values are NaN, and a numeric threshold raises `ValueError`.

**Returns**: `list[FunctionalAnnotation]`; `confidence` is `min(1, best_score / 100)`.

---

#### `classify_gene_families(genes, database="COG", reference_profiles=None, min_score=20.0) -> dict[str, list[str]]`

Assign genes (`{gene_id: protein_sequence}`) to families by scoring them against `reference_profiles` (same format as `hmm_db`; only profiles whose ID starts with the database prefix `COG:`, `KO:` (KEGG) or `PF` (Pfam) are used). Returns, per gene, the IDs of all profiles scoring at least `min_score`. When `reference_profiles` is None no family is assigned: every gene maps to an empty list (unannotated), because amino acid composition cannot identify a COG, KEGG or Pfam family.

---

### Pathway Analysis

#### `reconstruct_pathways(annotations, database="KEGG", pathway_definitions=None, min_completeness=0.0) -> list[PathwayResult]`

Map annotated KO/EC identifiers (`{gene_id: [KO or EC ids]}`) to pathways and compute completeness. Eight KEGG pathways are built in (glycolysis, TCA cycle, pentose phosphate, oxidative phosphorylation, arginine biosynthesis, nitrogen, methane and photosynthesis); `database` only filters them, so `database="MetaCyc"` without `pathway_definitions` returns nothing. Completeness is `total_matched / total_required`: the required KOs and ECs present for a definition without `reaction_steps`, or the steps with at least one annotated alternative for one with, so it always lies in [0, 1]. `redundancy` is the mean number of gene copies per matched KO and `confidence` is `completeness * min(1, redundancy)`.

**Returns**: `list[PathwayResult]` with completeness ≥ `min_completeness`, sorted by completeness descending.

---

#### `calculate_pathway_completeness(pathway, annotations) -> float`

Completeness of one `PathwayDefinition` against a set of observed KO/EC identifiers, by steps when it has `reaction_steps` (a step is complete when any alternative is present) and by required KOs/ECs otherwise.

---

#### `compare_pathway_profiles(samples, pathway_definitions=None, database="KEGG") -> dict[str, dict[str, float]]`

Reconstruct pathways for each sample (`{sample_id: {gene_id: [ids]}}`) and return `{pathway_id: {sample_id: completeness}}`.

---

#### `find_differential_pathways(group1_samples, group2_samples, pathway_definitions=None, min_diff=0.2) -> list[dict]`

Report the pathways whose mean completeness differs between two groups of samples by at least `min_diff` (no statistical test). Each dict has `pathway_id`, `pathway_name`, `group1_mean`, `group2_mean`, `difference`, `abs_difference` and `direction` (`"group1_enriched"` or `"group2_enriched"`); the list is sorted by `abs_difference`.

**Example**:
```python
from metainformant.metagenomics.functional.pathways import (
    PathwayDefinition, calculate_pathway_completeness, reconstruct_pathways,
)

annotations = {"gene1": ["K00844", "K00845"], "gene2": ["K01810"]}
for result in reconstruct_pathways(annotations)[:3]:
    print(f"{result.pathway_name}: {result.completeness:.0%}")

custom = PathwayDefinition(pathway_id="custom1", name="Custom", required_kos=["K00844", "K00845"])
print(calculate_pathway_completeness(custom, {"K00844", "K00845"}))  # 1.0
```

---

## Comparative

All functions live in `metainformant.metagenomics.comparative.differential_abundance`. `counts` is a samples × taxa matrix (a list of lists, a NumPy array or a pandas DataFrame), `groups` has one integer label per sample, and `taxa_names` one name per column.

#### `differential_abundance(counts, groups, taxa_names, method="aldex2_like", n_monte_carlo=128, seed=None) -> list[dict]`

**Purpose**: Identify taxa (or genes, pathways) whose abundance differs between two groups. `groups` must contain exactly two distinct labels.

**Supported methods**:

| method | Statistics | Compositional-aware? |
|--------|------------|---------------------|
| `"aldex2_like"` | Dirichlet Monte Carlo CLR (`n_monte_carlo` draws per sample, 0.5 pseudocount prior) averaged per sample, then Welch's t-test and Cohen's d | Yes |
| `"ancom_like"` | Pairwise log-ratio Welch tests; effect size = fraction of rejected comparisons W, p-value = 1 − W | Yes |
| `"simple_deseq"` | Median-of-ratios size factors, then Welch's t-test and Cohen's d on the normalized counts (no negative binomial model) | No |

**Parameters**:
- `method`: one of the above
- `n_monte_carlo` (int): Monte Carlo instances for `aldex2_like` (at least 1)
- `seed` (int | None): Seed of the Dirichlet draws (`aldex2_like`). The default `None` is unseeded, so p-values vary between runs; the other two methods are deterministic

**Returns**: list of dicts sorted by adjusted p-value, each with:
- `taxon` (str)
- `log2fc` (float): log2((mean_group1 + 1) / (mean_group2 + 1))
- `p_value` (float)
- `adjusted_p` (float): Benjamini-Hochberg adjusted p-value
- `effect_size` (float)
- `mean_group1`, `mean_group2` (float): means of the raw counts

**Example**:
```python
from metainformant.metagenomics.comparative.differential_abundance import differential_abundance

counts = [
    [100, 20, 5, 0, 50],
    [90, 25, 8, 2, 45],
    [110, 18, 3, 1, 55],
    [10, 80, 60, 40, 5],
    [15, 75, 55, 35, 8],
    [12, 85, 65, 45, 3],
]
groups = [0, 0, 0, 1, 1, 1]
taxa = ["Bacteroides", "Prevotella", "Ruminococcus", "Faecalibacterium", "Akkermansia"]
diff = differential_abundance(counts, groups, taxa, method="aldex2_like", seed=42)
significant = [r["taxon"] for r in diff if r["adjusted_p"] < 0.05]
```

---

#### `clr_transform(counts, pseudocount=0.5) -> list[list[float]]`

Centered log-ratio transform per sample: `log(x + pseudocount) − mean(log(x + pseudocount))`.

#### `indicator_species(counts, groups, taxa_names, n_permutations=999, seed=None) -> list[dict]`

IndVal indicator species analysis (specificity × fidelity) with a permutation p-value; works with more than two groups. Returns dicts with `taxon`, `indicator_value`, `p_value` and `associated_group`, sorted by indicator value. `seed` defaults to `None` (unseeded).

#### `effect_size_analysis(counts, groups, taxa_names) -> list[dict]`

LEfSe-style ranking: a Kruskal-Wallis p-value (SciPy, or an equivalent pure-Python tie-corrected computation for any number of groups) and a simplified LDA-like score: `log10(1 + standardized difference of CLR group means)` for two groups, and from the largest unstandardized pairwise difference of group means for more. Returns dicts with `taxon`, `lda_score`, `p_value` and `direction` (`"enriched_in_<label>"`), sorted by score.

#### `biomarker_discovery(counts, groups, taxa_names, method="random_forest", n_estimators=100, cv_folds=5, seed=None) -> dict`

Random forest feature importances on CLR-transformed data with cross-validated accuracy; falls back to |Cohen's d| ranking (`cv_accuracy` 0.0) without scikit-learn. `cv_folds` is capped by the size of the smallest group and raised to at least 2. `seed` defaults to `None` (unseeded). Returns `selected_taxa`, `importances`, `cv_accuracy` and `model_summary`.

---

## Visualization

All plot functions need matplotlib (`ImportError` otherwise). When `output_path` is given they save the figure (150 dpi; the format follows the file extension, and missing parent directories are created), close it and return the path as a string; otherwise they return the matplotlib `Figure`. For general-purpose plots (PCA scatter, volcano, heatmaps) see `metainformant.visualization`.

### Exported Functions (in `visualization/plots.py`)

| Function | Description |
|----------|-------------|
| `plot_krona_chart(taxonomy_data, output_path=None, title="Taxonomic Composition", max_depth=6, min_fraction=0.01)` | Static Krona-style sunburst of `{"A;B;C": abundance}` lineages |
| `plot_stacked_bar(abundance_data, level="phylum", output_path=None, title=None, top_n=15, sort_samples=True)` | Stacked bars of the top taxa per sample (`{sample: {taxon: abundance}}`), the rest as "Other" |
| `plot_rarefaction_curves(otu_table, step=100, output_path=None, title="Rarefaction Curves", n_iterations=10, seed=42)` | Rarefaction curve per sample (`{sample: {otu: count}}`) |
| `plot_ordination(distance_matrix, method="PCoA", output_path=None, title=None, groups=None, seed=42)` | PCoA or NMDS scatter of a nested-dict distance matrix, coloured by `groups` |
| `plot_alpha_diversity(samples, metrics=None, output_path=None, groups=None, title="Alpha Diversity")` | Bars or per-group boxplots of observed, Shannon, Simpson and Chao1 |
| `plot_heatmap(abundance_matrix, output_path=None, cluster=True, title="Abundance Heatmap", cmap="YlOrRd", top_n=30)` | Heatmap of the top taxa; `cluster=True` orders rows and columns by greedy nearest-neighbour similarity (no dendrogram) |

**`plot_ordination`**: `method` is case-insensitive. `"PCoA"` (default) double-centres the squared distances and takes the top two axes (labelled PC1 and PC2); `"NMDS"` runs `ordination(..., method="nmds")` from a random start drawn with `seed` (default 42; axes labelled NMDS1 and NMDS2); `seed` is ignored by PCoA. Any other `method` raises `ValueError`. With fewer than 3 samples the figure only reports that more are needed.

**`plot_alpha_diversity`** computes its own metrics from `{sample: {otu: count}}`: its `simpson` is the Gini–Simpson index `1 − Σ n(n − 1) / (N(N − 1))`, which differs from the `simpson` of `alpha_diversity` (Σ pᵢ²); names other than the four listed are not computed and plot as 0.

**Example**:
```python
from metainformant.metagenomics.visualization.plots import plot_ordination

distances = {
    "s1": {"s1": 0.0, "s2": 0.2, "s3": 0.7, "s4": 0.8},
    "s2": {"s1": 0.2, "s2": 0.0, "s3": 0.6, "s4": 0.7},
    "s3": {"s1": 0.7, "s2": 0.6, "s3": 0.0, "s4": 0.3},
    "s4": {"s1": 0.8, "s2": 0.7, "s3": 0.3, "s4": 0.0},
}
groups = {"s1": "control", "s2": "control", "s3": "treated", "s4": "treated"}
path = plot_ordination(
    distances, method="NMDS", groups=groups, seed=7, output_path="output/metagenomics/nmds.png"
)
```

---

## Algorithm Comparison

### Clustering: OTU vs. ASV

| | OTU Clustering | ASV Denoising |
|---|---------------|---------------|
| **Resolution** | 97% identity clusters (OTU) | Exact sequence variants (100%) |
| **Speed** | Slow: one pure-Python alignment per sequence–centroid pair that passes the pre-filter | Fast: no alignments; cost grows with unique sequences × ASVs |
| **Sensitivity** | Lumps true variants | Resolves single-nucleotide differences |
| **Error handling** | Noise absorbed by the identity threshold | Explicit error model; Poisson abundance test, Phred-based when quality scores are given |
| **Reproducibility** | Deterministic; depends on input order, sort key and threshold | Deterministic given the same data |
| **Best for** | Legacy comparison with 97% OTU tables, up to a few thousand sequences | Strain tracking, rare biosphere |


### Taxonomic Classification: Naive Bayes vs. k-mer vote ("blast")

| | Naive Bayes (RDP-style) | k-mer vote (`blast`) |
|---|-------------------|-------------------|
| **Reference data** | Reference sequences with lineages | Reference sequences with lineages |
| **Speed** | Rebuilds the k-mer tables from the whole reference for every query | Computes an 11-mer set for every reference for every query |
| **Confidence** | Bootstrap support per rank; abstains below `confidence_threshold` | Vote share per rank |
| **Unrelated queries** | Left `unclassified` | Left `unclassified` only when no 11-mer is shared |
| **Reproducibility** | Seeded bootstrap (`seed=42`) | Deterministic |
| **Best for** | Per-rank bootstrap support with abstention | Quick consensus of the closest references; small reference sets |

---

### Differential Abundance: Compositional vs. Standard

Compositional data (relative abundances) have sum constraint → spurious correlations. `differential_abundance` offers two compositional-aware methods and one count-normalization method:

- **`aldex2_like`**: CLR transform + Dirichlet Monte Carlo sampling + Welch's t-test
- **`ancom_like`**: Pairwise log-ratio tests of each taxon against every other taxon
- **`simple_deseq`**: Median-of-ratios size-factor normalization + Welch's t-test (no negative binomial model; not compositional-aware)

---

## Function Matrix

| Submodule | Function | Class / return | Output |
|-----------|----------|----------------|--------|
| amplicon | `cluster_otus` | `ClusteringResult` | OTU table |
| amplicon | `calculate_identity` | float | pairwise identity |
| amplicon | `filter_chimeras` | dict[str, bool] | chimera flags |
| amplicon | `estimate_error_rates` | `ErrorModel` | error model |
| amplicon | `denoise_sequences` | `DenoisingResult` | ASVs |
| amplicon | `merge_paired_reads` | dict[str, str] | merged reads |
| amplicon | `classify_taxonomy` | list[TaxonomyAssignment] | classified sequences |
| amplicon | `build_taxonomy_tree` | TaxonomyNode | tree |
| amplicon | `calculate_confidence` | dict[str, dict] | per-rank confidence summary |
| shotgun | `assemble_contigs` | list[Contig] | assembly |
| shotgun | `scaffold_contigs` | list[Scaffold] | scaffolds |
| shotgun | `calculate_assembly_stats` | AssemblyStats | metrics |
| shotgun | `calculate_tetranucleotide_freq` | list[float] | 256-D TNF profile |
| shotgun | `bin_contigs` | BinningResult | bins |
| shotgun | `assess_bin_quality` | list[GenomeBin] | marker-gene completeness and contamination |
| shotgun | `count_high_quality_bins` | int \| None | high-quality bin count |
| shotgun | `refine_bins` | list[GenomeBin] | refined bins |
| shotgun | `build_kmer_index` | KmerIndex | index |
| shotgun | `profile_community` | CommunityProfile | taxon abundances |
| shotgun | `calculate_relative_abundance` | dict[str, float] | relative abundances |
| diversity | `alpha_diversity` | dict | alpha value |
| diversity | `beta_diversity` | dict | distance matrix |
| diversity | `ordination` | dict | coords |
| diversity | `rarefaction_curve` | dict | curve data |
| diversity | `rarefy` | list[int] | subsampled counts |
| diversity | `permanova` | dict | test stats |
| functional | `predict_orfs` | list[ORF] | ORFs |
| functional | `annotate_genes` | list[FunctionalAnnotation] | profile hits |
| functional | `classify_gene_families` | dict[str, list[str]] | families per gene |
| functional | `reconstruct_pathways` | list[PathwayResult] | pathway completeness |
| functional | `calculate_pathway_completeness` | float | one pathway's completeness |
| functional | `compare_pathway_profiles` | dict[str, dict[str, float]] | completeness matrix |
| functional | `find_differential_pathways` | list[dict] | differing pathways |
| comparative | `differential_abundance` | list[dict] | per-taxon statistics |
| comparative | `clr_transform` | list[list[float]] | CLR matrix |
| comparative | `indicator_species` | list[dict] | IndVal results |
| comparative | `effect_size_analysis` | list[dict] | LEfSe-style ranking |
| comparative | `biomarker_discovery` | dict | importances and CV accuracy |
| visualization | `plot_krona_chart`, `plot_stacked_bar`, `plot_rarefaction_curves`, `plot_ordination`, `plot_alpha_diversity`, `plot_heatmap` | path or `Figure` | plots |

---

## Configuration

| Variable | Purpose | Default |
|----------|---------|--------|
| `META_OTU_THRESHOLD` | Identity threshold for `cluster_otus()` when `threshold=None`; an explicit argument wins. Must be a number in (0, 1] (`ValueError` otherwise); blank counts as unset | `0.97` |

No other environment variable is read; all other parameters are function arguments.

---

## Performance Tips

- **OTU clustering**: Keep `prefilter=True` (default); it skips only alignments that cannot reach the threshold and never changes the result.
- **Taxonomic profiling**: Build the `KmerIndex` once and reuse it across samples.
- **Diversity metrics**: All O(n_features × n_samples) or better; trivial for typical sizes (<10k features, <1000 samples).
- **Assembly**: Most expensive operation; consider running external assemblers (MEGAHIT) on HPC and only loading contigs into METAINFORMANT for downstream.

---

## Limitations

- **Assembly**: Not production-grade; use MEGAHIT/MetaSPAdes for real data
- **Binning**: Basic k-means. Bin quality needs marker gene annotations from an external search (HMMER, CheckM); the package does not search for markers
- **Annotation**: Six-frame ORF scan and ungapped scoring against a profile database you supply; no bundled databases and no E-values
- **Profiling**: k-mer LCA simple; Kraken2/Bracken faster and more accurate

The module provides a solid foundation and algorithmic reference; for production, integrate external tools where indicated.
