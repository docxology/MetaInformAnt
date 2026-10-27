# Getting Started with Metagenomics Analysis

Complete guide for running amplicon and shotgun metagenomics workflows with METAINFORMANT. Covers 16S/ITS amplicon processing and whole-metagenome shotgun analysis.

## Quick Start (10-Minute Tutorial)

### 1. Installation

```bash
cd /path/to/MetaInformAnt

# Create and activate virtual environment
uv venv
source .venv/bin/activate

# Install METAINFORMANT (numpy, scipy, scikit-learn and matplotlib are core
# dependencies and come with it)
uv pip install -e .
```

### 2. Amplicon Analysis: 16S Example

**Input**: Paired-end FASTQ files (or already merged single reads) of 16S V4 region.

Directory structure:
```
amplicon_data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

**Create analysis script** `run_amplicon.py`:

```python
#!/usr/bin/env python3
"""Amplicon metagenomics analysis demo."""

from collections import Counter

from metainformant.dna.io.fastq import read_fastq
from metainformant.metagenomics import (
    amplicon,
    diversity,
    visualization,
)

# ── 1. Load reads ───────────────────────────────────────────────────────
print("Loading reads...")
# For simplicity, assume reads are already merged; otherwise merge the pairs first
# with amplicon.asv_denoising.merge_paired_reads() (see Workflow 1 below).
# Here we just read single-end for demo. read_fastq() also handles .gz and
# returns {read_id: (sequence, quality_string)}.
fastq = read_fastq("amplicon_data/sample1_R1.fastq.gz")
reads = {read_id: seq for read_id, (seq, _qual) in fastq.items()}
print(f"  Loaded {len(reads)} reads")

# ── 2. Dereplicate ───────────────────────────────────────────────────────
print("Dereplicating sequences...")
counts = Counter(reads.values())
unique_seqs = {f"seq_{i}": seq for i, seq in enumerate(counts)}
abundances = {f"seq_{i}": cnt for i, cnt in enumerate(counts.values())}
print(f"  {len(unique_seqs)} unique sequences")

# ── 3. Denoise (ASV) ─────────────────────────────────────────────────────
print("Denoising to ASVs (DADA2-style)...")
# Pass the dereplicated sequences with their counts. Without an error model or
# quality scores, denoise_sequences() uses a default Phred-30 error model.
asv_result = amplicon.asv_denoising.denoise_sequences(unique_seqs, abundance=abundances)
asv_seqs = {f"asv_{i}": asv.sequence for i, asv in enumerate(asv_result.asvs)}
asv_counts = {f"asv_{i}": asv.abundance for i, asv in enumerate(asv_result.asvs)}
print(f"  {asv_result.num_asvs} ASVs")

# ── 4. Chimera check ─────────────────────────────────────────────────────
print("Filtering chimeras...")
chimeras = amplicon.otu_clustering.filter_chimeras(asv_seqs, abundance=asv_counts)
non_chimeric = {sid: seq for sid, seq in asv_seqs.items() if not chimeras.get(sid, False)}
print(f"  Removed {len(asv_seqs) - len(non_chimeric)} chimeras")

# ── 5. Taxonomic classification ──────────────────────────────────────────
print("Classifying taxonomy...")
# Build tiny local reference DB (in practice, use SILVA or Greengenes)
ref_db = {
    "ref1": "ACGTACGTACGTACGTACGTACGTACGTACGT",
    "ref2": "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG",
}
ref_tax = {
    "ref1": [("domain", "Bacteria"), ("phylum", "Actinobacteria")],
    "ref2": [("domain", "Bacteria"), ("phylum", "Bacteroidetes")],
}
classifications = amplicon.taxonomy.classify_taxonomy(
    non_chimeric, ref_db, ref_tax, method="naive_bayes", confidence_threshold=0.7
)
# Each assignment has .lineage [(rank, name), ...] and .confidence {rank: score};
# a rank below confidence_threshold, and every rank under it, is labelled "unclassified".
# The two toy references share almost no 8-mers with real reads, so the classifier
# abstains for nearly every ASV (use a real reference such as SILVA in practice).
n_classified = sum(1 for c in classifications if dict(c.lineage)["phylum"] != "unclassified")
print(f"  Classified {n_classified} ASVs to phylum level")

# ── 6. Build ASV table ───────────────────────────────────────────────────
# In practice: build samples × ASV count matrix
# Here we only have one sample "sample1"
asv_table = {asv_id: asv_counts[asv_id] for asv_id in non_chimeric}  # asv_id → count

# ── 7. Alpha diversity ───────────────────────────────────────────────────
print("Computing alpha diversity...")
counts_array = list(asv_table.values())
alpha = diversity.metrics.alpha_diversity(counts_array, metric="shannon")
print(f"  Shannon diversity: {alpha['value']:.3f}")

# ── 8. Beta diversity (if multiple samples) ──────────────────────────────
# dummy second sample for illustration: the same ASVs with reversed abundances
# (every sample must have the same number of ASVs)
counts2 = counts_array[::-1]
beta = diversity.metrics.beta_diversity([counts_array, counts2], metric="bray_curtis")
print(f"  Bray-Curtis dissimilarity between samples: {beta['distance_matrix'][0][1]:.3f}")

# ── 9. Visualization ─────────────────────────────────────────────────────
print("Generating plots...")
# Bar plot of top 10 ASVs. plot_stacked_bar takes {sample_id: {taxon: abundance}}
# and, when output_path is given, saves a PNG and returns its path.
png = visualization.plots.plot_stacked_bar(
    {"sample1": asv_table},
    level="ASV",
    top_n=10,
    output_path="output/metagenomics/amplicon_abundance.png",
)
print(f"  Saved {png}")

print("\nDone!")
```

**Example output** (from a simulated sample of 572 reads of 250 bp with 0.2% sequencing errors, drawn from four template sequences and one chimera of two of them; your numbers depend on your data):
```
Loading reads...
  Loaded 572 reads
Dereplicating sequences...
  242 unique sequences
Denoising to ASVs (DADA2-style)...
  5 ASVs
Filtering chimeras...
  Removed 1 chimeras
Classifying taxonomy...
  Classified 0 ASVs to phylum level
Computing alpha diversity...
  Shannon diversity: 1.275
  Bray-Curtis dissimilarity between samples: 0.423
Generating plots...
  Saved output/metagenomics/amplicon_abundance.png

Done!
```

---

### 3. Shotgun Metagenomics Example

**Input**: Raw shotgun FASTQ (single-end or paired).

```python-snippet
#!/usr/bin/env python3
"""Shotgun metagenomics analysis demo (assembly + profiling)."""

import random

from metainformant.metagenomics.shotgun import assembly, binning, profiling
from metainformant.metagenomics import functional

# ── 1. Load reads (placeholder) ─────────────────────────────────────────
# In practice, load FASTQ into an {id: sequence} dict, e.g.
# {rid: seq for rid, (seq, _qual) in read_fastq(path).items()} with
# metainformant.dna.io.fastq.read_fastq (its values are (sequence, quality) tuples),
# or use an external tool (seqtk, bioawk). Here we simulate error-free 100 bp reads
# (one every 25 bp) from two random 5 kb "genomes": ~4x coverage each.
rng = random.Random(1)
ref_genomes = {
    "genome1": "".join(rng.choices("ACGT", weights=[0.3, 0.2, 0.2, 0.3], k=5000)),
    "genome2": "".join(rng.choices("ACGT", weights=[0.2, 0.3, 0.3, 0.2], k=5000)),
}
reads = {
    f"{name}_read_{i}": genome[i:i + 100]
    for name, genome in ref_genomes.items()
    for i in range(0, len(genome) - 99, 25)
}
print(f"Loaded {len(reads)} reads")

# ── 2. Assemble contigs ──────────────────────────────────────────────────
print("Assembling contigs...")
contigs = assembly.assemble_contigs(
    reads,
    k_range=[21],  # single k for demo; the default depends on read length ([21, 33, 55] from 200 bp)
    min_contig_length=200,
)
stats = assembly.calculate_assembly_stats(contigs)
print(f"  Contigs: {stats.total_contigs}, Total bp: {stats.total_length:,}")
print(f"  N50: {stats.n50:,}, GC: {stats.gc_content:.1%}")

# ── 3. Binning ──────────────────────────────────────────────────────────
print("Binning contigs into MAGs...")
# Convert contigs to dict format
contig_dict = {c.contig_id: c.sequence for c in contigs}
# k-mer coverage estimated by the assembler (for multi-sample binning, pass
# per-sample coverage lists computed from read mapping instead)
coverage = {c.contig_id: c.coverage for c in contigs}
binning_result = binning.bin_contigs(
    contig_dict,
    coverage=coverage,
    method="combined",
    n_bins=2,  # or None to estimate from the assembly size
)
print(f"  Binned {binning_result.binned_contigs} / {binning_result.total_contigs} contigs")
# Bin quality needs single-copy marker gene annotations (contig ID -> marker names found on
# that contig), normally parsed from a marker search such as HMMER. Without them completeness,
# contamination and high_quality_bins are None.
print(f"  High-quality bins without annotations: {binning_result.high_quality_bins}")
# Hand-written stand-in for a marker search: the first contig carries every expected marker
# once, the second contig only half of them. The values below describe this stand-in, not the
# simulated sequences.
markers = list(binning.BACTERIAL_SINGLE_COPY_MARKERS)
marker_genes = {contigs[0].contig_id: markers, contigs[1].contig_id: markers[:18]}
assessed = binning.assess_bin_quality(binning_result.bins, marker_genes=marker_genes)
for b in assessed:
    print(f"  {b.bin_id}: completeness={b.completeness:.0%}, contamination={b.contamination:.0%}")
high_quality = binning.count_high_quality_bins(assessed)
print(f"  High-quality bins (completeness>=0.9, contamination<=0.05): {high_quality}")

# ── 4. Functional annotation ─────────────────────────────────────────────
print("Predicting ORFs in contigs...")
orfs = [
    orf
    for c in contigs
    for orf in functional.annotation.predict_orfs(c.sequence, sequence_id=c.contig_id, min_length=150)
]
print(f"  {len(orfs)} ORFs of at least 150 nt (ORFs that share a stop codon count separately)")

# ── 5. Taxonomic profiling (k-mer LCA) ──────────────────────────────────
print("Building k-mer index from reference genomes...")
# Tiny reference set: here, the genomes the reads were simulated from
taxonomy = {
    "genome1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
    "genome2": [("domain", "Bacteria"), ("phylum", "Bacteroidetes")],
}
kmer_index = profiling.build_kmer_index(ref_genomes, taxonomy, k=21)
print(f"  Indexed {len(kmer_index.kmer_to_taxon):,} unique k-mers")

# Profile a subset of reads
sample_reads = list(reads.values())[::2]
profile = profiling.profile_community(
    sample_reads,
    database=kmer_index,
    min_kmer_hits=2,
)
print(f"  Classified {profile.classified_reads}/{profile.total_reads} reads "
      f"({profile.classification_rate:.1%})")

# Print top taxa
for tax in profile.taxa[:5]:
    print(f"    {tax.name} ({tax.rank}): {tax.read_count} reads, "
          f"{tax.relative_abundance:.2%}")
```

**Expected output** (the reads are simulated with a fixed seed, so the numbers are reproducible; with real data they will differ):
```
Loaded 394 reads
Assembling contigs...
  Contigs: 2, Total bp: 9,900
  N50: 4,950, GC: 49.6%
Binning contigs into MAGs...
  Binned 2 / 2 contigs
  High-quality bins without annotations: None
  bin_0000: completeness=100%, contamination=0%
  bin_0001: completeness=50%, contamination=0%
  High-quality bins (completeness>=0.9, contamination<=0.05): 1
Predicting ORFs in contigs...
  86 ORFs of at least 150 nt (ORFs that share a stop codon count separately)
Building k-mer index from reference genomes...
  Indexed 9,960 unique k-mers
  Classified 197/197 reads (100.0%)
    Firmicutes (phylum): 99 reads, 50.25%
    Bacteroidetes (phylum): 98 reads, 49.75%
```

The quality values come entirely from the `marker_genes` annotation you supply, here a hand-written stand-in; with real data, build it from the output of a marker search (see [shotgun.md](shotgun.md)). The package never estimates completeness or contamination from sequence composition or bin size.

---

## System Requirements

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core |
| NumPy | ≥1.24 | PCoA eigendecomposition, Dirichlet Monte Carlo in `aldex2_like` (pure-Python fallbacks exist) |
| SciPy | ≥1.11 | Welch's t-test, Kruskal-Wallis, Poisson tail in `denoise_sequences` (pure-Python fallbacks exist; only Welch's t-test falls back to an approximation, a normal one) |
| scikit-learn | ≥1.3 | Random forest in `biomarker_discovery` (not used for binning, which has its own k-means) |
| matplotlib | ≥3.8 | `metagenomics.visualization.plots` |

All of these are installed with METAINFORMANT.

**External tools** (recommended but not required):
- **fastp** or **Trimmomatic**: adapter/quality trimming (preprocessing)
- **MEGAHIT** / **metaSPAdes**: assembly (production-scale)
- **Kraken2** / **Centrifuge**: taxonomic profiling (production)
- **CheckM** / **HMMER**: single-copy marker gene search for bin quality (`assess_bin_quality` and `refine_bins` take the markers found on each contig as input; they do not search for them)

METAINFORMANT does not wrap these tools: run them yourself (for example with Python's `subprocess`) and load their outputs into the plain dicts and lists the functions take.

---

## Detailed Analysis Workflows

### Workflow 1: Full Amplicon 16S from Demultiplexed FASTQ

This assumes you have per-sample FASTQ files from an Illumina MiSeq run (already demultiplexed).

**Step 1: Merging paired reads** (if applicable)

```python
from metainformant.dna.io.fastq import read_fastq
from metainformant.metagenomics.amplicon import asv_denoising


def load_reads(path):
    """FASTQ (optionally gzipped) -> {read_id: sequence}; quality strings are dropped."""
    return {read_id: seq for read_id, (seq, _qual) in read_fastq(path).items()}


forward = load_reads("sample_R1.fastq.gz")
reverse = load_reads("sample_R2.fastq.gz")
merged = asv_denoising.merge_paired_reads(
    forward, reverse,
    min_overlap=20,
    max_mismatch_ratio=0.2,
)
# merged: dict[id → merged_sequence]; only pairs that overlap are kept
```

**Step 2: Denoising all samples**

```python
all_reads = {}
for sample_id in sample_ids:
    fwd = load_reads(f"{sample_id}_R1.fastq.gz")
    rev = load_reads(f"{sample_id}_R2.fastq.gz")
    merged = asv_denoising.merge_paired_reads(fwd, rev)
    all_reads[sample_id] = merged

# Pool across samples for denoising (abundances are counted over the whole dataset)
pooled_reads = {}
for sid, seqs in all_reads.items():
    for read_id, seq in seqs.items():
        pooled_reads[f"{sid}_{read_id}"] = seq

asv_result = asv_denoising.denoise_sequences(
    pooled_reads,
    # quality_scores=...,  # optional: {pooled_read_id: [Phred scores]}; when given, the scores set the error probability of each mismatch, an error model is estimated from the reads, and each ASV gets a mean .quality
)
# asv_result.asvs: list of ASV objects (.sequence, .abundance, .quality) pooled over all samples
# asv_result.num_asvs, asv_result.reads_removed: number of ASVs, reads merged into other ASVs
```

**Step 3: Build sample × ASV table**

```python
import numpy as np

asv_sequences = [asv.sequence for asv in asv_result.asvs]  # list of strings
n_asvs = len(asv_sequences)
n_samples = len(sample_ids)

# Initialize count matrix
count_matrix = np.zeros((n_samples, n_asvs), dtype=int)

# Map ASV sequence → index
asv_index = {seq: i for i, seq in enumerate(asv_sequences)}

for sample_id, reads in all_reads.items():
    s_idx = sample_ids.index(sample_id)
    for seq in reads.values():
        if seq in asv_index:
            count_matrix[s_idx, asv_index[seq]] += 1
        # else: the read's sequence was merged into another ASV (or dropped), so it is not counted here

print(f"ASV table: {count_matrix.shape[0]} samples × {count_matrix.shape[1]} ASVs")
```

**Step 4: Taxonomic classification**

```python
from metainformant.metagenomics.amplicon import taxonomy

# Load reference database (SILVA, Greengenes) with your own loader (load_silva is not
# part of METAINFORMANT): ref_db is {reference_id: sequence} and ref_tax is
# {reference_id: [(rank, name), ...]}
ref_db, ref_tax = load_silva("silva_138.fasta", "silva_taxonomy.tsv")

# classify_taxonomy takes {sequence_id: sequence}, not a list
asv_ids = [f"ASV_{i}" for i in range(n_asvs)]
classifications = taxonomy.classify_taxonomy(
    sequences=dict(zip(asv_ids, asv_sequences)),
    reference_db=ref_db,
    reference_taxonomy=ref_tax,
    method="naive_bayes",
    confidence_threshold=0.7,
)
# One TaxonomyAssignment per ASV, in the order of asv_ids (.lineage, .confidence)
```

**Step 5: Alpha & beta diversity**

```python
from metainformant.metagenomics.diversity import metrics

# Alpha per sample (Shannon, and Simpson dominance for the output example below)
alpha_vals = []
simpson_vals = []
for sample_counts in count_matrix:
    result = metrics.alpha_diversity(
        abundances=sample_counts.tolist(),
        metric="shannon",
    )
    alpha_vals.append(result["value"])
    simpson = metrics.alpha_diversity(abundances=sample_counts.tolist(), metric="simpson")
    simpson_vals.append(simpson["value"])

# Beta (Bray-Curtis) distance matrix
beta = metrics.beta_diversity(
    samples=count_matrix.tolist(),
    metric="bray_curtis",
)
dist_mat = beta["distance_matrix"]

# Ordination (PCoA)
pcoa = metrics.ordination(dist_mat, method="pcoa", n_components=2)
coords = pcoa["coordinates"]  # list of [x, y] per sample
```

**Step 6: Differential abundance** (e.g., disease vs. control)

```python
from metainformant.metagenomics.comparative.differential_abundance import differential_abundance

# One group label per sample, in the order of sample_ids (exactly two distinct labels)
n_control, n_disease = 3, 3  # control samples first, then disease samples
groups = [0] * n_control + [1] * n_disease  # 0: control, 1: disease
diff = differential_abundance(
    counts=count_matrix.tolist(),  # samples × features
    groups=groups,
    taxa_names=asv_ids,
    method="ancom_like",  # pairwise log-ratios, compositional
)
# diff: list of dicts sorted by adjusted_p, each with
# "taxon", "log2fc", "p_value", "adjusted_p", "effect_size", "mean_group1", "mean_group2"
```

**Step 7: Visualization**

```python
import pandas as pd

from metainformant.metagenomics.visualization import plots
from metainformant.visualization.plots import volcano_plot

# Stacked bar (taxonomic composition at phylum level)
# Aggregate ASV counts to phylum using the classifications
phylum_abund = {sid: {} for sid in sample_ids}  # {sample_id: {phylum: count}}
for asv_idx, assignment in enumerate(classifications):
    phylum = dict(assignment.lineage)["phylum"]
    for s_idx, sid in enumerate(sample_ids):
        phylum_abund[sid][phylum] = phylum_abund[sid].get(phylum, 0) + int(count_matrix[s_idx, asv_idx])
plots.plot_stacked_bar(
    phylum_abund, level="phylum", output_path="output/metagenomics/phylum_composition.png"
)

# Ordination plot: takes a nested-dict distance matrix and computes its own 2-D ordination
# (PCoA by default; pass method="NMDS" for NMDS, whose random start uses seed=42 by default)
dist_by_sample = {
    sid_a: {sid_b: dist_mat[i][j] for j, sid_b in enumerate(sample_ids)}
    for i, sid_a in enumerate(sample_ids)
}
group_by_sample = {sid: ("Control" if i < n_control else "Disease") for i, sid in enumerate(sample_ids)}
plots.plot_ordination(
    dist_by_sample, groups=group_by_sample, output_path="output/metagenomics/pcoa.png"
)

# Volcano plot for the differential abundance results (shared plotting helper;
# metagenomics.visualization has no volcano plot)
volcano_plot(
    pd.DataFrame(diff),
    log2fc_col="log2fc",
    pval_col="adjusted_p",
    output_path="output/metagenomics/volcano.png",
)
```

---

### Workflow 2: Shotgun Metagenomics — Metagenome Assembly

This workflow demonstrates assembly and binning on a small synthetic dataset.

**Note**: Full shotgun assembly requires hundreds of millions of reads and GB of RAM. The demo below uses a tiny dataset; scale k-mer sizes (`k_range`) and memory for production.

```python
import random

from metainformant.metagenomics.shotgun import assembly, binning

# Load reads (external tool required for real dataset)
# reads = load_fastq("sample.fastq.gz")  # returns dict[id → sequence]
# (metainformant.dna.io.fastq.read_fastq returns {id: (sequence, quality)})

# For demo, simulate error-free reads from two random 20 kb genomes. The assembler
# needs non-repetitive sequence: reads from a periodic genome such as "ACGT" * 10000
# collapse into a tiny cyclic graph and produce no contigs.
def generate_mock_reads(n_reads=2500, read_len=150, genome_len=20000, gc=0.5, seed=1):
    rng = random.Random(seed)
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]  # A, C, G, T
    genome = "".join(rng.choices("ACGT", weights=weights, k=genome_len))
    reads = {}
    for i in range(n_reads):
        start = rng.randint(0, genome_len - read_len)
        reads[f"g{seed}_read_{i}"] = genome[start:start + read_len]
    return reads

reads = {**generate_mock_reads(gc=0.35, seed=1), **generate_mock_reads(gc=0.65, seed=2)}

# Assemble
contigs = assembly.assemble_contigs(
    reads,
    k_range=[21, 31],   # multiple k for better continuity
    min_contig_length=500,
)
stats = assembly.calculate_assembly_stats(contigs)
print(f"Assembly stats:\n"
      f"  N50 = {stats.n50:,} bp\n"
      f"  # contigs = {stats.total_contigs:,}\n"
      f"  Total length = {stats.total_length:,} bp\n"
      f"  GC = {stats.gc_content:.1%}")

# Bin
contig_dict = {c.contig_id: c.sequence for c in contigs}
coverage = {c.contig_id: c.coverage for c in contigs}

bins = binning.bin_contigs(
    contig_dict,
    coverage=coverage,
    method="combined",  # TNF + coverage
    n_bins=None,  # auto-detect
    min_contig_length=1000,
)

# Refine: needs marker gene annotations (contig ID -> marker names found on that contig) from a
# marker search such as HMMER; without them refine_bins raises ValueError. This hand-written
# stand-in gives the first contig every expected marker once and the second contig 30 of the 36,
# so the values below describe the stand-in, not the simulated sequences.
from metainformant.metagenomics.shotgun.binning import BACTERIAL_SINGLE_COPY_MARKERS, refine_bins
markers = list(BACTERIAL_SINGLE_COPY_MARKERS)
marker_genes = {contigs[0].contig_id: markers, contigs[1].contig_id: markers[:30]}
# Assesses the bins and applies the completeness and contamination thresholds
# (see shotgun.md for how contaminated bins are split or dropped)
refined = refine_bins(bins.bins, contig_dict, marker_genes=marker_genes)

print(f"Binning: {len(bins.bins)} bins, {bins.binned_contigs} binned contigs")
for i, genome_bin in enumerate(refined[:3]):
    print(f"  Bin {i}: {genome_bin.total_length:,} bp, "
          f"completeness={genome_bin.completeness:.2%}, "
          f"contamination={genome_bin.contamination:.2%}")
```

**Expected output** (the genomes are simulated with fixed seeds, so the numbers are reproducible; with real data they will differ):
```
Assembly stats:
  N50 = 19,996 bp
  # contigs = 2
  Total length = 39,963 bp
  GC = 50.3%
Binning: 2 bins, 2 binned contigs
  Bin 0: 19,996 bp, completeness=100.00%, contamination=0.00%
  Bin 1: 19,967 bp, completeness=83.33%, contamination=0.00%
```

Completeness and contamination come from the `marker_genes` stand-in, not from the simulated sequences. Both bins hold one contig, so a bin whose contamination exceeded `contamination_threshold` would be dropped rather than split.

---

## Data Formats

### Input Formats

The metagenomics functions take in-memory Python dicts and lists (`{id: sequence}`, count matrices), not files. These are the formats you will typically load from:

| Format | Module | How to load |
|--------|--------|-------------|
| FASTQ (single/paired) | amplicon, shotgun | `metainformant.dna.io.fastq.read_fastq` (gzipped or plain): `{read_id: (sequence, quality_string)}` |
| FASTA (merged reads, references) | amplicon, shotgun, functional | `metainformant.dna.sequence.core.read_fasta`: `{id: sequence}` |
| BAM/SAM (mapped) | shotgun | Not read by this module: compute per-contig coverage externally and pass it as `coverage={contig_id: value or [per-sample values]}` |

**FASTQ format** (expected; the quality string must be as long as the sequence):
```
@SRR123456.1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**Note**: `read_fastq` keeps the first word of each header as the read ID. It raises `ValueError` for a header without `@`, a sequence/quality length mismatch, or a truncated record; it stops silently at the first blank line, and it detects gzip only by the `.gz` file suffix. You can also use `bioawk`, `seqtk`, or your own reader.

---

### Output Formats

Diversity and comparative results are plain Python dicts or lists (of dicts, or of numbers for `rarefy` and `clr_transform`), as are the helper functions `filter_chimeras`, `merge_paired_reads`, `calculate_confidence`, `calculate_relative_abundance`, `classify_gene_families`, `compare_pathway_profiles` and `find_differential_pathways`; the main amplicon, shotgun and functional results (`ClusteringResult`, `DenoisingResult`, `TaxonomyAssignment`, `Contig`, `GenomeBin`, `CommunityProfile`, `PathwayResult`, ...) are dataclasses. Convert to:

- **TSV/CSV**: `pandas.DataFrame.to_csv()` for statistical software (R, Prism); a list of dicts such as `diff` converts directly with `pd.DataFrame(diff)`
- **JSON**: `json.dump(result, fh)` for the dict results; `json.dump(dataclasses.asdict(result), fh)` for dataclass results (use `TaxonomyNode.to_dict()` for taxonomy trees). `json.dump(result.__dict__)` fails for results that hold nested dataclasses.
- **BIOM**: HDF5-based format for QIIME compatibility (future export)

Example:
```python
import pandas as pd
alpha_df = pd.DataFrame([
    {"sample": s, "shannon": alpha_vals[i], "simpson": simpson_vals[i]}
    for i, s in enumerate(sample_ids)
])
alpha_df.to_csv("alpha_diversity.tsv", sep="\t", index=False)
```

---

## Common Pitfalls and Tips

### Pitfall 1: Passing proportions to count-based alpha diversity estimators

`chao1`, `ace` and `fisher_alpha` are defined on integer counts (singletons, doubletons, a log-series), so `alpha_diversity` raises `ValueError` for any abundance that is not an integer, such as the relative abundances `[0.7, 0.2, 0.1]`. They need **raw counts** (or rarefied counts); floats with an integer value, such as `3.0`, are accepted. `shannon`, `simpson`, `invsimpson`, `pielou_evenness` and `observed` normalize internally and work with counts or proportions. `alpha_diversity` takes a Python list or a NumPy array.

**Correct**:
```python
from metainformant.metagenomics.diversity.metrics import alpha_diversity

counts = [100, 50, 25]  # raw read counts
alpha = alpha_diversity(counts, metric="fisher_alpha")
print(alpha["value"])
```

**Wrong**:
```python
props = [0.5, 0.25, 0.125]  # proportions
alpha = alpha_diversity(props, metric="fisher_alpha")  # raises ValueError: needs integer counts
```

**Fine**: `alpha_diversity(props, metric="shannon")` returns the same value as counts in the same proportions (here 4:2:1).

---

### Pitfall 2: Forgetting to rarefy for beta diversity

Bray-Curtis and other abundance-based beta diversity metrics are **sensitive to sequencing depth**. Samples with different total reads are not comparable.

**Fix**: Rarefy all samples to even depth before computing distance:

```python
from metainformant.metagenomics.diversity import metrics

# rarefy() raises ValueError if depth exceeds a sample's total count, so let the
# shallowest sample set the depth
depth = int(count_matrix.sum(axis=1).min())
rarefied_counts = []
for sample_counts in count_matrix:
    rarefied = metrics.rarefy(sample_counts.tolist(), depth=depth, seed=42)
    rarefied_counts.append(rarefied)

beta = metrics.beta_diversity(rarefied_counts, metric="bray_curtis")
```

**Alternative**: Jaccard uses presence/absence only, so it ignores abundance differences, but it still depends on depth because rare taxa are only detected in deeper samples; rarefy first in either case.

---

### Pitfall 3: OTU clustering vs. ASV confusion

OTU clusters at 97% identity lump distinct sequences; ASVs are exact sequences. **Do not mix**: if you use ASVs, do not apply OTU clustering post-denoising.

---

## Next Steps

- Read [ARCHITECTURE.md](ARCHITECTURE.md) for system design details
- See [TUTORIALS.md](../TUTORIALS.md) for real-world study workflows
- Review [CAPABILITIES.md](CAPABILITIES.md) for complete function reference
- See the Configuration section of [ARCHITECTURE.md](ARCHITECTURE.md#configuration) for module settings (`META_OTU_THRESHOLD` and function parameters)
- Consult [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) for common issues

---

## External Resources

- **QIIME 2**: https://qiime2.org — Industry standard; consider for GUI/preprocessing
- **mothur**: https://mothur.org — Classic 16S pipeline (similar to OTU clustering here)
- **DADA2**: https://benjjneb.github.io/dada2/ — ASV denoising origin (R)
- **SILVA**: https://www.arb-silva.de — 16S/18S reference database
- **Greengenes**: http://greengenes.secondgenome.com — 16S database (legacy)
- **GTDB**: https://gtdb.ecogenomic.org — Genome taxonomy (for MAGs)

---

## Summary

Metagenomics analysis with METAINFORMANT can proceed via two main paths:

- **Amplicon**: Load reads → merge → denoise (ASV) or cluster (OTU) → classify → diversity
- **Shotgun**: Load reads → assemble → bin → annotate → profile → downstream stats

Start with the **amplicon workflow** if you have 16S/ITS data; use **shotgun** for whole metagenomes. Both can produce feature tables (samples × features count matrices that you assemble from the results) that feed into the `diversity`, `comparative`, and `visualization` modules.
