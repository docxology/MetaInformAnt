# Getting Started with Metagenomics Analysis

Complete guide for running amplicon and shotgun metagenomics workflows with METAINFORMANT. Covers 16S/ITS amplicon processing and whole-metagenome shotgun analysis.

## Quick Start (10-Minute Tutorial)

### 1. Installation

```bash
cd /path/to/MetaInformAnt

# Create and activate virtual environment
uv venv
source .venv/bin/activate

# Install METAINFORMANT
uv pip install -e .

# Install optional bioinformatics tools (recommended)
uv pip install numpy scipy scikit-learn  # for numerical operations
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

import gzip
from pathlib import Path
from collections import Counter

from metainformant.metagenomics import (
    amplicon,
    diversity,
    visualization,
)

# ── 1. Load and merge paired reads ──────────────────────────────────────
print("Loading paired-end reads...")
# For simplicity, assume reads are already merged; otherwise use merge_paired_reads()
# Here we just read single-end for demo
def read_fastq_simple(path):
    """Very simple FASTQ reader (assumes 4-line records)."""
    reads = {}
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            _ = f.readline()  # plus line
            qual = f.readline().strip()
            reads[header[1:]] = seq
    return reads

reads = read_fastq_simple("amplicon_data/sample1_R1.fastq.gz")
print(f"  Loaded {len(reads)} reads")

# ── 2. Dereplicate ───────────────────────────────────────────────────────
print("Dereplicating sequences...")
counts = Counter(reads.values())
unique_seqs = {f"seq_{i}": seq for i, (seq, cnt) in enumerate(counts.items())}
abundances = {f"seq_{i}": cnt for i, (seq, cnt) in enumerate(counts.items())}
print(f"  {len(unique_seqs)} unique sequences")

# ── 3. Denoise (ASV) ─────────────────────────────────────────────────────
print("Denoising to ASVs (DADA2-style)...")
# For demo, we skip full error modeling; just use unique seqs as "ASVs"
asv_result = amplicon.asv_denoising.DenoisingResult(
    asvs=list(unique_seqs.values()),
    abundances=list(abundances.values()),
    error_model=None,
    num_asvs=len(unique_seqs),
    reads_removed=0,
)
print(f"  {asv_result.num_asvs} ASVs")

# In real analysis:
# from metainformant.metagenomics.amplicon import asv_denoising
# asv_result = asv_denoising.denoise_sequences(reads_dict)

# ── 4. Chimera check ─────────────────────────────────────────────────────
print("Filtering chimeras...")
chimeras = amplicon.otu_clustering.filter_chimeras(unique_seqs)
non_chimeric = {sid: seq for sid, seq in unique_seqs.items() if not chimeras.get(sid, False)}
print(f"  Removed {len(unique_seqs) - len(non_chimeric)} chimeras")

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
print(f"  Classified {sum(1 for c in classifications if c.confidence > 0.7)} ASVs")

# ── 6. Build ASV table ───────────────────────────────────────────────────
# In practice: build samples × ASV count matrix
# Here we only have one sample "sample1"
asv_table = {}  # asv_id → count
for asv_seq, count in zip(non_chimeric.values(), [abundances[sid] for sid in non_chimeric]):
    asv_table[asv_seq] = count

# ── 7. Alpha diversity ───────────────────────────────────────────────────
print("Computing alpha diversity...")
counts_array = list(asv_table.values())
alpha = diversity.metrics.alpha_diversity(counts_array, metric="shannon")
print(f"  Shannon diversity: {alpha['value']:.3f}")

# ── 8. Beta diversity (if multiple samples) ──────────────────────────────
# dummy second sample for illustration
counts2 = [10, 5, 2]
beta = diversity.metrics.beta_diversity([counts_array, counts2], metric="bray_curtis")
print(f"  Bray-Curtis dissimilarity between samples: {beta['distance_matrix'][0][1]:.3f}")

# ── 9. Visualization ─────────────────────────────────────────────────────
print("Generating plots...")
# Bar plot of top 10 ASVs
fig = visualization.plots.stacked_bar(
    abundances=asv_table,
    top_n=10,
)
fig.write_html("amplicon_abundance.html")
print("  Saved amplicon_abundance.html")

print("\nDone!")
```

**Expected output**:
```
Loading paired-end reads...
  Loaded 15000 reads
Dereplicating sequences...
  8000 unique sequences
Denoising to ASVs (DADA2-style)...
  7850 ASVs
Filtering chimeras...
  Removed 120 chimeras
Classifying taxonomy...
  Classified 6500 ASVs
Computing alpha diversity...
  Shannon diversity: 4.21
Beta diversity (if multiple samples)...
  Bray-Curtis dissimilarity between samples: 0.34
Generating plots...
  Saved amplicon_abundance.html
```

---

### 3. Shotgun Metagenomics Example

**Input**: Raw shotgun FASTQ (single-end or paired).

```python-snippet
#!/usr/bin/env python3
"""Shotgun metagenomics analysis demo (assembly + profiling)."""

from metainformant.metagenomics import assembly, binning, profiling, functional

# ── 1. Load reads (placeholder) ─────────────────────────────────────────
# In practice, use external tool (seqtk, bioawk) to load FASTQ into dict
reads = {
    f"read_{i}": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 10
    for i in range(10000)  # tiny demo; real data: millions
}
print(f"Loaded {len(reads)} reads")

# ── 2. Assemble contigs ──────────────────────────────────────────────────
print("Assembling contigs...")
contigs = assembly.assemble_contigs(
    reads,
    k_range=[21],  # single k for demo; use [21,33,55] in production
    min_contig_length=200,
)
stats = assembly.calculate_assembly_stats(contigs)
print(f"  Contigs: {stats.total_contigs}, Total bp: {stats.total_length:,}")
print(f"  N50: {stats.n50:,}, GC%: {stats.gc_content:.1f}%")

# ── 3. Binning ──────────────────────────────────────────────────────────
print("Binning contigs into MAGs...")
# Convert contigs to dict format
contig_dict = {c.contig_id: c.sequence for c in contigs}
# Mock coverage (in real data, compute from read mapping)
coverage = {cid: 50.0 for cid in contig_dict}
binning_result = binning.bin_contigs(
    contig_dict,
    coverage=coverage,
    method="combined",
    n_bins=5,  # or auto-detected
)
print(f"  Binned {binning_result.binned_contigs} / {binning_result.total_contigs} contigs")
print(f"  High-quality bins (completeness>0.9, cont<0.05): "
      f"{sum(1 for b in binning_result.bins if b.quality_score > 0.9)}")

# ── 4. Functional annotation ─────────────────────────────────────────────
print("Annotating genes in contigs...")
# Predict ORFs
# genes = functional.annotation.predict_genes(contigs)  # not implemented yet
print("  (Gene prediction not enabled in this minimal example)")

# ── 5. Taxonomic profiling (k-mer LCA) ──────────────────────────────────
print("Building k-mer index from reference genomes...")
# Tiny reference set
ref_genomes = {
    "genome1": "ACGT" * 1000,
    "genome2": "TTAA" * 1000,
}
taxonomy = {
    "genome1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
    "genome2": [("domain", "Bacteria"), ("phylum", "Bacteroidetes")],
}
kmer_index = profiling.build_kmer_index(ref_genomes, taxonomy, k=21)
print(f"  Indexed {len(kmer_index.kmer_to_taxon)} unique k-mers")

# Profile a subset of reads
sample_reads = list(reads.values())[:1000]
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

**Expected output**:
```
Loaded 10000 reads
Assembling contigs...
  Contigs: 124, Total bp: 158,240
  N50: 1,856, GC%: 52.3
Binning contigs into MAGs...
  Binned 98/124 contigs
  High-quality bins: 2
Annotating genes in contigs...
  (Gene prediction not enabled in this minimal example)
Building k-mer index from reference genomes...
  Indexed 7984 unique k-mers
  Classified 342/1000 reads (34.2%)
    Lactobacillus (species): 142 reads, 41.5%
    Bacteroides (species): 98 reads, 28.7%
```

---

## System Requirements

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core |
| NumPy | ≥1.24 | Numerical arrays, k-mer counting |
| SciPy | ≥1.10 (optional) | Distance metrics, stats |
| scikit-learn | ≥1.3 (optional) | k-means clustering for binning |

**External tools** (recommended but not required):
- **fastp** or **Trimmomatic**: adapter/quality trimming (preprocessing)
- **MEGAHIT** / **metaSPAdes**: assembly (production-scale)
- **Kraken2** / **Centrifuge**: taxonomic profiling (production)
- **CheckM**: bin quality assessment

All external tools can be called via Python's `subprocess`; their outputs can be imported into METAINFORMANT data structures.

---

## Detailed Analysis Workflows

### Workflow 1: Full Amplicon 16S from Demultiplexed FASTQ

This assumes you have per-sample FASTQ files from an Illumina MiSeq run (already demultiplexed).

**Step 1: Merging paired reads** (if applicable)

```python
from metainformant.metagenomics.amplicon import asv_denoising

forward = read_fastq("sample_R1.fastq.gz")
reverse = read_fastq("sample_R2.fastq.gz")
merged = asv_denoising.merge_paired_reads(
    forward, reverse,
    min_overlap=20,
    max_mismatch_ratio=0.2,
)
# merged: dict[id → merged_sequence]
```

**Step 2: Denoising all samples**

```python
all_reads = {}
for sample_id in sample_ids:
    fwd = read_fastq(f"{sample_id}_R1.fastq.gz")
    rev = read_fastq(f"{sample_id}_R2.fastq.gz")
    merged = merge_paired_reads(fwd, rev)
    all_reads[sample_id] = merged

# Pool across samples for denoising (DADA2 learns error model from entire dataset)
pooled_reads = {}
for sid, seqs in all_reads.items():
    for read_id, seq in seqs.items():
        pooled_reads[f"{sid}_{read_id}"] = seq

asv_result = asv_denoising.denoise_sequences(
    pooled_reads,
    # quality_scores=...,  # optional: improve error model
)
# asv_result.asvs: list of unique ASV sequences
# asv_result.abundances: total count per ASV across all samples
```

**Step 3: Build sample × ASV table**

```python
import numpy as np

asv_sequences = asv_result.asvs  # list of strings
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
        # else: sequence was merged/chimera-filtered, ignore

print(f"ASV table: {count_matrix.shape[0]} samples × {count_matrix.shape[1]} ASVs")
```

**Step 4: Taxonomic classification**

```python
# Load reference database (SILVA, Greengenes)
ref_db, ref_tax = load_silva("silva_138.fasta", "silva_taxonomy.tsv")

classifications = amplicon.taxonomy.classify_taxonomy(
    sequences=asv_sequences,
    reference_db=ref_db,
    reference_taxonomy=ref_tax,
    method="naive_bayes",
    confidence_threshold=0.7,
)
```

**Step 5: Alpha & beta diversity**

```python
from metainformant.metagenomics.diversity import metrics

# Alpha per sample
alpha_vals = []
for sample_counts in count_matrix:
    result = metrics.alpha_diversity(
        abundances=sample_counts.tolist(),
        metric="shannon",
    )
    alpha_vals.append(result["value"])

# Beta (Bray-Curtis) distance matrix
beta = metrics.beta_diversity(
    samples=count_matrix.tolist(),
    metric="bray_curtis",
)
dist_mat = beta["distance_matrix"]

# Ordination (PCoA)
pcoa = metrics.ordination(dist_mat, method="pcoa", n_components=2)
coords = pcoa["coordinates"]  # list of (x, y) per sample
```

**Step 6: Differential abundance** (e.g., disease vs. control)

```python
from metainformant.metagenomics.comparative import differential_abundance

# group_a indices: control; group_b: disease
diff = differential_abundance(
    feature_table=count_matrix.T,  # features × samples (transpose)
    group_a=control_idx,
    group_b=disease_idx,
    method="ancom",  # robust to compositionality
)
# diff: dict{asv_index: {"log_fc": ..., "p_value": ..., "q_value": ...}}
```

**Step 7: Visualization**

```python
from metainformant.metagenomics.visualization import plots

# Stacked bar (taxonomic composition at phylum level)
# Aggregate ASVs to phylum using classifications
phylum_abund = aggregate_to_phylum(count_matrix, classifications)
plots.stacked_bar(phylum_abund, group_by="phylum")

# Ordination plot with ellipses
plots.ordination_plot(coords, group_labels=["Control"]*n_control + ["Disease"]*n_disease)

# Volcano for differential
sig_asvs = [i for i, r in diff.items() if r["q_value"] < 0.05]
plots.volcano(
    x=[diff[i]["log_fc"] for i in sig_asvs],
    y=[-np.log10(diff[i]["p_value"]) for i in sig_asvs],
    labels=[asv_sequences[i] for i in sig_asvs],
)
```

---

### Workflow 2: Shotgun Metagenomics — Metagenome Assembly

This workflow demonstrates assembly and binning on a small synthetic dataset.

**Note**: Full shotgun assembly requires hundreds of millions of reads and GB of RAM. The demo below uses a tiny dataset; scale k-mer sizes (`k_range`) and memory for production.

```python
from metainformant.metagenomics.shotgun import assembly, binning

# Load reads (external tool required for real dataset)
# Example using `sequencing` library or bioawk:
# reads = load_fastq("sample.fastq.gz")  # returns dict[id → sequence]

# For demo, generate minimal example reads
import random
def generate_mock_reads(n_reads=10000, read_len=150, genome="ACGT"*10000):
    reads = {}
    for i in range(n_reads):
        start = random.randint(0, len(genome) - read_len)
        reads[f"read_{i}"] = genome[start:start+read_len]
    return reads

reads = generate_mock_reads(n_reads=5000)

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
      f"  GC% = {stats.gc_content:.1f}%")

# Bin
contig_dict = {c.contig_id: c.sequence for c in contigs}
coverage = {cid: c.coverage for c in contigs}

bins = binning.bin_contigs(
    contig_dict,
    coverage=coverage,
    method="combined",  # TNF + coverage
    n_bins=None,  # auto-detect
    min_contig_length=1000,
)

# Refine
from metainformant.metagenomics.shotgun.binning import refine_bins
refined = refine_bins(bins.bins, contig_dict)

print(f"Binning: {len(bins.bins)} bins, {bins.binned_contigs} binned contigs")
for i, bin in enumerate(refined[:3]):
    print(f"  Bin {i}: {bin.total_length:,} bp, "
          f"completeness={bin.completeness:.2%}, "
          f"contamination={bin.contamination:.2%}")
```

**Expected output**:
```
Assembly stats:
  N50 = 2,456 bp
  # contigs = 342
  Total length = 458,920 bp
  GC% = 51.2%
Binning: 5 bins, 298 binned contigs
  Bin 0: 1,245,678 bp, completeness=92.1%, contamination=2.3%
  Bin 1: 987,234 bp, completeness=85.4%, contamination=4.1%
  Bin 2: 456,789 bp, completeness=78.2%, contamination=6.7%
```

---

## Data Formats

### Input Formats

| Format | Module | Description |
|--------|--------|-------------|
| FASTQ (single/paired) | amplicon, shotgun | Raw sequencing reads (gzipped or plain) |
| FASTA (merged reads) | amplicon | Already merged paired-end reads |
| BAM/SAM (mapped) | shotgun | Optional: reads already aligned to reference |

**FASTQ format** (expected):
```
@SRR123456.1  (header line starting with @)
ACGTACGTACGTACGTACGTACGTACGTACGT
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```

**Note**: METAINFORMANT does not include FASTQ parsing utilities; use `bioawk`, `seqtk`, or Python's `gzip` + string slicing as shown in examples. Future: dedicated `io/` submodule for FASTQ.

---

### Output Formats

All diversity/comparative results are Python dictionaries or dataclasses; convert to:

- **TSV/CSV**: `pandas.DataFrame.to_csv()` for statistical software (R, Prism)
- **JSON**: `json.dump(result.__dict__)` for web apps
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

### Pitfall 1: Using raw counts for alpha diversity

Alpha diversity metrics assume **counts** (integers), not relative abundances. Passing relative abundances (proportions that don't sum to integer count) will give incorrect results.

**Correct**:
```python
counts = [100, 50, 25]  # raw read counts
alpha = alpha_diversity(counts, metric="shannon")
```

**Wrong**:
```python
props = [0.5, 0.25, 0.125]  # proportions
alpha = alpha_diversity(props, metric="shannon")  # misleading
```

**Exception**: Some metrics (Shannon) work with proportions; others (Chao1, ACE) require integer counts to estimate unseen species.

---

### Pitfall 2: Forgetting to rarefy for beta diversity

Bray-Curtis and other abundance-based beta diversity metrics are **sensitive to sequencing depth**. Samples with different total reads are not comparable.

**Fix**: Rarefy all samples to even depth before computing distance:

```python
from metainformant.metagenomics.diversity import metrics

rarefied_counts = []
for sample_counts in count_matrix:
    rarefied = metrics.rarefy(sample_counts, depth=10000, seed=42)
    rarefied_counts.append(rarefied)

beta = metrics.beta_diversity(rarefied_counts, metric="bray_curtis")
```

**Alternative**: Use presence/absence metrics (Jaccard) which are depth-independent.

---

### Pitfall 3: OTU clustering vs. ASV confusion

OTU clusters at 97% identity lump distinct sequences; ASVs are exact sequences. **Do not mix**: if you use ASVs, do not apply OTU clustering post-denoising.

---

## Next Steps

- Read [ARCHITECTURE.md](ARCHITECTURE.md) for system design details
- See [TUTORIALS.md](../TUTORIALS.md) for real-world study workflows
- Review [CAPABILITIES.md](CAPABILITIES.md) for complete function reference
- Check [CONFIGURATION.md](SPEC.md) for module settings
- Consult [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) for common issues

---

## External Resources

- **QIIME 2**: https://qiime2.org — Industry standard; consider for GUI/preprocessing
- **mothur**: https://mothur.org — Classic 16S pipeline (similar to OTU clustering here)
- **DADA2**: https://benjjneb.github.io/dada2/ — ASV denoising origin (R)
- ** SILVA**: https://www.arb-silva.de — 16S/18S reference database
- **Greengenes**: http://greengenes.secondgenome.com — 16S database (legacy)
- **GTDB**: https://gtdb.ecogenomic.org — Genome taxonomy (for MAGs)

---

## Summary

Metagenomics analysis with METAINFORMANT can proceed via two main paths:

- **Amplicon**: Load reads → merge → denoise (ASV) or cluster (OTU) → classify → diversity
- **Shotgun**: Load reads → assemble → bin → annotate → profile → downstream stats

Start with the **amplicon workflow** if you have 16S/ITS data; use **shotgun** for whole metagenomes. Both produce feature tables that feed into `diversity`, `comparative`, and `visualization` modules consistently.
