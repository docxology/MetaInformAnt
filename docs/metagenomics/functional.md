# Functional Annotation and Pathway Reconstruction

Gene prediction, functional annotation via profile HMM search, gene family classification, metabolic pathway reconstruction, and cross-sample pathway comparison.

## Key Concepts

**ORF prediction** scans all six reading frames for open reading frames defined by start and stop codons. Supports partial ORFs at contig edges, which are common in metagenomic assemblies. Uses bacterial genetic code (NCBI Table 11) by default.

**HMM-based annotation** scores predicted protein sequences against position-specific scoring matrices (PSSMs) representing gene families from COG, KEGG, and Pfam databases. In production, HMMER profiles would replace the simplified scoring.

**Pathway reconstruction** maps annotated KO/EC identifiers to known metabolic pathways (KEGG) and scores pathway completeness as the fraction of required reactions present. Built-in definitions cover glycolysis, TCA cycle, pentose phosphate, oxidative phosphorylation, nitrogen metabolism, methane metabolism, and others.

## Data Models

### Annotation

- `ORF`: `sequence_id`, `start`, `end`, `strand`, `nucleotide_seq`, `protein_seq`, `frame`, `partial`.
- `HMMHit`: `query_id`, `target_id`, `target_name`, `score`, `e_value`, `bias`, `domain_score`, `domain_e_value`.
- `FunctionalAnnotation`: `orf_id`, `gene_families`, `kegg_orthologs`, `cog_categories`, `pfam_domains`, `ec_numbers`, `best_hit`, `confidence`.

### Pathways

- `PathwayDefinition`: `pathway_id`, `name`, `database`, `required_kos`, `required_ecs`, `optional_kos`, `reaction_steps`, `category`, `subcategory`.
- `PathwayResult`: `pathway_id`, `pathway_name`, `completeness`, `matched_kos`, `missing_kos`, `redundancy`, `confidence`.

## Function Reference

### Gene Prediction

#### `predict_orfs(sequence, sequence_id="seq", min_length=100, allow_partial=True, translation_table=11) -> List[ORF]`

Predict open reading frames in all six reading frames. Returns ORFs sorted by position with nucleotide and translated protein sequences.

### Annotation

#### `annotate_genes(sequences, hmm_db=None, e_value_threshold=1e-5, min_score=25.0) -> List[FunctionalAnnotation]`

Annotate protein sequences against an HMM database. Profile IDs prefixed with `KO:` map to KEGG, `COG:` to COG categories, and `PF` to Pfam domains. EC numbers are extracted from annotation strings.

#### `classify_gene_families(genes, database="COG", reference_profiles=None, min_score=20.0) -> Dict[str, List[str]]`

Classify genes into functional families. When no HMM profiles are available, falls back to amino acid composition heuristics.

### Pathway Reconstruction

#### `reconstruct_pathways(annotations, database="KEGG", pathway_definitions=None, min_completeness=0.0) -> List[PathwayResult]`

Map functional annotations to metabolic pathways and score completeness. Uses built-in KEGG pathway definitions when no custom set is provided. Confidence score combines completeness with gene copy redundancy.

#### `calculate_pathway_completeness(pathway, annotations) -> float`

Score a single pathway's completeness given a set of observed KO/EC annotations.

#### `compare_pathway_profiles(samples, pathway_definitions=None) -> Dict[str, Dict[str, float]]`

Compare pathway completeness across multiple samples. Returns a matrix of pathway_id -> {sample_id: completeness}.

#### `find_differential_pathways(group1_samples, group2_samples, min_diff=0.2) -> List[Dict]`

Identify pathways differentially present between two groups based on mean completeness difference. Returns pathway_id, group means, difference, and direction.

## Usage Examples

```python
from metainformant.metagenomics import (
    predict_orfs, annotate_genes, classify_gene_families,
    reconstruct_pathways, compare_pathway_profiles,
)

# Predict ORFs from a contig
contig = "ATGAAAGCGTTTCGATCGATCGATCG" * 20 + "TAA"
orfs = predict_orfs(contig, sequence_id="contig_001", min_length=50)
print(f"Predicted {len(orfs)} ORFs")

# Annotate predicted proteins
proteins = {orf.sequence_id: orf.protein_seq for orf in orfs}
annotations = annotate_genes(proteins)

# Classify into COG families
families = classify_gene_families(proteins, database="COG")

# Reconstruct metabolic pathways
ko_annotations = {"gene1": ["K00844", "K00845"], "gene2": ["K01810"]}
pathways = reconstruct_pathways(ko_annotations)
for pw in pathways[:5]:
    print(f"{pw.pathway_name}: {pw.completeness:.0%} complete")

# Cross-sample pathway comparison
samples = {"sample_a": ko_annotations, "sample_b": {"gene3": ["K01647"]}}
profiles = compare_pathway_profiles(samples)
```

## Built-in KEGG Pathways

| ID | Name | Required KOs |
|----|------|-------------|
| map00010 | Glycolysis / Gluconeogenesis | 11 |
| map00020 | TCA cycle | 12 |
| map00030 | Pentose phosphate pathway | 7 |
| map00190 | Oxidative phosphorylation | 14 |
| map00220 | Arginine biosynthesis | 8 |
| map00910 | Nitrogen metabolism | 9 |
| map00680 | Methane metabolism | 8 |
| map00195 | Photosynthesis | 8 |

## Configuration

Environment variable prefix: `META_`

## Related Modules

- `metainformant.metagenomics.amplicon` -- OTU/ASV input for functional profiling
- `metainformant.metagenomics.shotgun` -- assembly and binning
- `metainformant.metagenomics.diversity` -- community diversity metrics
- `metainformant.metagenomics.comparative` -- differential analysis
