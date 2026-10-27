# Functional Annotation and Pathway Reconstruction

Gene prediction, functional annotation by scoring against position-specific profiles (a simplified stand-in for profile HMM search), gene family classification, metabolic pathway reconstruction, and cross-sample pathway comparison.

## Key Concepts

**ORF prediction** scans all six reading frames for open reading frames defined by start and stop codons. Supports partial ORFs at contig edges, which are common in metagenomic assemblies. `translation_table` selects the NCBI genetic code (default 11, bacterial; `ValueError` for an id that is not an NCBI table). The table sets the amino-acid assignments and the stop codons (TGA codes for tryptophan, not a stop, in table 4) and, through its own start codons, which of ATG, GTG and TTG can open an ORF: all three for tables 11 and 4, ATG and TTG for table 1. The rare alternative starts (CTG, ATT, ATC, ATA) are not used. Every in-frame start codon upstream of a stop produces its own ORF if that ORF reaches `min_length`, so ORFs that share a stop codon are nested and overlap. With `allow_partial=True`, position 0 of frame 0 also counts as an implicit start, and each (strand, start, end) is still reported once even when the sequence begins with a start codon.

**HMM-based annotation** scores predicted protein sequences against position-specific scoring matrices (PSSMs) representing gene families from COG, KEGG, and Pfam databases. No profile database is bundled: pass your own as `hmm_db` (see `annotate_genes`). In production, HMMER profiles would replace the simplified ungapped scoring. The scoring matrices have no calibrated null model, so no E-value or composition bias is computed: hit `e_value`, `bias` and `domain_e_value` are NaN, and hits are filtered by `min_score` alone.

**Pathway reconstruction** maps annotated KO/EC identifiers to known metabolic pathways (KEGG) and scores pathway completeness. A definition without `reaction_steps` is scored as the fraction of its `required_kos` and `required_ecs` that are present. A definition with `reaction_steps` (a list of alternatives per step) is scored by steps, in `reconstruct_pathways` and in `calculate_pathway_completeness` alike: a step is present when any of its alternatives is annotated, and completeness is the fraction of steps present, so it never exceeds 1 (`required_kos` and `required_ecs` are then informational only). Built-in definitions (glycolysis, TCA cycle, pentose phosphate, oxidative phosphorylation, nitrogen metabolism, methane metabolism, and others) have no reaction steps.

## Data Models

### Annotation

- `ORF`: `sequence_id`, `start`, `end`, `strand`, `nucleotide_seq`, `protein_seq`, `frame`, `partial`, plus the properties `length_nt` and `length_aa`.
- `HMMHit`: `query_id`, `target_id`, `target_name`, `score`, `e_value`, `bias`, `domain_score`, `domain_e_value`, `query_start`, `query_end`, `target_start`, `target_end`, `description`. `e_value`, `bias` and `domain_e_value` are NaN (not computed); `query_start`, `query_end`, `target_start` and `target_end` are placeholders (the whole query, no target coordinates), not alignment bounds.
- `FunctionalAnnotation`: `orf_id`, `gene_families`, `kegg_orthologs`, `cog_categories`, `pfam_domains`, `ec_numbers`, `best_hit`, `description`, `confidence`.

### Pathways

- `PathwayDefinition`: `pathway_id`, `name`, `description`, `database`, `required_kos`, `required_ecs`, `optional_kos`, `reaction_steps`, `category`, `subcategory`, plus the property `total_required`.
- `PathwayResult`: `pathway_id`, `pathway_name`, `completeness`, `matched_kos`, `matched_ecs`, `missing_kos`, `missing_ecs`, `total_required`, `total_matched`, `redundancy`, `confidence`.

## Function Reference

The `metainformant.metagenomics` package exposes only submodules. The functions live in `metainformant.metagenomics.functional.annotation` (`predict_orfs`, `annotate_genes`, `classify_gene_families`) and `metainformant.metagenomics.functional.pathways` (`reconstruct_pathways`, `calculate_pathway_completeness`, `compare_pathway_profiles`, `find_differential_pathways`).

### Gene Prediction

#### `predict_orfs(sequence, sequence_id="seq", min_length=100, allow_partial=True, translation_table=11) -> List[ORF]`

Predict open reading frames in all six reading frames. `min_length` is in nucleotides. Returns ORFs sorted by position with nucleotide and translated protein sequences; reverse-strand ORFs report forward-strand coordinates. All ORFs from one call share `sequence_id`, so build unique keys (for example from the index or the coordinates) before using them as gene IDs. Raises `ValueError` if `translation_table` is not an NCBI genetic code id.

### Annotation

#### `annotate_genes(sequences, hmm_db=None, e_value_threshold=None, min_score=25.0) -> List[FunctionalAnnotation]`

Annotate protein sequences (`dict` of gene ID -> protein sequence) against an HMM database, given as `{profile_id: {position: {amino_acid: score}}}`. A profile hits when its best ungapped alignment score, averaged over the matched positions, reaches `min_score`. Profile IDs prefixed with `KO:` map to KEGG (`kegg_orthologs`, prefix removed), `COG:` to COG categories (`cog_categories`, prefix removed), and `PF` to Pfam domains (`pfam_domains`); any other ID goes to `gene_families`. EC numbers are extracted from `gene_families` entries that contain `EC:<a.b.c.d>`. Without `hmm_db` every annotation is empty. `e_value_threshold` must stay `None`: E-values cannot be computed for position-specific scoring matrices without a calibrated null model, so a numeric threshold raises `ValueError`, and hit E-values are NaN. `confidence` is `min(1, best_score / 100)`.

#### `classify_gene_families(genes, database="COG", reference_profiles=None, min_score=20.0) -> Dict[str, List[str]]`

Classify genes (`dict` of gene ID -> protein sequence) into functional families, where `database` is `"COG"`, `"KEGG"` or `"Pfam"`. With `reference_profiles` (same format as `hmm_db`; only profiles whose ID starts with the database prefix are used) it returns the IDs of all profiles scoring at least `min_score` (an empty or non-matching profile set gives empty lists). When `reference_profiles` is None no family is assigned: every gene maps to an empty list (unannotated), because a family can only be assigned by comparison with reference profiles and amino acid composition cannot identify one.

### Pathway Reconstruction

#### `reconstruct_pathways(annotations, database="KEGG", pathway_definitions=None, min_completeness=0.0) -> List[PathwayResult]`

Map functional annotations to metabolic pathways and score completeness. `annotations` is a dict of gene ID -> list of KO (`K` plus five digits) and EC identifiers, not a list of `FunctionalAnnotation` objects. Uses built-in KEGG pathway definitions (`database` only filters these; no MetaCyc pathways are built in) when no custom set is provided. Returns the pathways sorted by completeness, descending. `completeness` is `total_matched / total_required`: for a definition without `reaction_steps` these are the numbers of required KOs and ECs that are present and required; for a definition with `reaction_steps` they are the numbers of steps with at least one annotated alternative and of steps. It therefore lies in [0, 1]. For step-based definitions `matched_kos` and `matched_ecs` list the annotated alternatives of the present steps and `missing_kos` and `missing_ecs` list every alternative of the absent steps. `redundancy` is the mean number of gene copies per matched KO, and `confidence` is `completeness * min(1, redundancy)`.

#### `calculate_pathway_completeness(pathway, annotations) -> float`

Score a single pathway's completeness. `pathway` is a `PathwayDefinition` object and `annotations` is a set of observed KO/EC identifiers. If the definition has `reaction_steps`, a step counts as complete when any of its alternatives is present; otherwise the fraction of `required_kos` and `required_ecs` present is returned.

#### `compare_pathway_profiles(samples, pathway_definitions=None, database="KEGG") -> Dict[str, Dict[str, float]]`

Compare pathway completeness across multiple samples (`samples` is a dict of sample ID -> the `annotations` dict taken by `reconstruct_pathways`). Returns a matrix of pathway_id -> {sample_id: completeness}.

#### `find_differential_pathways(group1_samples, group2_samples, pathway_definitions=None, min_diff=0.2) -> List[Dict]`

Identify pathways differentially present between two groups based on mean completeness difference (no statistical test). Each group is a dict of sample ID -> annotations dict. Returns dicts with `pathway_id`, `pathway_name`, `group1_mean`, `group2_mean`, `difference`, `abs_difference`, and `direction` (`"group1_enriched"` or `"group2_enriched"`), sorted by `abs_difference`.

## Usage Examples

```python
from metainformant.metagenomics.functional.annotation import (
    annotate_genes, classify_gene_families, predict_orfs,
)
from metainformant.metagenomics.functional.pathways import (
    compare_pathway_profiles, reconstruct_pathways,
)

# Predict ORFs from a contig
contig = "ATGAAAGCGTTTCGATCGATCGATCG" * 20 + "TAA"
orfs = predict_orfs(contig, sequence_id="contig_001", min_length=50)
print(f"Predicted {len(orfs)} ORFs")

# Annotate predicted proteins (ORFs from one contig share sequence_id, so key by index)
proteins = {f"{orf.sequence_id}_orf{i}": orf.protein_seq for i, orf in enumerate(orfs)}
annotations = annotate_genes(proteins)  # no hmm_db given: every annotation is empty

# Without reference_profiles no family is assigned: every gene maps to an empty list
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
| map00020 | Citrate cycle (TCA cycle) | 12 |
| map00030 | Pentose phosphate pathway | 7 |
| map00190 | Oxidative phosphorylation | 14 |
| map00220 | Arginine biosynthesis | 8 |
| map00910 | Nitrogen metabolism | 9 |
| map00680 | Methane metabolism | 8 |
| map00195 | Photosynthesis | 8 |

## Configuration

This module reads no environment variables (`annotation.py` and `pathways.py` define `_ENV_PREFIX = "META_"`, but no function uses it).

## Related Modules

- `metainformant.metagenomics.amplicon` -- OTU/ASV input for functional profiling
- `metainformant.metagenomics.shotgun` -- assembly and binning
- `metainformant.metagenomics.diversity` -- community diversity metrics
- `metainformant.metagenomics.comparative` -- differential analysis
