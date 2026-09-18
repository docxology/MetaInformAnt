# Downstream Analytical Contracts (RNA campaign)

> Status: pre-freeze. The evidence manifest is not frozen; all cross-species
> analysis remains descriptive-only (see
> [METHODS_LITERATURE_ALIGNMENT.md](METHODS_LITERATURE_ALIGNMENT.md),
> "Claim boundaries and freeze gates"). The boundary statement, verbatim:
> **Descriptive-only statistics in cross-species outputs until the evidence
> manifest freezes; inferential tests (Wilcoxon/chi²-style) are gated and
> labeled for post-freeze use only.** No p-values, confidence intervals, or
> significance language may appear in cross-species outputs.

This page defines the identifiers, defaults, and replication requirements that
any consumer of RNA-campaign matrices must declare before analysis, and the
checklist that gates promotion from descriptive reporting to biological
inference. Every requirement names the source file and function that enforces
or implements it. Guard tests live in
`tests/rna/test_promotion_guards.py`.

The machine-readable contract record is
`src/metainformant/rna/analysis/statistics_contract.py`
(`AnalysisProvenance`, fail-closed validators; required by
`projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md`).

## 1. Canonical identifiers

### Sample
A **sample** is one sequencing run (SRA run accession, `SRR…`), keyed by the
`sample_id` token used by the amalgkit workspace:

- acquisition and quantification are keyed by
  `{out_dir}/quant/{sample_id}/{sample_id}_abundance.tsv`
  (`src/metainformant/rna/amalgkit/_amalgkit_impl.py` merge path-resolution
  header; `docs/rna/amalgkit/PATH_RESOLUTION.md`);
- matrix columns are sample keys after `merge`
  (`docs/rna/amalgkit/steps/07_merge.md`: "transcript_id × SRR…" columns);
- per-sample currency is decided by quantification/metadata provenance
  sidecars, not by file existence
  (`src/metainformant/rna/engine/provenance.py::is_current_quantification`,
  `::is_current_metadata`).

### Condition
A **condition** is a declared contrast level on sample metadata. The
within-species exploratory machinery aligns condition labels to matrix columns
positionally and refuses mismatched lengths
(`src/metainformant/rna/analysis/expression_analysis.py::_align_conditions_to_counts`,
`::differential_expression`). Contrast levels used in the gated comparative
layer must be declared before analysis
(`src/metainformant/rna/analysis/inferential_comparative.py::ComparativeDesign`
— `reference_level`, `treatment_level`; undeclared levels refuse).

### Replicate
The **biological replicate unit** is declared per analysis, never inferred:
`src/metainformant/rna/analysis/statistics_contract.py::AnalysisProvenance.replicate_unit`
(a required, placeholder-refusing field; rendered as
`analysis_provenance_replicate_unit`). The comparative layer additionally
requires the declared `study_col` grouping and a minimum number of
observations per study
(`src/metainformant/rna/analysis/inferential_comparative.py::ComparativeDesign.min_observations_per_study`).

### Species
A **species** is the `Genus_species` token declared by the per-species config
(`species_list` in `config/amalgkit/amalgkit_<genus_species>.yaml`; pattern
documented in [CLADE_REPLICATION.md](CLADE_REPLICATION.md) §3). Species labels
must be unique and non-placeholder wherever they index an analysis: columns of
the orthology presence table and leaves of the species tree
(`src/metainformant/rna/analysis/statistics_contract.py::validate_orthology_profile_invariants`,
`::validate_species_tree_invariants`).

### Orthology
The **orthology unit** is the orthogroup, bridged to per-species transcripts by
`src/metainformant/rna/analysis/ortholog_mapping.py::build_orthogroup_bridge`
(versioned inputs via `::OrthologySourceMetadata`, `::MappingArtifactManifest`;
per-species coverage audited by `::audit_species_retention` with threshold
`min_retention`). The canonical analysis view is the orthogroup × species 0/1
presence table (`::orthology_presence_table`), which must satisfy the
declared invariants
(`src/metainformant/rna/analysis/statistics_contract.py::validate_orthology_profile_invariants`).
Gene-to-ortholog aggregation for expression is declared per call
(`src/metainformant/rna/analysis/cross_species.py::build_ortholog_map`,
`::map_expression_to_orthologs`, `aggregation` default `"mean"`).

### Matrix
A **matrix** is a features × samples table whose axis meaning depends on the
production stage; consumers must name the stage:

| Stage | Feature axis | Source |
| --- | --- | --- |
| merge | transcripts (`{Scientific_Name}_tc.tsv` counts, `_tpm.tsv` TPM, `_eff_len.tsv` effective lengths) | `docs/rna/amalgkit/steps/07_merge.md` |
| wsfilter | transcripts after within-species filtering (exclusions recorded) | `docs/rna/amalgkit/steps/09_wsfilter.md`; step runner `src/metainformant/rna/steps.py::STEP_RUNNERS["wsfilter"]` |
| cstmm | orthologs under cross-species TMM | `docs/rna/amalgkit/steps/08_cstmm.md` |
| finalize | transcripts — the required downstream endpoint | `docs/rna/amalgkit/steps/10_finalize.md`; runner `src/metainformant/rna/steps.py::STEP_RUNNERS["finalize"]` |

Until a validated transcript-to-gene mapping exists for a species matrix,
outputs remain feature-level (claim boundary 3,
[METHODS_LITERATURE_ALIGNMENT.md](METHODS_LITERATURE_ALIGNMENT.md)). Cross-species
results must declare their role explicitly; unlabeled results are refused
(`src/metainformant/rna/analysis/statistics_contract.py::result_role`).

## 2. Normalization defaults

- **Cross-sample within-species comparison:** TPM
  (`src/metainformant/rna/analysis/expression_core.py::normalize_counts`,
  `method="tpm"`; gene lengths required). Tau computation uses log2(TPM) with
  the lowest-10%-mean rule
  (`src/metainformant/rna/analysis/tissue_specificity.py::compute_tau`,
  `log2=True`, `lowest_fraction=0.10`).
- **Differential-expression-style machinery:** raw estimated counts as input,
  not normalized values (`docs/rna/amalgkit/steps/07_merge.md`, "Counts vs
  TPM"); size-factor normalization via median of ratios is available as
  `method="median_ratio"`
  (`src/metainformant/rna/analysis/expression_core.py::estimate_size_factors`).
  `normalize_counts` defaults to `method="cpm"`; every caller must name the
  method explicitly in provenance.
- **Cross-species comparison:** TMM normalization over single-copy orthologs
  via `cstmm` (`docs/rna/amalgkit/steps/08_cstmm.md`; step runner
  `src/metainformant/rna/steps.py::STEP_RUNNERS["cstmm"]`).
- The exact `--norm`/batch settings of each finalize run must be recorded with
  the run; outputs must not be described as biologically normalized or
  batch-free unless those properties were established by the declared method
  and diagnostics (`docs/rna/amalgkit/steps/10_finalize.md`, "Completion
  evidence").

## 3. Batch handling

- The amalgkit `--batch` flag is an HPC-array indexing device, **not** a
  biological batch declaration (`docs/rna/amalgkit/steps/06_quant.md`
  parameter table; `docs/rna/amalgkit/steps/09_wsfilter.md` `--batch no`).
- Batch effects may be *detected and reported descriptively*
  (`src/metainformant/rna/analysis/qc_filtering.py::detect_batch_effects`,
  methods `kruskal`/`silhouette`/`pvca`); detection output is QC evidence and
  carries no inferential claim.
- No batch-correction claim may be made unless the correction was applied by a
  declared method with recorded diagnostics
  (`docs/rna/amalgkit/steps/10_finalize.md`).

## 4. Missing-data policy

- Compared cells must be finite: cross-species alignment refuses missing or
  non-finite values in the compared cells and rejects ambiguous duplicate
  labels
  (`src/metainformant/rna/analysis/cross_species.py::_validate_finite_expression`,
  `::_validate_expression_labels`).
- Orthology presence is explicit: absence is encoded as `0`; missing values
  and out-of-{0,1} values are refused, never coerced
  (`src/metainformant/rna/analysis/statistics_contract.py::validate_orthology_profile_invariants`).
- Observations for the comparative layer must resolve missingness at the
  caller level; missing values in design-used columns refuse fail-closed
  (`src/metainformant/rna/analysis/inferential_comparative.py::_validate_observations`).
- Tissue coverage: unmatched-tissue tau values must be masked or reported with
  their denominators (release rule 2,
  [METHODS_LITERATURE_ALIGNMENT.md](METHODS_LITERATURE_ALIGNMENT.md));
  per-tissue completeness is computed by
  `src/metainformant/rna/analysis/conservation_profiles.py::compute_per_tissue_completeness`.
- Cohort denominators (included/excluded) are declared on the provenance
  record (`AnalysisProvenance.cohort_included_count` / `.cohort_excluded_count`)
  and derived from the data by
  `src/metainformant/rna/analysis/cohort_accounting.py::build_cohort_funnel`
  (stage order in `::STAGE_NAMES`, rendering via `::render_funnel_lines`).
- A halted or impossible analysis is recorded explicitly as
  `analysis_role="stopped"`/`"unavailable"`; such a record must not declare
  any field implying results exist
  (`src/metainformant/rna/analysis/statistics_contract.py::NON_ANALYSIS_ROLES`,
  enforced in `::validate_analysis_provenance`).

## 5. Filtering defaults

- Within-species filtering is owned by `wsfilter` after `merge`; the run
  records the output table, excluded samples/features, parameters, and
  row/column counts (`docs/rna/amalgkit/steps/09_wsfilter.md`, "Validation").
- Exploratory low-expression filtering defaults: genes with mean count below
  `min_count=10` in fewer than `min_samples=2` samples
  (`src/metainformant/rna/analysis/expression_core.py::filter_low_expression`).
- Tissue-specificity filtering defaults: drop the lowest-10%-mean-expression
  features before tau
  (`src/metainformant/rna/analysis/tissue_specificity.py::filter_low_expression`).
- Cross-species fingerprint profiles require at least
  `MIN_FINGERPRINT_FEATURES = 100` valid features per species (module
  constant; `min_valid_features` parameter of
  `src/metainformant/rna/analysis/cross_species.py::compute_fingerprint_divergence_matrix`
  and `::compute_fingerprint_stability`).
- Orthology coverage defaults: every orthogroup must map to at least
  `max(min_species_per_orthogroup, ceil(min_species_fraction × n_species))`
  species, defaults `2` and `0.5`
  (`src/metainformant/rna/analysis/statistics_contract.py::validate_orthology_profile_invariants`).
- Deviations from these defaults are allowed only when recorded as declared
  parameters in the analysis provenance.

## 6. Multiplicity defaults

- **Descriptive lane (current, pre-freeze): no multiplicity procedure.**
  `multiple_testing_family`, `multiple_testing_method`, and
  `tested_feature_count` must be `None` (rendered `not-applicable`) exactly
  when `analysis_role="descriptive"`; the symmetry is enforced fail-closed in
  both directions
  (`src/metainformant/rna/analysis/statistics_contract.py::validate_analysis_provenance`).
- **Inferential lane (post-freeze, gated):** the allowed procedures are
  `bh-fdr`/`benjamini-hochberg` and `bonferroni`; an inferential record must
  declare the family, the procedure, and a positive `tested_feature_count`
  matching the supplied family exactly. The only adjustment path is
  `src/metainformant/rna/analysis/statistics_contract.py::declared_inferential_bh_fdr`,
  which re-validates the contract and refuses role, procedure, or
  family-size mismatches; the comparative layer routes its adjustment through
  it (`src/metainformant/rna/analysis/inferential_comparative.py::run_inferential_comparative_analysis`).
- The exploratory within-species adjustment helper
  (`src/metainformant/rna/analysis/expression_analysis.py::adjust_pvalues`,
  default `method="bh"`) is exploratory machinery; any inferential use of its
  outputs is out of scope until the promotion checklist below is satisfied.

## 7. Replicate requirements

- The biological replicate unit must be declared and non-placeholder
  (`AnalysisProvenance.replicate_unit`; enforced by
  `validate_analysis_provenance`).
- The comparative layer requires at least 2 studies (heterogeneity needs ≥ 2)
  and the declared per-study minimum observations
  (`ComparativeDesign.min_studies` ≥ 2,
  `.min_observations_per_study` ≥ 1; validated in
  `src/metainformant/rna/analysis/inferential_comparative.py::_validate_design`).
- Gated group comparisons require at least 2 valid observations per group
  (`src/metainformant/rna/analysis/tissue_specificity.py::wilcoxon_duplication_specificity`).
- Species-level designs require at least 2 species; the species tree needs ≥ 2
  uniquely labeled leaves and an explicit, provenance-declared rootedness
  statement — topology alone cannot establish biological rooting
  (`src/metainformant/rna/analysis/statistics_contract.py::validate_species_tree_invariants`).
- Descriptive fingerprint summaries replicate over feature resampling with the
  declared seed and resampling count
  (`cross_species.compute_fingerprint_stability`; `AnalysisProvenance.random_seed`,
  `.resampling_count`).

## 8. Promotion checklist (descriptive → biological inference)

A descriptive-stage artifact may be promoted to a biological-inference claim
only when **all** of the following hold. Each item names its enforcing source.

1. **Complete metadata.** The metadata and selection records match the exact
   current contract for every work tree the analysis consumes
   (`src/metainformant/rna/engine/provenance.py::is_current_metadata`); cohort
   denominators are declared on the record
   (`AnalysisProvenance.cohort_included_count` / `.cohort_excluded_count`).
2. **Finalized matrices.** Every consumed species workspace is finalized with
   matching feature/sample axes and passes the completion evidence checks
   (`docs/rna/amalgkit/steps/10_finalize.md`, "Completion evidence"; sidecar
   currency via
   `src/metainformant/rna/engine/provenance.py::is_current_downstream` with the
   required step set including `"finalize"`).
3. **Current provenance.** Quantification, metadata, and downstream sidecars
   are current and bound to the run
   (`src/metainformant/rna/engine/provenance.py::write_quant_provenance`,
   `::write_metadata_provenance`, `::write_downstream_provenance`,
   `::quantification_contract_id`); the analysis declares a validated,
   non-placeholder `AnalysisProvenance`
   (`src/metainformant/rna/analysis/statistics_contract.py::validate_analysis_provenance`)
   with declared artifact paths and a data-root snapshot id
   (`AnalysisProvenance.artifact_paths`, `.data_root_snapshot_id`).
4. **Review.** The metadata-harmonization review state is a real declared
   value, never a placeholder
   (`AnalysisProvenance.metadata_harmonization_review`; enforced by
   `validate_analysis_provenance` — `"pending"`, `"todo"`, `""`, etc. refuse);
   species-tree provenance (source and branch-length scale) is declared for
   tree-dependent analyses
   (`AnalysisProvenance.species_tree_source`,
   `.species_tree_branch_length_scale`).
5. **Frozen evidence manifest.** Every inferential entry point requires an
   explicit `evidence_manifest_frozen=True` affirmation at the call site;
   the default is refusal
   (`statistics_contract.declared_inferential_bh_fdr`;
   `inferential_comparative.require_inferential_contract`;
   `tissue_specificity.wilcoxon_duplication_specificity`).
6. **Declared inferential contract.** The contract is declared
   `analysis_role="inferential"` with a BH-family procedure and a
   `tested_feature_count` equal to the presented family
   (`inferential_comparative.require_inferential_contract`,
   `statistics_contract.declared_inferential_bh_fdr`).
7. **Invariant preconditions.** Orthology presence
   (`statistics_contract.validate_orthology_profile_invariants`) and species
   tree (`statistics_contract.validate_species_tree_invariants`) invariants
   pass; ortholog retention audit is at or above the declared threshold
   (`ortholog_mapping.audit_species_retention`).

Items 5–7 are enforced mechanically; items 1–4 are enforced fail-closed
wherever they are declared (placeholder or stale provenance refuses to render
or run) and are confirmed by owner review against the evidence manifest. A
descriptive-stage artifact cannot satisfy this checklist, and the public APIs
make the boundary mechanical: see
`tests/rna/test_promotion_guards.py`.

## 9. What promotion does not mean

- A green engineering gate never establishes cohort completion or biological
  inference (TODO.md preamble): passing tests and a frozen manifest license
  the *gated* APIs; they do not by themselves make a scientific claim.
- Promotion does not rewrite history: the descriptive artifacts and their
  `attrs["role"] == "descriptive"` labels remain descriptive
  (`cross_species.compute_fingerprint_divergence_matrix`,
  `cross_species.compute_fingerprint_stability`); inferential output is a new,
  separately labeled result (`role="inferential"`, `gate="post-freeze"`).
- No result-derived manuscript claim may be made until outputs, hashes, and
  run records exist and pass the evidence-manifest gate (claim boundary 5,
  [METHODS_LITERATURE_ALIGNMENT.md](METHODS_LITERATURE_ALIGNMENT.md)).
