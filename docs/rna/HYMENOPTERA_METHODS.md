# Hymenoptera methods and analysis boundaries

The 27-species campaign uses a frozen inventory, configuration, reference index,
and quantification contract. Keep these bindings stable while processing;
changing a reference or quantifier creates a new analysis cohort.

## Tool contract

Amalgkit is pinned to **0.16.60**, tag `v0.16.60`, source revision
`c656a52aacdcee6fd3bf7e8031769ca957204ebc`. This is the latest upstream release
verified on 2026-10-05 ([release](https://github.com/kfuku52/amalgkit/releases/tag/v0.16.60)).
The AWS worker also verifies Kallisto 0.52.0 and SRA Toolkit 3.4.1. The current downstream
Amalgkit stages are Python implementations; legacy R installation helpers
are not prerequisites for this pinned workflow. A newer
upstream development commit is not silently substituted into the campaign.

## Resolving archive metadata gaps

An optional `source_resolutions.json` supplement binds the original inventory
hash to preserved NCBI experiment XML. It requires matching run identity,
taxonomy, RNA-Seq strategy, public loaded data, positive spot/base counts, and
a primary public SRA file. Source-file HEAD checks confirm the declared sizes.
Hardened XML parsing rejects DTDs, entity expansion, external entities, and
malformed evidence. Resolution cannot add samples or alter a reference/configuration. Worker
partitions overlay recovered counts into a new metadata file and record its
hash, the original metadata hash, and the source-evidence hash.

For these runs, scheduling uses a modeled raw-file reservation of SRA bytes
plus three times base count plus 512 bytes per spot. This is a conservative
resource model, not measured FASTQ size; free-space floors, admission limits,
and job deadlines still apply. The original frozen ENA records remain intact.

## Numerical input and result contracts

- CPM, TPM, RPKM, quantile, and median-ratio normalization reject non-finite
  or negative counts. TPM/RPKM reject infinite or nonpositive supplied lengths.
  The existing missing-length median replacement remains a documented
  approximation; it is not measured transcript length.
- The standard median-ratio estimator uses genes positive in every sample.
  When no such gene exists, the existing library-size fallback remains an
  approximation, not equivalence to the full DESeq2 estimator.
- Differential-expression methods require at least two samples per group,
  finite nonnegative counts, and unique feature/sample IDs. Welch's test with
  one constant group uses the variable group's degrees of freedom.
- BH/Bonferroni correction rejects infinite or out-of-range probabilities.
  NaN identifies an unscoreable test and is retained; it occupies a place in
  the declared family, conservatively counted as one during adjustment.
- PCA rejects infinities, duplicated axes, invalid dimensions, singleton
  samples, and matrices without between-sample variation. Missing-value
  mean imputation requires `missing="mean"` and is recorded in the result.
  Oversized component requests retain the existing dimension cap. PC signs
  use a deterministic loading anchor; reconstruction and explained variance
  are checked against independent numerical implementations.
- `WithinSpeciesOrchestrator` requires an explicit `run` or `sample` metadata
  column, complete expression-sample coverage, and unique metadata IDs.
  Metadata may include additional unquantified runs. Its PCA defaults to
  log2 CPM; choose `normalization="log2"` for already normalized inputs.
  PCA and exploratory DE sidecars bind input hashes, retained IDs, and methods.
  Calculation failures propagate from `run_all()`.

## Statistical interpretation

Run accessions do not prove independent biological replication. The convenience
Welch comparison is exploratory. It does not adjust for study, batch,
biological-individual identity, pairing, confounding, or phylogeny. The
`deseq2_like` utility is not the DESeq2 package and does not supply its complete
empirical-Bayes inference pipeline.

Use the existing `statistics_contract`, `inferential_comparative`, and
`phylogenetic_comparative` interfaces for declared estimands, replicate units,
multiplicity families, sensitivity analyses, and validated species trees.
Gene-level comparisons require a verified transcript-to-gene/orthology bridge;
feature fingerprints alone cannot establish gene conservation or phylogeny.

## Figures

Native cross-species fingerprints retain their descriptive role, 0–2 distance
scale, and average-linkage clustering. Feature-resampling sensitivity plots
show recorded interval endpoints and the original estimate separately: a point
outside its percentile interval must not extend the interval. Invalid interval
bounds or non-finite values are rejected before figure output.

Operational coverage plots describe the frozen eligible denominator and S3
receipt inventory at one checkpoint. They are not final expression results.
Real-data diagnostic PCA uses a named subset with explicit input/selection
provenance; it cannot stand in for the completed cohort.

## Validation and downstream release

Regression tests include invalid-input controls, the unreplicated `p=0`
counterexample, the sparse size-factor counterexample, the unequal-group Welch
oracle, SciPy BH correction, scikit-learn PCA variance, SVD reconstruction,
real-file metadata boundaries, and real matplotlib interval geometry.

Quant locking is the first gate. Then restore and verify every eligible receipt,
run `merge → wsfilter → finalize → sanity` for the same snapshot, record downstream
provenance, harmonize metadata, validate orthology/design/tree inputs, regenerate
figures and captions, and require the project evidence manifest without
`--allow-missing`. Until those gates pass, complete-cohort statistics and
manuscript release remain unavailable.

See [durable quantification](DURABLE_QUANT.md),
[analysis APIs](../../src/metainformant/rna/analysis/README.md), and the
[project](../../projects/hymenoptera_amalgkit/README.md).

## Thin project orchestration

The project delegates finalized-matrix validation and mean profiles to
`expression_io`, membership/retention/duplicate audits to `ortholog_diagnostics`,
mean-expression orthogroup distances to `ortholog_profiles`, descriptive Wilson
intervals to `counting_statistics`, and JSON contract parsing to `statistics_io`.
The numerical definitions remain in the parent package; cohort YAMLs, selected-root
paths, manifests, captions, command arguments and project evidence decisions remain
in the nested repository. Its method-adapter validator resolves actual parent
callables and rejects shadowing implementations, unavailable symbols and broken
Markdown fragments.

Mapped mean-profile distances correlate expression across shared orthogroups,
selecting the first recorded transcript per mapping cell. Unsupported or constant
pairs remain NaN with their overlap counts; malformed declared matrices fail.
This is distinct from correlating each gene across aligned samples and from
comparing native feature fingerprints.

See the [project ownership map](../../projects/hymenoptera_amalgkit/doc/02_workflow/04_metainformant_methods.md)
and its [dated fleet checkpoint](../../projects/hymenoptera_amalgkit/doc/01_infrastructure/completion_checkpoint_20261006.json).
The earlier coverage snapshots retain their observation times and counts. Six-worker
operation and a 64-vCPU quota are capacity facts, not evidence of sustained linear
throughput or completed inference.
