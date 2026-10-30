# Changelog

All notable changes to METAINFORMANT are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Physical-file verification snapshots (`core.io.verification.BoundReads`) and
  conservative recovery process observations (`rna.engine.process_observation`),
  with real file-change, PID-identity, denied-inspection and timeout controls.
- Local RNA acquisition frontier (`metainformant.rna.engine.acquisition_frontier`,
  opt-in via `AMALGKIT_LOCAL_FRONTIER_SIZE`): separate acquisition, recovery and
  quantification pools bounded by a planned-identity frontier and frozen
  raw-size reservations, an fsynced transition ledger, STOP-aware draining, and
  reuse of current-contract quantifications while acquisition slots are busy.
  Supported by per-sample acquisition admission
  (`AMALGKIT_PIPELINE_ACQUISITION_SLOTS`), size-band ordering
  (`acquisition_order`), verified ENA staging and promotion
  (`rna.retrieval.ena_staging`), shared ENA connection budgets
  (`rna.retrieval.ena_transport`), strict FASTQ record validation, and bounded
  ENA metadata retries inside the first-pass deadline.
- Read-only campaign observation (`rna.engine.campaign_observation`,
  `campaign_observation_report`, `scripts/rna/observe_acquisition.py`):
  content-bound capture of a running campaign's manifest and journals,
  restore-bound new-output counts, and an isolated report bundle with an
  acquisition figure (`visualization.plots.rna_campaign`).
- `metainformant.rna.analysis.statistics_contract`: predeclared, fail-closed
  analysis-provenance records (descriptive/inferential role separation,
  role-conditional multiple-testing semantics), BH-FDR helper, gated
  inferential wrapper, orthology profile invariants, and species-tree
  invariants with caller-declared rootedness provenance. The hymenoptera
  cross-species runner validates provenance before writing any artifact.
  Extended additively: optional fail-closed-when-declared bindings for the
  data-root snapshot id, cohort inclusion/exclusion denominators, per-artifact
  paths, metadata-harmonization review state, and species-tree source plus
  branch-length scale; new non-analysis roles `stopped`/`unavailable`
  (statistical_analysis_plan sections 3, 4, 8, 9) that cannot carry
  result-implying fields, rendered as additive `analysis_provenance_*`
  summary lines.
- `metainformant.rna.engine.preflight`: mandatory start-of-run campaign
  preflight. `StreamingPipelineOrchestrator.run_all()` now probes data-root
  write access (create/rename/stat/unlink) and bare `amalgkit` PATH resolution
  before any discovery or scheduling work, refusing to start into a broken
  environment (observed 2026-09-03: a producer without external-volume write
  access marked 7,291 tasks failed; another failed every quant batch on a
  missing `amalgkit` binary). Also available standalone via
  `python -m metainformant.rna.engine.preflight`.
- `metainformant.rna.engine.progress_db.classify_sample_error`: durable
  failure classes for stored sample errors, separating environmental damage
  (write denial, missing tool) from transfer and quantification failures;
  surfaced as `db_failure_classes` in the hymenoptera campaign status report.
- `metainformant.rna.analysis.inferential_comparative`: gated predeclared
  inferential comparative analysis — study-aware OLS contrast fits,
  DerSimonian-Laird random-effects heterogeneity summaries, bootstrap
  confidence intervals, and BH-FDR exclusively via the gated
  `declared_inferential_bh_fdr` wrapper; registered sensitivity runner
  (`leave_one_study_out`, `exclude_covariate`) with directional agreement
  checks.
- `dna.population.normalized_fay_wu_h` (Zeng et al. 2006) and an optional
  `outgroup` for `calculate_fay_wu_h` and `neutrality_test_suite`.
- `metainformant.rna.analysis.ortholog_mapping` extended additively: a
  fail-closed per-species retention audit (`audit_species_retention`,
  `DEFAULT_MIN_RETENTION`) and the versioned `MappingArtifactManifest` v1
  with fail-closed write/read round-trips. The existing behavior of the
  retention audit columns is unchanged.
- Documentation call and link gates (`quality.doc_verification`):
  `find_stale_doc_calls` binds every quoted call on an imported
  `metainformant` name against the callee's current signature, and
  `find_broken_doc_links` requires relative links and images to resolve with
  exact case and anchors into Markdown files to match a GitHub heading slug or
  an explicit id. Both run in `scripts/verify_documentation_code.py`, and the
  repository is held at zero findings. `generate_core_api_reference.py --check`
  reports a stale core API reference.

### Changed
- RNA acquisition anchors its first-pass service deadline at admission, before
  synchronous logging. Admission queue time remains outside that budget.
- Unknown pharmacogenomic alleles now produce Indeterminate calls, unsupported
  genes raise instead of borrowing CYP2D6 bands, and missing ACMG frequencies
  cannot supply rarity evidence. Unscored diplotypes retain a missing score.
- NCBI genome adapters use Datasets v2 report, sequence-report and download
  endpoints, with explicit request errors, pagination and annotation selection.
  Package downloads preserve existing archives on failed or interrupted transfers;
  repeated assembly accessions across report pages raise an integrity error.
- DDM likelihoods use a first-passage density; sequence evolution honours its
  mutation rate; population simulation adds an ancestral recombination graph
  and demographic sampling rather than heuristic summary substitutions.
- Multi-trait colocalization enumerates causal-sharing configurations in log
  space, with explicit prior and study-model assumptions.
- Long-read workflows share SAM-compatible MM/ML decoding, preserve unknown
  modification probabilities, reject stale MN sequence lengths, and retain phasing keys through JSON.
- MGF parsing rejects truncated records; multiomics CCA shares a ridge core
  and rejects unregularized rank degeneracy; Cox and log-rank calculations
  handle ties and group covariance explicitly.
- Information metrics respect logarithm units, enrichment tests validate
  their background and correction choices, and failed model evaluations
  remain failures instead of zero scores. Association covariance uses the design
  pseudoinverse directly, and pure-Python R-squared rejects one-sample folds.
- Single-cell and spatial methods validate graph weights, diffusion spectra,
  missing coordinates and undefined statistics. Requested algorithms are no
  longer silently replaced by unrelated fallbacks.
- Structural-variant processing uses both translocation breakends; protein
  structure analysis uses explicit contact, surface and alignment algorithms.
  Missing structure fields remain missing.
- Domain documentation and the active backlog now describe these contracts;
  unresolved provenance, approximation and live-service work remains open.
- The blocking dependency audit covers every package in `uv.lock` (all extras,
  dependency groups and environment markers): `security_checks.sh` exports a
  PEP 751 `pylock.toml` and audits it in place against OSV, which carries both
  GitHub and PyPA advisories. It previously scanned only the default export
  through a temporary pip environment and missed dev-extra advisories such as
  nbconvert's. A locked dependency that cannot be audited blocks unless
  baselined.
- Repository-wide typing completion: `mypy --config-file pyproject.toml
  src/metainformant` reports 0 errors in 750 files under mypy 2.4.0, both with
  every extra installed and in CI's dev-only environment (404 errors in 121
  files before this release's final typing pass). Includes third-party import
  overrides for unstubbed libraries, removal of stale `# type: ignore`
  comments, and typing fixes across every domain package. CI and
  `uv_quality.sh` now run the zero-error budget gate
  (`scripts/quality/check_mypy_budget.py`) over all of `src/metainformant`
  instead of only `rna`.
- The blocking full-suite coverage floor rises from 58% to 70% (hosted CI
  measures 71.02%).
- Documentation verification checks every quoted `metainformant` import and
  every attribute chain on an imported module or class against the installed
  package (`quality.doc_verification.find_stale_doc_imports`,
  `find_stale_doc_attributes`, also run by `scripts/verify_documentation_code.py`).
  Documentation across the domains was rewritten to the real APIs; examples
  for functions that do not exist were removed rather than kept as
  aspirational text.
- `TODO.md` keeps only open work.
- Scripts configure logging with `configure_logging_from_env()`, so they
  report progress at INFO by default; `CORE_LOG_LEVEL` overrides it.
- Interface changes that come with the fixes below:
  `calculate_fay_wu_h` returns a `None` p-value without an outgroup;
  `two_epoch_effective_size` now means the coalescent effective size;
  `generate_two_populations(n_sites=...)` counts polymorphic sites and
  defaults to Watterson's expectation; `canonical_correlation` no longer
  accepts scikit-learn keyword arguments; `cross_validate_biological` raises
  on metric failures; degenerate linear GWAS fits return `status: "error"`;
  `SurvivalPredictor.predict` returns restricted mean survival time;
  `EventSequencePredictor` no longer lists `lstm`;
  `generate_synthetic_life_events` raises `TypeError` for unsupported
  options; and the AntWiki species index raises `NotImplementedError`.
- `interpret_neutrality_results` (and so `neutrality_test_suite`'s
  `neutrality_summary`) labels Tajima's D or Fu and Li's D* as selection or
  demography only when its p-value is below `alpha` (default 0.05);
  non-significant statistics read `neutral_evolution`. Results without
  p-values keep the sign-based labels.
- `classify_orthogroups` and `join_expression_with_orthology` vectorized
  (~8.9x and ~12.6x on the campaign-scale benchmark, outputs bit-identical).
- Cross-species figures are provenance-honest: descriptive-only annotations,
  data-derived species denominators, explicit not-a-species-tree labeling.
- Static boundary gate: the MCP tool layer is a recognized adapter surface;
  eqtl/popgen orchestration modules moved under sanctioned `workflow/`
  packages; math/simulation plot helpers moved into `visualization/analysis/`.
- Subprocess test timeouts raised for cold-import headroom on contended
  external-drive checkouts (examples runner, life-events CLI).
- The full-suite coverage floor rises from 70.00% to 74.90%, about one point
  under the 75.91% hosted measurement (Test Suite run 40594851520).
- Interface changes from the defect fixes below (callers may notice):
  - `generate_clinical_report(drugs=None)` analyses every CPIC drug (`drugs=[]`
    requests none); `load_allele_definitions` defaults `source=None` and
    rejects contradictory arguments.
  - Metabolomics `read_csv` keeps missing cells as NaN and raises on other
    non-numeric text; normalization refuses NaN; `ppm_tolerance` must be
    positive; unknown `ion_mode` and undersized `background_size` raise.
    `correct_batch_combat` raises on single-sample batches.
  - Metagenomics bin quality is `None` without marker annotations,
    `refine_bins` requires them, numeric `e_value_threshold` raises, Naive
    Bayes abstains below its confidence threshold, and chao1, ACE and Fisher's
    alpha reject proportions.
  - Epigenome peak files store p and q as `-log10` (probabilities in memory);
    `find_tf_binding_sites` and `identify_cpg_islands` require sequence.
  - Single-cell `leiden_clustering` raises without leidenalg unless
    `allow_kmeans_fallback=True`; `select_hvgs`, `compute_neighbors` and
    `compute_diffusion_map` return copies unless `inplace=True`; the `run_*`
    external-tool helpers raise `NotImplementedError`. Xenium and MERFISH
    loaders raise when coordinates are missing, and Visium raises when no
    barcode matches. `compute_ligand_receptor_interactions` needs gene names,
    and every `ripley_k` value moves with exact edge correction.
  - Translocation calls carry `chrom2`/`pos2`; insertions are exempt from size
    filtering; `survivor_merge` combines genotypes and qualities.
  - The long-read adapter table drops its two 8 bp entries; custom adapters
    shorter than the false-positive minimum are matched exactly.
  - Population sequence generators are coalescent: same-seed outputs differ
    from before, and requests too dense for infinite sites raise. Negative RNA
    fold changes divide.
  - `effective_size_sex_ratio` no longer adjusts by `sex_ratio`;
    `fisher_exact_test` and `compare_statistics` raise instead of returning
    p = 1.0; `calculate_confidence_intervals` estimates intervals for a
    dictionary of observation lists and rejects point values.
  - `granger` and `grnboost2` are no longer advertised by GRN inference;
    `enrich_genes` adds a `significant` column; `random_state=None` means fresh
    entropy throughout networks, ontology, ML and information;
    `OllamaConfig.max_retries` counts attempts and must be at least 1.
  - Failed GWAS and eQTL tests carry NaN statistics, are written as `NA` and
    are excluded from multiple-testing correction; `download_variant_data`
    supports only `summary_stats`.
  - AlphaFold helpers that need the remote service raise
    `NotImplementedError`; `parse_alphafold_confidence` returns one pLDDT per
    residue. Quality metrics raise for unsupported data types.

### Removed
- The generated Cursor skill wrappers (`.cursor/`), their generator
  (`scripts/package/generate_cursor_skills.py`), its parity test and
  `docs/CURSOR_SKILLS.md`. The `cursorrules/` guidelines remain.
- Lowercase duplicates of `docs/SETUP.md`, `docs/agents/SPEC.md`,
  `docs/simulation/AGENTS.md` and `docs/spatial/INTEGRATION.md`. Each was
  byte-identical to its uppercase twin, and on case-insensitive filesystems
  the pair collided.
- `scripts/rna/discover_species.py` and `scripts/rna/setup_genome.py`, which
  needed functionality the package no longer provides. The genome guides now
  document `rna.amalgkit.genome_prep` directly.

### Fixed
- `find_orfs` (`dna.expression.translation`, `dna.sequence.core`) matched start
  codons at any offset. Every ORF was reported once per frame of its strand,
  up to three times, with only one frame label correct. Starts are now
  codon-aligned within each frame.
- 41 broken documentation links. All 30 index links in the generated core API
  reference were dead on GitHub because it used `{#id}` heading attributes.
  Nine anchors left out GitHub's double hyphen for `&` and `+`, and two links
  pointed at a missing directory and at untracked campaign output.
- Pharmacogenomic phenotypes contradicted CPIC: `classify_phenotype` was
  wrong on 32 of 78 CPIC diplotypes (CYP2D6 `*1/*4` read normal, CYP2C19
  `*1/*1` rapid) and gaps in the bands read as ultrarapid. Allele values now
  match CPIC's allele functionality records (CYP2D6 `*9`, `*14`, `*41`; CYP2C9
  `*5`, `*8`; SLCO1B1 `*5`, `*15`; UGT1A1 `*37`), and the metabolism facade no
  longer carries a second, conflicting table.
- Metabolomics: unparseable cells became 0.0; enrichment returned p = 0.0
  when `background_size` was too small; ComBat depended on batch order and
  shrank by a constant. ComBat is now parametric empirical Bayes (Johnson et
  al. 2007).
- Metagenomics invented results: bins scored complete with no marker search,
  Naive Bayes classified random sequence at confidence 1.0, and functional
  annotation invented COG/KEGG/Pfam labels. The ASV Poisson tail lost
  precision and merged real variants; Fisher's alpha, pathway completeness and
  the OTU prefilter were wrong.
- Epigenome peak files did not round-trip, overlap detection lost 28% of
  nested pairs, enrichment overflowed, and `find_tf_binding_sites`, CpG-island
  detection and the ChIP/ATAC integration were simulated or placeholders.
- Single-cell and spatial: diffusion pseudotime was Euclidean distance (now
  DPT, matching scanpy), size factors were inverted, three steps mutated
  their input, the radius graph never built, mutual kNN weights doubled,
  Xenium and MERFISH plots swapped axes, and loaders invented `[0, 0]`
  coordinates. The `singlecell` extra required the wrong Louvain package.
  Ligand-receptor scoring hashed gene symbols onto columns, deconvolution
  "marker correlations" were a variance proxy, `ripley_k` skipped edge
  correction at 100 points or fewer, `seurat_v3` HVG selection was advertised
  but rejected, and domain selection called an elbow heuristic a gap statistic.
- Structural variants: translocations used breakpoints on two chromosomes as
  one interval, insertions never passed the size filter, breakpoint
  refinement and gene-overlap queries missed data, `dosage_gain` was
  unreachable, and CBS took 37 s for 8,000 bins (now 0.22 s, identical
  segments).
- Long reads: 8 bp adapters trimmed about 200 bp from adapter-free reads and
  split random reads as chimeras; haplotagging tagged every read HP=0; the
  "bayesian" methylation model was capped at 2/3; pysam ML tags crashed.
- Population genetics: the sequence generators used a star genealogy and
  founder resampling, so neutral samples had Tajima's D near -2.6 and the
  bottleneck and expansion signatures were reversed. They are now a seeded
  Kingman coalescent, and the neutrality tests offer simulated nulls. The LD
  generator now hits its r² target.
- Math and core: effective sizes asserted on integers and jumped at a sex
  ratio of 0.5; the SciPy-free Fisher test, LD-decay rate and DDM timing were
  rough approximations; `compare_statistics` ran Welch's t-test for any test
  type; dictionary confidence intervals were an invented ±10% band; the
  downloader double-counted bytes and failed on complete files.
- Networks, ontology, ML and information: a fabricated pathway p-value, a GSEA
  null that permuted weights instead of gene sets, Stirling log-factorials,
  lost ontology hierarchies, an approximate "Von Neumann" entropy, broken
  default validation strategies, and global random seeding.
- Visualization and GWAS plots: keyword collisions, inflated genome-wide
  positions (34 s for 20,000 variants, now 0.27 s), random Q-Q expected values,
  inverted sequence logos, and plots that fabricated annotation categories, a
  heritability and p-values.
- GWAS analysis, phenotype and protein: failed association tests were
  reported as beta 0 and p = 1; pipeline outputs lacked MAF; the cis-eQTL
  window was ignored; cross-omic p-values were approximations and Spearman
  ignored ties; protein motif, pocket and AlphaFold functions were
  placeholders; an eQTL script labelled synthetic positions as real. The
  mixed model caught a monomorphic SNP only when the BLAS rejected the
  nearly singular system, so on Linux it reported such a SNP as tested
  (p near 1); a rank check now rejects it on every platform.

- Neutrality tests (`dna.population`) were wrong. Tajima's D subtracted the
  whole-locus Watterson's theta from the per-site nucleotide diversity, so a
  four-sequence example gave -6.24 instead of 0.59, and the variance was
  floored at 0.01. Fu and Li's D* and F* were unstandardized counts with no
  variance. The Fay and Wu's H p-value divided a per-site H by a per-locus
  standard deviation and was always near 1. All are now standardized linear
  functions of the site frequency spectrum with exact neutral variances from
  Fu's (1995) covariance (Achaz 2009), which reproduce the published closed
  forms (Tajima 1989; Fu and Li 1993 for D*; Zeng et al. 2006) and are checked
  against coalescent simulation for F*. Only columns scored in every sequence
  are used. The H p-value now comes from Zeng et al.'s normalized H and needs
  an outgroup (`calculate_fay_wu_h(seqs, outgroup)`,
  `neutrality_test_suite(seqs, outgroup=...)`, new `normalized_fay_wu_h`);
  without one it is `None`.
- `calculate_summary_statistics(genotype_matrix=...)` decoded 0/1/2 dosages as
  allele pairs `(g // 2, g % 2)`, so homozygous-alternate genotypes counted as
  heterozygous. Only dosage 1 is heterozygous now; other values are missing.
- Effective-population-size functions returned placeholders:
  `exponential_growth_effective_size` returned 0.8 or 0.9 times the initial
  size, `two_epoch_effective_size` ignored the time since the change, and
  `demography.bottleneck_effective_size` used invented weights. They now
  return Wright's harmonic mean over the growth trajectory or bottleneck
  generations, and the coalescent effective size
  N_recent + (N_ancient - N_recent) exp(-T / 2N_recent) for two epochs.
  `statistics.bottleneck_effective_size` omitted the generation count from
  the harmonic mean.
- Benjamini-Hochberg q-values in `networks.analysis.pathway` enforced
  monotonicity from the wrong end, so one enriched pathway made every null
  pathway significant.
- Plot functions given axes inside a matplotlib SubFigure crashed when saving,
  and four phenotype plots returned the SubFigure instead of its figure; they
  now resolve the root figure (`visualization.config.figures.root_figure`).
- Pure-Python fallbacks in `structural_variants` used the upper median for
  merged and called breakpoints and insert-size statistics, returned a zero
  insert-size spread for a single observation, and skipped the median
  centering and GC correction of CNV log2 ratios.
- GWAS: the *P. barbatus* pipelines called the QC filter with a stale
  signature and mismatched threshold names, then analyzed the unfiltered
  variants. They read significance from a key the corrections do not return,
  so no variant was ever flagged, left Bonferroni-adjusted p-values above 1,
  and matched phenotypes to samples by position. Phenotypes now align by
  sample ID, and missing, extra or duplicate IDs are refused. `run_gwas` coerces quoted YAML
  QC thresholds to numbers. LDSR reported the slope's standard error times
  the LD-score standard deviation as the intercept's standard error; it now
  uses the OLS intercept standard error. `plot_qq_overlay` crashed on NumPy
  p-value arrays.
- Spatial: `gearys_c` raised for COO, DIA and BSR sparse weights, and
  `correlation_mapping` used Spearman for any method other than Pearson; it
  now refuses unknown methods. The 3-D protein structure plot passed
  `projection` to `plt.subplots`, which raises.
- `scripts/cloud/sync_cloud_down.sh` disabled SSH host-key checking; it now
  uses `accept-new`, and the gcloud location comes from `GOOGLE_CLOUD_SDK_BIN`
  instead of a hard-coded home directory.
- `scripts/menu/run_menu.py` imported removed modules and could not start.
- Sixty-two scripts under `scripts/` failed at import because the package
  had moved or stopped re-exporting the names they used; they now run against
  the current APIs, checked with `--help` and small real runs. Fixing them
  exposed scripts that never flagged significant association results, swapped
  the Lotka-Volterra predator parameters, and selected regression features
  with a classification F-test. `download_honeybee_variants.py` prints its
  download guide instead of fetching whole SRA projects, and options backed by
  removed functions fail with an explicit message. A slow test resolves every
  `metainformant` import under `scripts/`; `discover_species.py` and
  `setup_genome.py` remain listed as needing functionality the package no
  longer provides.
- `simulation.models.popgen.generate_two_populations(fst=...)` produced
  F_ST near 0.5 whatever the target. It now follows the Balding-Nichols model,
  whose Hudson F_ST equals the target in expectation, and by default draws
  Watterson's expected number of polymorphic sites from `theta`.
- `visualization.config.conventions.save_figure_deterministic` could not
  write SVG, JPEG, TIFF, WebP or PGF because it passed PNG metadata to every
  backend; SVG output is now byte-reproducible. `circular_tree_plot` raised
  for trees with more than one level of internal nodes, and the matplotlib
  fallback of `interactive_scatter` failed on categorical colour columns.
- `networks.analysis.pathway.pathway_enrichment_analysis(method="hypergeometric")`
  divided by zero when the query covered a whole pathway, did not restrict
  pathways to the background, and left a single result without
  `corrected_p_value`.
- Statistics reported as medians took the upper middle value for even counts
  in structural-variant consensus calls, long-read summaries, multi-omics
  pathway concordance, network summaries, pathway activity scores, FASTQ
  quality, AntWiki phenotype distributions, life-event group comparisons and
  LDSR's genomic-control lambda. Network summaries without NumPy also
  reported a zero degree standard deviation.
- GWAS association: the logistic test used missing genotypes (-1) as
  dosages and always reported `converged: True`; it now fits complete cases,
  reports `n_missing_dropped`, and reports statsmodels' actual convergence.
  The linear test computed single-predictor p-values from the normal
  distribution instead of Student's t, and reported singular designs,
  constant genotypes and too few samples as successful fits with beta 0 and
  p 1; these are now errors. QC never applied `min_qual`, and the filtered
  VCF writer printed `None` for missing fields. `filter_variants_bcftools`
  passed its keep-conditions to `bcftools filter -e`, which removed the
  passing sites; it now uses `-i`.
- GWAS fine-mapping paired a beta-scale prior variance with V = 1 on the
  Z scale, which flattened every approximate Bayes factor: a z = 7.5 signal
  among 299 nulls gave a 285-variant 95% credible set, and `eqtl_coloc`
  returned PP.H4 near 0.003 whether or not the traits shared a causal
  variant. Wakefield's ABF now uses V = SE^2 (or sample size with MAF and
  case fraction, as in coloc) with coloc's default priors, in log space so
  strong signals stay finite; Z-only input uses an explicit `z_scale` model.
  The credible set behind the interactive fine-mapping plot used ABF = 1/p
  and now delegates to the same implementation, and Z-to-p conversion no
  longer returns 0 beyond |z| = 8.3.
- `write_summary_statistics` wrote `Q_FDR=1` and `MAF=0` when results lacked
  them; Q_FDR is now Benjamini-Hochberg over the results, MAF comes from
  optional genotypes, and any other missing value is `NA`.
  `interactive_manhattan` drew every variant on chromosome 1, and GWAS plots
  failed when the output directory did not exist.
- `multiomics.analysis.integration.canonical_correlation` ignored
  `regularization`; it is now closed-form ridge CCA. Single-cell
  `compute_pca(use_hvgs=True)` used every gene. `cross_validate_biological`
  turned metric failures into null scores; failures now raise and
  `compare_classifiers` never picks a failed method.
- `life_events` methods were placeholders. `learn_event_embeddings` returned
  random vectors nudged by event-name keywords and ignored its window, epoch
  and algorithm settings; it now trains skip-gram or CBOW word2vec with
  negative sampling. `shapley_values` returned random numbers; it now
  computes exact Shapley values (seeded permutation sampling above ten
  events) that sum to the change in prediction. The simple predictor labelled
  nearly every sequence as class 1 and lost its model on reload, LSTM and GRU
  models reported themselves fitted without torch, unknown outcome
  relationships silently became random, and JSON loading dropped sequence
  metadata and event fields. `MultiTaskPredictor` predicted each task's
  training mean for every sequence; it now fits ridge and logistic heads on
  event counts and learned embeddings. `SurvivalPredictor` reported itself
  fitted without torch and predicted a constant 365 days; it now fits a Cox
  model by Newton-Raphson on the Breslow partial likelihood (matching
  statsmodels), an exponential model or Kaplan-Meier, and predicts restricted
  mean survival time. `generate_synthetic_life_events` silently ignored six
  options; it applies the rare-event, missing-data and temporal-noise options
  and rejects the rest. `survival.cox_ph_model` used a diagonal Newton step
  and sequential ties, giving standard errors about a third of the true ones
  with correlated covariates; it now shares the Breslow Cox fit with full
  inverse-information standard errors.
- The AntWiki scraper never consulted robots.txt despite its `check_robots`
  setting; it now enforces robots rules and Crawl-delay on every request and
  redirect. Its species index was a hardcoded list of five species and now
  raises `NotImplementedError`, and output filenames derived from remote page
  titles could escape the output directory. OBO `def:` parsing left the
  closing quote and dbxrefs in definitions; they are parsed and written back
  correctly, with xrefs kept in `metadata["definition_xrefs"]`.
- `phenotype.run_phewas` computed p-values for 30 or fewer degrees of
  freedom from a homemade series that returned 0.0065 where Student's t
  gives 0.46; it now uses SciPy's t distribution. Simpson's evenness divided
  by the Gini-Simpson index, so a perfectly even community scored 1/(S-1)
  instead of 1.
- Sequence repeat and palindrome scans in
  `information.analyze_sequence_information` took super-linear time (115 s
  for a 500 bp tandem repeat); they now return the same results in under a
  second. `epigenome.load_methylation_cov` misread Bismark coverage columns,
  `dna.phylogeny.bootstrap_support` mutated its input tree, `parse_vcf` could
  not read gzipped VCFs, spatially variable gene detection rejected sparse
  matrices, HTML output without Plotly failed after creating files, and
  config workflows marked unexecuted processing steps and failed downloads
  as completed.
- `core.io.paths.ensure_directory` returned `None` while the same-named
  `core.io.ensure_directory` returns the directory; both now return it.
  Simulation seeding passed `randint(0, 2**32)`'s inclusive upper bound to
  `np.random.seed`, which rejects it; that one draw now maps to 0, leaving
  every other seeded stream unchanged.
: a blank accession
  requested InterPro's unbounded listing of every UniProt-annotated entry,
  which outlasted the 30 s read timeout and failed the network test lane on
  every main push. Blank input now returns `[]` without a request, and
  accessions are URL-quoted into a single path segment.
- `core.utils.symbols`: every scan root inside a checkout shared one index
  cache file, so indexing one subtree could return another root's cached
  index. Cache files are now keyed by scan root.
- `security_checks.sh`: a scanner that produced no report (for example
  pip-audit failing to provision its environment) was treated as skipped and
  the blocking gate passed, a stale report from an earlier run could stand in
  for a failed scan, and bandit's own read or parse errors never blocked
  because bandit still exits 0. Each scanner now clears its report first, and
  the unit-tested `security_evaluate.py` evaluates only the scanners that ran:
  a missing, unreadable, empty or partial report, any bandit scan error or
  zero-line scan, and a pip-audit report that misses a locked package all
  block. pip-audit baseline locations now match package names exactly.
- RNA reference alias manifest is content-deterministic: removed the
  `generated_at` wall-clock field from the hashed payload so an unchanged
  reference produces a byte-identical manifest across orchestrator restarts.
  Previously every restart changed the manifest digest and quarantined the
  species' entire prior quantification for redundant re-download and
  re-quantification (observed 2026-08-30: 1,620 apis_mellifera samples
  re-queued by a single restart). Added a byte-idempotency regression test.
- Scientific-correctness fixes: Miller-Madow entropy bias correction now adds
  `(k-1)/(2n*ln2)` bits (was subtracting a nats-scale term, overcorrecting
  sign and scale in `_plugin_entropy_estimator`, `_miller_madow_entropy_estimator`,
  `bias_correction`, and `effective_sample_size_correction`); Panzeri-Treves
  uniform fallback produces counts consistently; simulation
  `generate_genotype_matrix` samples one allele frequency per SNP (not per
  individual); GWAS logistic regression raises instead of fabricating
  standard errors when statsmodels is unavailable; popgen
  `genotype_structure_analysis` transposes to sites-by-individuals before
  HWE testing; `joint_pca` layer weights now apply after standardization;
  multiomics similarity imputes per-column means; `joint_nmf` defaults to a
  fixed random seed.
- Shared sequence primitives (`reverse_complement`, `validate_dna_sequence`,
  complement tables) moved to `metainformant.core.sequence` so non-DNA
  domains use them without violating the cross-domain boundary gate; seven
  duplicated local implementations consolidated; per-call `str.maketrans`
  table construction hoisted to module constants.
- Dead code removed (~2.6k lines): unreferenced `visualization/analysis`
  plot modules (`math_plots`, `simulation_plots`), unreferenced GWAS
  plotting/config/SRA helpers, dead private helpers, and duplicated
  `_parse_function_signature` / `sha256_file` copies.
- Tests: the conftest auto-`network` marker heuristic matches whole words
  (22 pure-CPU tests no longer misfiled into `-m network`); unseeded
  `np.random` visualization tests seeded; wall-clock ratio assertions in
  core parallel/TUI tests replaced with deterministic correctness checks;
  tautological InterPro/Entrez `assert True` catch-alls replaced with
  real response assertions; bare "should not raise" validation tests now
  assert returned values.
- CI/infra: pre-commit hook versions aligned with pyproject pins
  (black 26.5.1, isort 8.0.1, flake8 7.3.0); `pytest-xdist` detection fixed
  (`import xdist`); uv cache keyed on `uv.lock`; dead tag-push trigger
  removed from `release.yml` (Release events are the publish path);
  `uv_quality.sh` aligned with the blocking CI quality gate; Dockerfile
  pins `amalgkit=0.16.60`, drops standalone kallisto/fastp layers shadowed
  by the conda env, and installs the full CRAN package set the dependency
  verifier expects.
- Docs: CONTRIBUTING install extras fixed (`.[dev]`; `test`/`docs` never
  existed), duplicate CHANGELOG Unreleased sections merged, stale MCP
  module description updated to the implemented stdio server, measured
  README/AGENTS and plotting-module counts, Mermaid nav-graph node ids
  repaired in `docs/index.md`.

### Security
- dnspython 2.9.0 (CVE-2026-114605) and pymongo 4.19.0 (CVE-2026-106435),
  with matching security floors in `[tool.uv] constraint-dependencies`.
- Resolved all 60 open Dependabot advisories in `uv.lock`: pillow 13.0.0,
  certifi 2026.10.23, filelock 4.2.4, nbconvert 7.17.2, setuptools 84.2.0 and
  virtualenv 21.14.13 (also clearing GHSA-cx6f-vjhg-g3cc, a prompt-escaping
  regression introduced in 21.14.9), plus the dev pin black 26.10.1, mirrored
  in pre-commit. `[tool.uv] constraint-dependencies` holds security floors for
  these transitive packages so a re-lock cannot regress into an advisory
  range; the floors bound resolution only and do not change package metadata.

## [0.4.0] - 2026-08-09

### Breaking

- Removed deprecated compatibility import paths, legacy test-wrapper entrypoints,
  and positional cross-species alignment behavior. See
  `docs/MIGRATION_0.4.md` for the supported replacements.

### Added

- Explicit NCBI contact/anonymous-mode policy and campaign-local SRA subprocess
  environments.
- Provenance-qualified RNA readiness, cleanup authorization, and release gates.

## [0.2.7] - 2026-03-01

### Added

- **`scripts/rna/fetch_pending_metadata.py`** — batch metadata fetcher for SRA runs not yet in the pipeline manifest; handles rate limiting and incremental writes
- **`scripts/rna/find_missing_tissues.py`** — scans all species configs and tissue patches to report bioprojects with incomplete tissue coverage; used for targeted tissue-patch research
- **`scripts/rna/extract_bioproject_details.py`** — fetches NCBI BioProject title and SRA run counts for a list of accessions; used alongside `find_missing_tissues.py`
- **`scripts/rna/report_completed.py` / `report_completed.sh`** — reports quantified sample counts per species from disk, distinct from log-based counts; useful for progress monitoring independent of the orchestrator log
- **`scripts/rna/link_indices.sh`** — symlinks shared kallisto indices across species output dirs to avoid redundant index builds
- **`scripts/gwas/run_pbarbatus_gwas.py`** — end-to-end GWAS runner for *Pogonomyrmex barbatus*: downloads SRA reads, aligns with BWA, calls variants with BCFtools, and runs association with real phenotype data

### Changed

- **`.gitignore`** — added `*.heartbeat.json`, `**/.downloads/`, and `*.safely_removed` to prevent pipeline sentinel files from polluting git history
  - `*.heartbeat.json`: transient JSON progress files written during downloads; write-only, never read back, auto-recreated by active jobs
  - `*.safely_removed`: marker files confirming FASTQ deletion after successful kallisto quantification; redundant given that a non-empty `*_abundance.tsv` is the canonical proof of work
- **Pipeline housekeeping**: removed 3,898 existing heartbeat files and 3,055 `.safely_removed` markers from disk (∼7k stale sentinel files cleaned up)

### Pipeline Status (as of this release)

- *Apis mellifera*: **2,661 samples quantified** (of 5,715 total; ~3,900 auto-skipped for size >4B bases, pending second pass)
- *Acromyrmex echinatior*: **3 samples quantified** (pipeline just started)
- All other 21 species: queued — pipeline processes species sequentially

## [0.2.6] - 2026-02-20


### Fixed

- t-SNE perplexity auto-adjustment: threshold corrected to `(n_samples-1)/3` per sklearn
- `compute_diffusion_map`: removed `sigma=1.0` from `eigsh` that caused singular matrix crash; added diagonal regularization and dense fallback
- `compute_neighbors`: `n_pcs` now stored in neighbor params for downstream inspection
- Disk space: cleaned `.uv-cache` (476M) and ran `git gc --prune=all` (315MB recovered)

## [0.2.5] - 2026-02-20

### Changed

- Version bump to 0.2.5 reflecting accumulated test hardening
- CHANGELOG updated with v0.2.3–v0.2.4 entries
- ENA mention added to `scripts/rna/AGENTS.md`

## [0.2.4] - 2026-02-20

### Fixed

- `SingleCellData` diffusion map: `len(data)` → `data.n_obs` (TypeError fix)
- `ENADownloader` re-exported from `rna.retrieval.__init__`
- RNA test config YAML paths corrected (`blue/` → `output/`)
- `run_workflow.py` startup: deferred heavy imports for fast `--help` (<1s)
- Early-exit for no-args/`--list-configs` before heavy imports

## [0.2.3] - 2026-02-20

### Fixed

- **Protein**: 20 relative → absolute imports in `orchestration.py`; `read_taxon_ids` return type `List[str]`
- **Visualization**: 7 plot function re-exports; config self-reference; 22 analysis re-exports; phylo import fix; `WONG` palette constant
- **Quality**: 6 standalone contamination functions rewritten with per-sequence-index returns; `per_base_quality` quantiles guard for <2 data points; `generate_contamination_report` METAINFORMANT header
- **RNA**: `calculate_md5`, `clean_stagnant_file`, `verify_gzip_integrity` implemented in `ena_downloader.py`
- **Phenotype**: `extract_phenotypes_from_events`, `aggregate_temporal_phenotypes`, `map_events_to_traits` rewritten for `life_events.core.events` compatibility; `Measurement`/`MorphometricProfile` re-exported from `morphological/__init__`; backward-compat `career_changes`/`health_score` keys
- **Singlecell**: PCA `n_components` auto-clamp to `min(n_samples-1, n_features)`

## [0.2.2] - 2026-02

### Added

- Comprehensive repo-wide documentation audit and consistency fixes
- Specialized modules: longread, metagenomics, pharmacogenomics, spatial, structural_variants
- Enhanced Mermaid architecture diagrams with correct node references

### Changed

- Updated all documentation to use consistent plot count (70+)
- Synchronized version references across README, CHANGELOG, PAI

### Fixed

- Broken Mermaid node IDs in README data flow diagrams
- Stale version references (0.2.0 → 0.2.2) across all docs

## [0.2.1] - 2026-02

### Added

- Streaming RNA-seq orchestrator with size-ordered scheduling
- Tiered fallback recovery for failed pipeline steps
- Concurrency-hardened workspace initialization (mkdir race conditions)
- Clean-restart operational protocols for large-scale transcriptomic runs

### Changed

- Pipeline optimized for 16 parallel workers with multi-threaded scaling
- Improved subprocess timeout handling (30-minute cap on downstream calls)

### Fixed

- Pipeline stalls from missing subprocess timeouts
- FASTQ extraction directory structure handling
- Permission errors on external volume log writes

## [0.2.0] - 2025-12

### Added

- Complete UV package management integration (replaces pip)
- Enhanced RNA workflow with automated genome indexing
- Comprehensive GWAS visualization suite (manhattan, QQ, regional plots)
- Information theory module with syntactic and semantic measures
- Life events analysis with sequence modeling
- Single-cell analysis with Leiden clustering
- Multi-omics integration workflows
- Progress tracking with real-time TUI
- Intelligent caching for expensive computations

### Changed

- Python 3.11+ minimum requirement (previously 3.10+)
- Reorganized configuration system with environment variable overrides
- Improved error handling with context-aware messages
- Enhanced parallel processing with better thread management

### Fixed

- Import errors reduced from ~225 to ~63 (72% improvement)
- Test collection success rate improved to 87%
- Memory optimization for large datasets

### Removed

- pip-based installation (use `uv` exclusively)
- Legacy configuration formats

## [0.1.0] - 2025-06

### Added

- Initial release with core bioinformatics modules
- DNA analysis: sequences, alignment, population genetics, phylogeny
- RNA analysis: amalgkit integration, workflow orchestration
- Protein analysis: sequences, structures, AlphaFold integration
- GWAS pipeline: association testing, QC, visualization
- Network analysis: community detection, centrality measures
- Visualization: 57+ plot types
- Quality control: FASTQ analysis, contamination detection
- Machine learning: classification, regression, feature selection
- Mathematical biology: coalescent, population genetics theory

### Notes

- Real-implementation policy established: all tests use real implementations
- UV-only package management policy established
- Output directory isolation policy established
