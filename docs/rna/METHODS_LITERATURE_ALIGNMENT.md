# Methods Literature Alignment — Round 3 (2026-09-01)

> Lane R3_SYN of the 2026-09-01 MetaInformAnt / hymenoptera_amalgkit improvement
> round. This document aligns four published insect/bilaterian transcriptome
> methods papers with the current MetaInformAnt pipeline and the
> hymenoptera_amalgkit child manuscript, records what is adopted, planned, or
> rejected, and states the claim boundaries that apply to each adopted method.
>
> **Status:** design/alignment document. It contains no results from the live
> data root. The evidence manifest is not frozen; all cross-species analysis
> remains descriptive-only (see [claim boundaries](#claim-boundaries-and-freeze-gates)).

## Sources (verified 2026-09-01 via Crossref + Europe PMC)

1. **[R3-P1]** Yuan H, Jia H, Zhou X, Li H, Wu K. (2026). Comprehensive
   stage- and tissue-specific transcriptome of the global ecosystem service
   insect, marmalade hoverfly *Episyrphus balteatus*. *Scientific Data* 13:789.
   DOI: [10.1038/s41597-026-07148-9](https://doi.org/10.1038/s41597-026-07148-9).
   PMC13216633.
2. **[R3-P2]** Xu H, Colgan TJ. (2025). Localized tissue-specific gene
   expression and gene duplications are important sources of social morph
   differences in a social bumblebee. *Molecular Biology and Evolution*
   42(6):msaf063. DOI:
   [10.1093/molbev/msaf063](https://doi.org/10.1093/molbev/msaf063).
   PMC11968646. Study species *Bombus terrestris* is in the configured
   27-species roster.
3. **[R3-P3]** Leader DP, Naseem MT, Halberg KV. (2024). BeetleAtlas: an
   ontogenetic and tissue-specific transcriptomic atlas of the red flour beetle
   *Tribolium castaneum*. *Journal of Molecular Biology* 436(10):168520. DOI:
   [10.1016/j.jmb.2024.168520](https://doi.org/10.1016/j.jmb.2024.168520).
4. **[R3-P4]** Mantica F, Iñiguez LP, Marquez Y, et al. (2024). Evolution of
   tissue-specific expression of ancestral genes across vertebrates and
   insects. *Nature Ecology & Evolution* 8:1140–1153. DOI:
   [10.1038/s41559-024-02398-5](https://doi.org/10.1038/s41559-024-02398-5).
   Pipeline code: [github.com/fedemantica/bilaterian_GE](https://github.com/fedemantica/bilaterian_GE)
   (readme inspected 2026-09-01).

The transferable core shared by these papers is: **(a)** tau tissue-specificity
index; **(b)** orthology-class stratification (one-to-one vs multicopy /
paralog retention); **(c)** duplication–specialization coupling summaries;
**(d)** phylogeny-aware tissue-specificity gain/loss accounting; and
**(e)** atlas-style stage × tissue × species visualization.

## Paper-by-paper alignment

### R3-P1 — Yuan+ 2026, *E. balteatus* stage/tissue transcriptome

**Methods inventory (from the paper's abstract/data resource):** 30 RNA-seq
libraries spanning 3 life stages (egg, larva, pupa) and 7 adult tissues;
de novo transcriptome assembly (85,676 unigenes, N50 = 1,028 bp, 97.9% BUSCO);
functional annotation of 53.08% of unigenes; DESeq2 differential expression
across stage and tissue contrasts; public data release as a community resource.

**What our pipeline already covers:**

- Acquisition → quantification → merge → filter → finalize → sanity for public
  SRA/ENA RNA-seq, per species:
  `projects/hymenoptera_amalgkit/doc/02_workflow/`,
  `projects/hymenoptera_amalgkit/docs/manuscript/methods.md` (nine-stage
  contract), `projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`.
- Kallisto transcript-level quantification with per-run provenance sidecars:
  `src/metainformant/rna/analysis/validation.py`,
  `projects/hymenoptera_amalgkit/scripts/record_downstream_provenance.py`.
- Count normalization including TPM:
  `src/metainformant/rna/analysis/expression_core.py::normalize_counts`
  (`_normalize_tpm`).
- Exploratory differential expression machinery (DESeq2-like, t-test,
  Wilcoxon) with p-value adjustment:
  `src/metainformant/rna/analysis/expression_analysis.py::differential_expression`.
- Tissue/stage-resolved study design precedent (the paper's 3 stages × 7
  tissues) mirrors our metadata harmonization targets
  (`docs/manuscript/statistical_analysis_plan.md` §3).

**What we adopt (mapped to Round-3 lanes):**

- The stage/tissue-resolved release structure — data dictionary, figure plan,
  per-tissue summaries — informs the atlas visualization grammar [R3-ATLAS].
- Its demonstration that a multi-tissue public-data resource is publishable as
  a *descriptive* resource (assembly metrics + DE gene sets, no cross-species
  inferential claim) matches our two-stage release design: descriptive first,
  conditional biological layer only after freeze.

**What we explicitly reject (and why):**

- De novo transcriptome assembly as a per-species step. Our pipeline
  deliberately quantifies against configured, validated reference transcriptomes
  with SHA-256-bound indices
  (`projects/hymenoptera_amalgkit/docs/manuscript/methods.md`,
  "Reference-target integrity"). Adding per-species de novo assembly would
  create a second, non-comparable feature space across the 27 species and
  break the transcript-identifier lineage that
  `projects/hymenoptera_amalgkit/scripts/build_ortholog_bridge.py` depends on.

**Claim-boundary notes:** the paper's resource is single-species; nothing in it
authorizes cross-species statements. Its DE analyses are within-species and
contrast-based — admissible here only post-freeze under the harmonized-contrast
rules of `statistical_analysis_plan.md` §5.1.

### R3-P2 — Xu & Colgan 2025, *B. terrestris* tissue-specificity and duplications

**Methods inventory (from the paper and the round brief):** TPM
normalization; removal of the lowest-10%-mean-expression genes; log2 transform
before tau computation (via `tispec`); tissue-specificity index tau per Yanai
et al. 2005; orthology classification of genes into single-copy orthologs vs
multicopy/paralog groups; comparison of morph-biased vs unbiased genes for tau
and duplication enrichment; dN/dS comparison for biased protein-coding genes;
Wilcoxon/chi-square style tests for distribution differences.

**What our pipeline already covers:**

- TPM normalization (`expression_core.py::normalize_counts`, method="tpm").
- Low-expression filtering
  (`expression_core.py::filter_low_expression`) — thresholds differ from the
  paper's lowest-10%-mean rule, which we adopt explicitly (see below).
- Orthogroup bridge construction OrthoDB → NCBI Protein → NCBI RNA → Kallisto
  transcripts:
  `projects/hymenoptera_amalgkit/scripts/build_ortholog_bridge.py`;
  transcript–orthogroup tables in
  `src/metainformant/rna/analysis/ortholog_mapping.py::build_transcript_orthogroup_table`.
- Ortholog expression mapping and conservation:
  `src/metainformant/rna/analysis/cross_species.py::map_expression_to_orthologs`,
  `cross_species.py::compute_expression_conservation`.
- Differential machinery (Wilcoxon) exists in
  `expression_analysis.py::differential_expression` for gated reuse.

**What we adopt (mapped to Round-3 lanes):**

- **Tau tissue-specificity, exactly per the paper's protocol:** TPM input,
  drop the lowest-10%-mean-expression genes, log2 transform, then
  tau = sum(1 − x_i / max(x_i)) / (n − 1) over tissues with max > 0.
  Zero-expression genes are excluded; single-tissue matrices are undefined
  (not 0). → [R3-CORE]
  `src/metainformant/rna/analysis/tissue_specificity.py` (planned path,
  Round-3 core lane).
- **Orthology-class stratification:** classify genes as single-copy ortholog
  vs multicopy within each orthogroup via the bridge schema, and stratify tau
  distributions by class. → [R3-CORE] join helpers in the same module.
- **Duplication–specialization coupling (descriptive):** medians/quantiles of
  tau by orthology class in the default output. → [R3-CORE] and [R3-COMP]
  coupling tables (`tissue_specificity_evolution.py`, planned path).
- **Gated inferential reuse:** Wilcoxon rank-sum comparisons of tau
  distributions across classes (as in the paper) are implemented as a separate,
  explicitly-gated function with post-freeze labeling — never in default
  descriptive output. → [R3-CORE] gated API (planned path).

**What we explicitly reject (and why):**

- dN/dS estimation for biased-gene sets. Our scope this round is expression;
  coding-sequence evolution would require codon alignments and a validated
  ortholog-to-CDS mapping we do not currently build
  (`build_ortholog_bridge.py` produces protein/RNA links, not CDS
  alignments). Rejected as out of scope for the freeze; may be proposed later
  as its own gated extension.
- Morph-biased (queen vs worker) DE contrasts at this stage: the paper's core
  contrast requires a harmonized queen/worker design that the current cohort
  manifest does not yet guarantee (`statistical_analysis_plan.md` §3, §8
  stopping rules). The *tissue-specificity × orthology* half of the paper is
  adopted; the *morph-bias* half remains a post-freeze plan item.

**Claim-boundary notes:** tau is a descriptive index of expression
concentration across the tissues sampled in a given matrix. It is
sampling-dependent: a species with few available tissues yields tau values
that are not comparable to tau from tissue-rich species. Any cross-species tau
table must therefore report the per-species tissue count and mask
uncomparable pairs (adopted in the [R3-ATLAS] heatmap conventions).
*B. terrestris* itself is in the roster, so this paper's protocol is directly
transferable to at least one configured species.

### R3-P3 — Leader+ 2024, BeetleAtlas

**Methods inventory (from the paper's abstract):** quantitative expression
database for 9 adult + 7 larval tissues and 4 embryonic stages of
*T. castaneum*; web application exposing per-gene total expression and tissue
enrichment, isoform-level data; cross-species query via *Drosophila* gene
identifiers; five search modes (tissue-high expression, larva–adult
difference, stage peaks, functional category, profile similarity).

**What our pipeline already covers:**

- Expression matrices, feature statistics, and quality summaries per species
  (`run_cross_species_analysis.py` outputs: `feature_statistics.tsv`,
  `profile_quality.tsv`).
- Ortholog-linked gene lookup across species via the bridge (above).
- Descriptive figure set with captions and accessibility records
  (`figure_captions.md`, `figure_accessibility.md`).

**What we adopt (mapped to Round-3 lanes):**

- **Atlas grammar** — organized expression heatmaps with tissue/stage
  hierarchy, ranked tissue-specificity summaries, and cross-species ortholog
  profile displays. → [R3-ATLAS]
  `src/metainformant/rna/analysis/atlas_plots` (planned path, Round-3 atlas
  lane): (a) species × tissue tau heatmap with deterministic fixed ordering;
  (b) per-orthogroup cross-species expression-profile small multiples;
  (c) tau-vs-orthology-class summary strip plots (descriptive quantiles).
- **Query-mode thinking** (tissue-high genes, profile-similar genes) as a
  candidate post-freeze feature for the child project's exploratory interface;
  recorded as a plan item, not built this round.

**What we explicitly reject (and why):**

- A web application / interactive database deliverable. This round's contract
  is a deterministic, evidence-manifest-bound figure module
  (matplotlib, MPLBACKEND=Agg, byte-stable outputs); an interactive resource is
  a separate infrastructure decision with its own hosting/validation gates.
- Isoform-level atlas displays: our finalized matrices are the current
  finalize-stage outputs; isoform resolution would require re-running
  quantification-side tooling and is not part of the frozen-cohort contract.

**Claim-boundary notes:** BeetleAtlas itself performs within-species
enrichment/contrast queries; its cross-species reach is identifier translation,
not comparative inference. We keep that boundary: atlas figures are
descriptive displays of the current snapshot cohort and support no
significance or convergence claim.

### R3-P4 — Mantica+ 2024, tissue-specificity evolution across bilaterians

**Methods inventory (from the paper's abstract and the bilaterian_GE readme):**
8 tissues × 20 bilaterian species; orthology by Broccoli preliminary runs plus
curated best-ancestral ortholog selection using sequence similarity and
expression correlation within orthogroups; tau for all protein-coding genes in
each species; a **symmetric phylogeny** treatment (duplications create
vertebrate/insect sub-trees rather than a single species tree); systematic
inference of tissue-specificity gains and losses per tissue across the
phylogeny; association of gains with gene duplication and expression
specialization of one copy; robustness checks (alternative ortholog
quantification, sample randomization across meta-samples, orthogroup
randomization).

**What our pipeline already covers:**

- Orthogroup-level tables joining expression to orthology
  (`ortholog_mapping.py::build_transcript_orthogroup_table`,
  `cross_species.py::map_expression_to_orthologs`).
- Descriptive expression conservation per ortholog
  (`cross_species.py::compute_expression_conservation`).
- Species-tree provenance requirements already encoded in the analysis plan
  (`statistical_analysis_plan.md` §5.3 requires a versioned species tree for
  phylogenetic association).

**What we adopt (mapped to Round-3 lanes):**

- **Tissue-specificity state assignment:** a gene/orthogroup is
  "tissue-specific" in a species when tau exceeds the threshold convention
  adopted from the paper-2 protocol (recorded per run; not a universal
  constant). → [R3-COMP]
  `src/metainformant/rna/analysis/tissue_specificity_evolution.py` (planned
  path).
- **Gain/loss counting as descriptive parsimony counting** on a supplied
  Newick tree: per-node counts of tissue-specificity gains and losses per
  orthogroup, reported as counts only — no rate model, no ancestral
  probabilistic inference, no significance. → [R3-COMP] (planned path).
- **Duplication–specialization coupling table:** orthogroups whose gain events
  coincide with multicopy class membership, summarized descriptively. →
  [R3-COMP] + [R3-CORE].
- **Species-tree provenance rule:** the 27-species tree source (TimeTree /
  OrthoDB / hymenoptera literature) must be recorded with version and
  branch-length scale in the design doc — [R3-COMP]
  `doc/05_cross_species/TISSUE_SPECIFICITY_EVOLUTION_DESIGN.md` in the child
  repo (planned path; link below).

**What we explicitly reject (and why):**

- The paper's OU-comparison / model-based validation of gains (bilaterian_GE
  readme, Extended Data Fig. 7 machinery): model-fitting is inferential and is
  barred pre-freeze.
- "Best-ancestral ortholog" re-selection via expression-correlation filtering
  of orthogroups: a powerful but assumption-laden curation step; adopting it
  would redefine our orthology space mid-campaign. We keep the current bridge
  policy (one-to-one preferred, one-to-many retained only under the declared
  aggregation rule of `statistical_analysis_plan.md` §4) and record the
  paper's approach as an alternative for a future gated sensitivity lane.
- Whole-genome-duplication narratives: vertebrate WGD context does not
  transfer to the Hymenoptera panel; we adopt only the duplication–expression
  specialization *coupling method*, not the paleopolyploidy interpretation.

**Claim-boundary notes:** parsimony gain/loss counts are descriptive
character-counting on a supplied tree; they do not estimate rates, dates, or
selection. The symmetric-phylogeny treatment matters for us only insofar as
multicopy orthogroups may need per-copy state assignment; this is documented in
the [R3-COMP] design doc (planned path below). Any cross-species statement of
"conservation of tissue-specificity" remains barred until the evidence
manifest freezes.

## Cross-cutting adoption map

| Method element | Source paper(s) | Destination (planned paths) | Round-3 lane | Status |
|---|---|---|---|---|
| tau (TPM, drop-lowest-10%-mean, log2) | R3-P2 | `src/metainformant/rna/analysis/tissue_specificity.py` | R3-CORE | planned |
| Orthology-class (1:1 vs multicopy) stratification | R3-P2, R3-P4 | same, plus bridge schema from `build_ortholog_bridge.py` | R3-CORE | planned |
| Duplication–specialization coupling (descriptive) | R3-P2, R3-P4 | `tissue_specificity.py` + `tissue_specificity_evolution.py` | R3-CORE+R3-COMP | planned |
| Gated Wilcoxon (post-freeze only) | R3-P2 | `tissue_specificity.py` gated API | R3-CORE | planned |
| Gain/loss parsimony counting on Newick | R3-P4 | `tissue_specificity_evolution.py` | R3-COMP | planned |
| Species-tree provenance contract | R3-P4 | `doc/05_cross_species/TISSUE_SPECIFICITY_EVOLUTION_DESIGN.md` (child) | R3-COMP | planned |
| Atlas heatmap / small-multiple grammar | R3-P1, R3-P3 | `src/metainformant/rna/analysis/atlas_plots` | R3-ATLAS | planned |
| Deterministic figure contract (Agg, fixed ordering) | R3-P3 | `atlas_plots` + repo chart conventions | R3-ATLAS | planned |
| TPM / low-expression filtering (already present) | R3-P1, R3-P2 | `expression_core.py` | existing | covered |
| Orthogroup bridge (already present) | R3-P2, R3-P4 | `build_ortholog_bridge.py`, `ortholog_mapping.py` | existing | covered |
| Ortholog expression conservation (already present) | R3-P4 | `cross_species.py` | existing | covered |
| Within-species DE machinery (already present, gated) | R3-P1 | `expression_analysis.py` | existing, gated | covered |
| Profile conservation over matched tissues (Spearman) | R3-P2, R3-P4 | `conservation_profiles.py::compute_profile_conservation` | R4-RNA | implemented (2026-09-01) |
| TPM distribution summary artifacts | R3-P1, R3-P3 | `conservation_profiles.py::compute_tpm_distribution_summary` | R4-RNA | implemented (2026-09-01) |
| Per-tissue completeness + pairwise tissue-overlap tables | R3-P2 (claim boundary 2) | `conservation_profiles.py::compute_per_tissue_completeness`, `::compute_tissue_overlap_summary` | R4-RNA | implemented (2026-09-01) |

## Explicitly rejected elements (consolidated)

| Rejected element | Source | Reason |
|---|---|---|
| Inferential tests on cross-species distributions (pre-freeze) | R3-P2 | Evidence manifest not frozen; hard project boundary |
| dN/dS on biased-gene sets | R3-P2 | No validated CDS alignment lane; out of scope |
| De novo per-species assembly | R3-P1 | Breaks cross-species feature-space comparability |
| Interactive atlas web app | R3-P3 | Separate infrastructure decision; not this round |
| Isoform-level atlas displays | R3-P3 | Outside the finalized-matrix contract |
| OU-model validation of specificity gains | R3-P4 | Model-based inference barred pre-freeze |
| Best-ancestral-ortholog re-selection | R3-P4 | Would redefine orthology space mid-campaign |
| WGD narrative transfer | R3-P4 | Vertebrate-specific context |

## Claim boundaries and freeze gates

1. **Descriptive-only until freeze.** All Round-3 outputs (tau tables,
   orthology-class summaries, parsimony counts, atlas figures) are descriptive
   statistics of the current snapshot cohort. No p-values, confidence
   intervals, or significance language may appear in cross-species outputs.
   Wilcoxon/chi-square style machinery may exist in code but must be exported
   behind an explicit post-freeze gate and labeled as such.
2. **Sampling-dependence of tau.** Cross-species tau comparisons require
   per-species tissue counts and comparable tissue sets; unmatched-tissue tau
   values must be masked or reported with their denominators. This follows
   from the paper-2 protocol and is restated here as a release rule.
3. **Feature vs gene vs ortholog.** Until a validated transcript-to-gene
   mapping exists for a species matrix, outputs remain feature-level
   (`docs/manuscript/methods.md`, "Native expression-fingerprint analysis").
   Ortholog-level statements require the bridge coverage audit of
   `statistical_analysis_plan.md` §4.
4. **Tree provenance.** Gain/loss counting requires a versioned, cited species
   tree with branch-length scale and taxon-pruning rule recorded in the
   [R3-COMP] design doc; an expression dendrogram is never a substitute.
5. **No manuscript-results promotion.** None of the adopted methods may be
   cited as executed in the manuscript until their outputs, hashes, and run
   records exist and pass the evidence-manifest gate
   (`analysis_readiness.md`, "Manuscript release gate").

## Integration plan for the child manuscript

Child-side integration proceeds as a dated addendum in
`docs/manuscript/ANALYSIS_NOTES_R3_2026-09-01.md` in the child repo (see that
file), covering: (a) tau + orthology-class methods as the planned descriptive
layer; (b) the gated inferential extension path; (c) links to this alignment
document and to the [R3-COMP] design doc. The manuscript Methods section file
(`docs/manuscript/methods.md`) was reviewed on 2026-09-01: it contains no
future-facing analysis-plan subsection (that role is played by
`statistical_analysis_plan.md`), so per lane rules the manuscript structure is
unchanged and the addendum route is used.

## Cross-references

- Child addendum: `projects/hymenoptera_amalgkit/docs/manuscript/ANALYSIS_NOTES_R3_2026-09-01.md`
- Comparative design doc [R3-COMP]:
  `projects/hymenoptera_amalgkit/doc/05_cross_species/TISSUE_SPECIFICITY_EVOLUTION_DESIGN.md`
  — planned path per the Round-3 comparative lane brief; if it does not exist
  when this document is committed, it is expected at that path in a
  subsequent Round-3 commit, and the link above is the durable pointer.
- Analysis-plan boundary source: `projects/hymenoptera_amalgkit/docs/manuscript/statistical_analysis_plan.md`
- Release-gate source: `projects/hymenoptera_amalgkit/docs/manuscript/analysis_readiness.md`
