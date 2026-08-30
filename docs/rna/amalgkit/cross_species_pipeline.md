# Cross-species RNA analysis

Cross-species analysis is a separate evidence lane from per-species
quantification. It must start only after each included species has a complete
and reconciled merge/filter/finalization record.

## Optional Amalgkit branch

The installed Amalgkit environment exposes:

- `cstmm` for ortholog-aware cross-species normalization;
- `csfilter` for cross-species filtering when its input contract is met.

Example configuration:

```yaml
steps:
  cstmm:
    out_dir: /Volumes/external_drive/Data/amalgkit/cross_species
    orthogroup_table: /Volumes/external_drive/Data/amalgkit/cross_species/orthogroups.tsv
    dir_busco: /Volumes/external_drive/Data/amalgkit/cross_species/busco
  csfilter:
    out_dir: /Volumes/external_drive/Data/amalgkit/cross_species
    input_dir: /Volumes/external_drive/Data/amalgkit/cross_species/cstmm
```

Run `amalgkit cstmm --help` and `amalgkit csfilter --help` in the exact
environment used for the analysis. Save the command output, ortholog table
hash, input species list, feature intersection, and resulting row/column
counts.

## Native project analysis

The Hymenoptera project includes a native Python cross-species analysis for
descriptive expression fingerprints. It writes a manifest, pairwise summary,
matrix, and figures. These outputs are not a substitute for inferential
meta-analysis: they contain no p-values, confidence intervals, phylogenetic
model, or multiple-testing correction unless those are explicitly added and
validated.

```bash
cd projects/hymenoptera_amalgkit
uv run python scripts/prepare_cross_species_inputs.py \
  --data-root /Volumes/external_drive/Data/amalgkit
uv run python scripts/run_cross_species_analysis.py \
  --data-root /Volumes/external_drive/Data/amalgkit \
  --output-dir output/cross_species
```

Read [the project cross-species documentation](../../../projects/hymenoptera_amalgkit/doc/05_cross_species/README.md)
for inclusion rules and evidence interpretation.

## Minimum release gate

Do not release a cross-species figure until:

1. species and sample inclusion are explicit;
2. each matrix is tied to the finalization manifest;
3. ortholog mapping and feature intersection are recorded;
4. missingness and filtering are reported;
5. captions state whether results are descriptive or inferential;
6. all figures regenerate from the recorded inputs.
