# `amalgkit csfilter`

`csfilter` is an optional cross-species filtering command. It is separate from
the required per-species `wsfilter` → `finalize` chain and must only run when
the installed CLI and the supplied cross-species input contract are verified.

## Required inputs

- finalized or explicitly declared cross-species expression inputs;
- a sample/species map with stable identifiers;
- ortholog mapping if features are not already shared;
- recorded filtering parameters and tool version.

## Evidence rule

The command's output is an input to cross-species analysis, not a significance
test. Retain the input manifest, feature intersection, excluded rows, and
species/sample counts. The Hymenoptera project uses its native analysis
script for descriptive fingerprints and does not silently label those results
as an Amalgkit cross-species correlation result.
