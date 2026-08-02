# Workflow execution summary template

This page is a recording template, not a completion claim. Fill it from the
current run manifest and evidence generator.

## Run identity

- Date/time (UTC):
- Git commit:
- MetaInformAnt version:
- Amalgkit version:
- Python version:
- Data root:
- Configuration hash:
- Operator/host:

## Scope

- Configured species:
- Materialized species:
- Included species:
- Excluded or blocked species and reason:
- Metadata query and selection rules:
- Tissue normalization file/hash:
- Ortholog mapping file/hash, if used:

## Current workflow evidence

| Step | Required evidence | Count/status |
|---|---|---|
| metadata | source metadata table | |
| select | selected metadata and exclusion summary | |
| getfastq | validated reads or failure records | |
| integrate | integrated metadata | |
| quant | valid abundance files | |
| merge | merged expression tables | |
| wsfilter | filtered tables and exclusions | |
| finalize | final table and sample map | |
| sanity | integrity report | |

## Cross-species lane

- Analysis type: descriptive / inferential
- Input species and sample counts:
- Feature intersection:
- Normalization:
- Model and uncertainty:
- Multiplicity correction:
- Sensitivity analyses:
- Figure/table manifests:

## Sign-off

The result is publication-ready only when the report, manuscript variables,
figures, captions, and evidence manifest all regenerate from the same inputs
and the residual limitations are stated.
