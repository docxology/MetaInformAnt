# Integration

Cross-domain adapters for RNA orchestration. RNA engine modules may import only
`metainformant.core` and their own domain; when they need another domain, they
go through a module here, the approved boundary surface checked by
`scripts/quality/check_module_boundaries.py`.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports the `campaign_figures` adapter |
| `campaign_figures.py` | Re-exports `plot_acquisition_observation` from `metainformant.visualization.plots.rna_campaign` for the campaign observation report |

## Key Functions

| Function | Description |
|----------|-------------|
| `campaign_figures.plot_acquisition_observation()` | Render the acquisition-observation figure written into an observation bundle by `rna.engine.campaign_observation_report.write_observation_bundle()` |
