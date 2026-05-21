# Specification: spatial

## Scope

Documentation for spatial transcriptomics workflows, including spot/cell
coordinate input, neighborhood analysis, spatial autocorrelation, niche
analysis, communication, deconvolution, and visualization.

## Boundaries

- Source package: `src/metainformant/spatial/`
- User docs: `docs/spatial/`
- Runtime outputs: `output/spatial/` or a user-specified output directory

## Interface Policy

- Prefer canonical imports from `metainformant.spatial`.
- Use shared I/O and path helpers from `metainformant.core.io`.
- Cross-domain integration should be documented as an adapter or workflow,
  not as an implicit dependency from core utilities.
