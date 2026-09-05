# PAI - config/life_events

## Context & Intent

`config/life_events/` configures the life-events temporal analysis module (life course and event sequence modeling). Per `AGENTS.md`: "Life events temporal analysis configuration for life course modeling." It holds one file, `life_events_template.yaml`, the documented template consumed by `scripts/life_events/run_life_events_analysis.py` and backed by `src/metainformant/life_events/`.

## Virtual Hierarchy

- Parent: `config/` (repo-wide YAML config layer; per root `SPEC.md` configs can be overridden by environment variables with domain prefixes).
- Sibling template dirs: `config/singlecell/`, `config/networks/`, `config/multiomics/` — same template conventions.
- Downstream consumers: `scripts/life_events/run_life_events_analysis.py` (thin wrapper, invoked with `--config config/life_events/life_events_template.yaml` per `README.md`) and `src/metainformant/life_events/`.

## Maintenance Notes

From `AGENTS.md` (binding rules):

- "Validate with schema before committing new configs."
- "Follow REAL IMPLEMENTATION policy — tests use real config files."
- "Use `uv` for dependency management."
- "Environment overrides use the life-events namespace documented by the module."

## AI Workflows

- **New config**: start from `life_events_template.yaml`; its sections are `embedding` (skipgram/cbow, `embedding_dim` 100, window 5), `model` (`model_type`: embedding/simple/lstm; task classification/regression), `workflow`, `output`, plus top-level `work_dir`, `log_dir`, `threads`.
- **Run**: `python3 scripts/life_events/run_life_events_analysis.py --input data/life_events/sequences.json --config config/life_events/life_events_template.yaml` (exact command from `README.md`).
- **Caveats to preserve**: `lstm` requires PyTorch and falls back if unavailable; keep `random_state` set for reproducibility.
- **Scope discipline**: edit YAML here only; never touch `scripts/life_events/` or `src/metainformant/life_events/` from config tasks.
