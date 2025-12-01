# Workflow Consolidation Status

Status document for METAINFORMANT RNA workflow script consolidation.

## Consolidation Summary

**Status**: Complete (December 2025)  
**Result**: 20+ scripts consolidated into 4-5 thin orchestrators

## Completed Actions

### Phase 1: Logic Migration to src/

All business logic moved from scripts to `src/metainformant/rna/`:

| Module | Functions Added | Purpose |
|--------|-----------------|---------|
| `orchestration.py` | `run_workflow_for_species()`, `cleanup_unquantified_samples()`, `monitor_workflows()` | Workflow orchestration |
| `cleanup.py` | `cleanup_partial_downloads()`, `fix_abundance_naming_for_species()` | Cleanup utilities |
| `discovery.py` | `search_species_with_rnaseq()`, `get_genome_info()`, `generate_config_yaml()` | Species discovery |
| `genome_prep.py` | `verify_genome_status()`, `orchestrate_genome_setup()` | Genome preparation |
| `monitoring.py` | `assess_all_species_progress()`, `initialize_progress_tracking()` | Progress monitoring |

### Phase 2: Thin Orchestrators Created

Scripts now only parse arguments and call src methods:

| Script | Purpose | src Function |
|--------|---------|--------------|
| `run_workflow.py` | End-to-end workflow | `orchestration.run_workflow_for_species()` |
| `setup_genome.py` | Genome setup | `genome_prep.orchestrate_genome_setup()` |
| `discover_species.py` | Species discovery | `discovery.*` functions |
| `check_environment.py` | Environment validation | `environment.validate_environment()` |

### Phase 3: Obsolete Scripts Removed

Deleted 20+ scripts:
- Workflow orchestration: `orchestrate_workflows.py`, `run_multi_species.py`, etc.
- Genome setup: `orchestrate_genome_setup.py`, `verify_genomes_and_indexes.py`, etc.
- Discovery: `discover_ant_species_with_rnaseq.py`, etc.
- Environment: `check_r_dependencies.py` (merged)

### Phase 4: Documentation Updated

- `scripts/rna/README.md` - Updated with new structure
- `scripts/rna/AGENTS.md` - Documents consolidation history
- `scripts/rna/amalgkit/README.md` - Updated references

## Current Script Structure

```
scripts/rna/
├── run_workflow.py          # Main orchestrator (recommended)
├── setup_genome.py          # Genome setup
├── discover_species.py      # Species discovery
├── check_environment.py     # Environment validation
├── amalgkit/
│   ├── run_amalgkit.sh      # Bash workflow runner
│   └── verify_workflow.sh   # Workflow validation
└── (utility scripts)
```

## Design Principles

1. **Thin Scripts**: Scripts only parse args and call src methods
2. **Single-Species Sequential**: Focus on one species at a time
3. **Hybrid Approach**: metainformant for complex logic, direct CLI for simple steps
4. **All Logic in src**: No business logic in scripts directory
5. **Backward Compatibility**: New orchestrators support same use cases

## Benefits Achieved

- Single source of truth for workflow logic
- Easier maintenance and testing
- Consistent behavior across all uses
- Clear separation of concerns
- Reduced code duplication

## Related Documentation

- [scripts/rna/README.md](../README.md) - Script organization
- [scripts/rna/AGENTS.md](../AGENTS.md) - Consolidation history
- [docs/rna/workflow.md](../../../docs/rna/workflow.md) - Workflow documentation

---

*Consolidation completed December 2025 by Code Assistant Agent.*

