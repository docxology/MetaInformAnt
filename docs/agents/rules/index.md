# Agent Directives: Module Rules

## Role

Module-specific Agent Rules for consistent development patterns across METAINFORMANT. These files serve as the authoritative source for coding conventions, domain-specific patterns, and integration guidelines for each of the 28 bioinformatics modules.

## Contents

Each Rule file contains domain-specific guidelines:

| Rule File | Domain |
|-----------|--------|
| `core.md` | Core infrastructure (I/O, config, logging, parallel, workflow) |
| `dna.md` | DNA sequence analysis, alignment, phylogeny |
| `rna.md` | RNA-seq, amalgkit workflow orchestration |
| `gwas.md` | GWAS pipelines, association testing |
| `protein.md` | Protein sequence/structure (AlphaFold, UniProt, InterPro) |
| `epigenome.md` | Epigenomic analysis (methylation, ChIP-seq, ATAC-seq) |
| `singlecell.md` | Single-cell RNA-seq analysis |
| `spatial.md` | Spatial transcriptomics (Visium, MERFISH, Xenium) |
| `multiomics.md` | Multi-omic integration |
| `networks.md` | Biological network analysis |
| `ml.md` | Machine learning pipelines, LLM integration |
| `math.md` | Population genetics theory (coalescent, FST, LD) |
| `information.md` | Information theory (entropy, mutual info) |
| `ontology.md` | Gene Ontology, semantic similarity |
| `phenotype.md` | Trait analysis |
| `ecology.md` | Community diversity, ecological metadata |
| `simulation.md` | Synthetic data generation, agent-based models |
| `quality.md` | Quality control metrics, batch effects |
| `visualization.md` | Plotting (70+ types), publication figures |
| `longread.md` | Long-read sequencing (PacBio/Nanopore) |
| `metagenomics.md` | Microbiome analysis (amplicon, shotgun) |
| `structural_variants.md` | CNV/SV detection, annotation |
| `pharmacogenomics.md` | Clinical pharmacogenomics (CPIC, PharmGKB) |
| `metabolomics.md` | Metabolomics (MS, pathway mapping) |
| `menu.md` | Interactive CLI menu and discovery system |

**Total**: 28 module rule files.

## Usage

These rules are **automatically loaded** by AI agents (Cursor, Hermes) when working in corresponding module directories (`src/metainformant/{module}/` or `docs/{module}/`).

### Workflow

1. **Start with Core Rules** (`rules/core.md`):
   - I/O via `metainformant.core.io` (never `import json`)
   - Logging via `metainformant.core.utils.logging`
   - Package management with `uv`
   - Output to `output/{module}/`
   - Zero-mocking policy enforcement

2. **Then consult module-specific rule** for your domain:
   - Domain-specific import patterns
   - Typical workflow structure
   - Domain-specific utilities
   - Configuration variables and env prefixes
   - Cross-module dependencies

### Example

Working on `src/metainformant/rna/`? Read `rules/rna.md` first.

## Key Patterns Enforced

All rule files reinforce these universal project standards:

| Standard | How Enforced |
|----------|--------------|
| **UV Package Manager** | All examples use `uv add`, `uv run`, never `pip` |
| **No Mocks** | Tests explicitly forbid `@patch`, `Mock()` — use real implementations |
| **Core I/O** | All file ops via `metainformant.core.io` (`load_json`, `dump_json`, etc.) |
| **Structured Logging** | `from metainformant.core.utils.logging import get_logger` |
| **Type Hints** | Python 3.11+, full annotations |
| **Output Discipline** | Write only to `output/{module}/`, never `src/` |
| **Absolute Imports** | `from metainformant.dna import align` (not `..core` relative) |

## Coordination with Other Agents

### Module-to-Module Interactions

Your module may depend on other modules. Rules for safe cross-module coordination:

**Lazy imports** (optional dependencies):

```python
def optional_feature():
    try:
        from metainformant.rna import amalgkit
        return amalgkit.do_work()
    except ImportError:
        raise ImportError(
            "rna module required — install with: uv add metainformant[rna]"
        )
```

**Shared core dependencies only**: Modules should not import from other *domain* modules directly. Domain-level coordination happens through:
- **Shared config files** (YAML/JSON passed to workflow manager)
- **CLI composition** (script calls another module's CLI)
- **Pipeline phase ordering** (BasePipelineManager encodes dependency order)

**Example**: `gwas` may invoke `rna` scripts via `subprocess.run()` rather than Python import.

### Agent Coordination Patterns

For details on how multiple agents/phases collaborate:

- **[Orchestration](../ORCHESTRATION.md)** — BasePipelineManager API
- **[Workflows](../MULTI_AGENT_WORKFLOWS.md)** — Real-world multi-phase patterns
- **[Communication](../COMMUNICATION_PROTOCOLS.md)** — Metadata passing, file-based handoff

## Testing Requirements

All rule files specify testing patterns:

```markdown
## Testing
- Use `tmp_path` fixture for file isolation
- Use real implementations only — NO MOCKING (see ../../../tests/NO_MOCKING_POLICY.md)
- Test at least one happy path and one error path per public function
- Include integration test of end-to-end workflow (optional but encouraged)
```

## Best Practices (Quick Reference)

```text
 DO:
  - Use metainformant.core.io for all file I/O
  - Use resource_aware_workers() for parallel execution
  - Write to output/{module}/ with timestamped subdirs
  - Document environment variables with module prefix: {MODULE}_KEY=value
  - Update this rule file when adding new features

 DON'T:
  - Import json/csv directly — use metainformant.core.io
  - Use unittest.mock.patch for external calls
  - Write to src/ or repository root
  - Hardcode thread counts — compute resource_aware_workers()
  - Assume other modules' internal APIs are stable
```

## Cross-References

| Topic | Link |
|-------|------|
| Agent coordination hub | [../README.md](../README.md) |
| Universal agent directives | [../AGENTS.md](../AGENTS.md) |
| Orchestration API | [../../ORCHESTRATION.md](../ORCHESTRATION.md) |
| Multi-agent workflows | [../../MULTI_AGENT_WORKFLOWS.md](../MULTI_AGENT_WORKFLOWS.md) |
| Communication patterns | [../../COMMUNICATION_PROTOCOLS.md](../COMMUNICATION_PROTOCOLS.md) |
| Safety & error handling | [../SAFETY.md](../SAFETY.md) |
| Operational best practices | [../BEST_PRACTICES.md](../BEST_PRACTICES.md) |
| No-mocking test policy | [../../../tests/NO_MOCKING_POLICY.md](../../../tests/NO_MOCKING_POLICY.md) |

## Maintenance

**Updating a rule file**:
1. Edit `docs/agents/rules/{module}.md`
2. Run test suite to ensure no broken cross-links: `uv run python scripts/package/generate_cursor_skills.py --check`
3. Regenerate Cursor skills if needed: `uv run python scripts/package/generate_cursor_skills.py --check`
4. Commit with message: "docs(agents): update {module} rules for ..."

**Adding a new module**:
1. Create `docs/agents/rules/newmodule.md` using template from existing rule
2. Add entry to `rules/index.md` contents table
3. Add to `docs/index.md` module overview matrix
4. Add to `src/metainformant/` with `AGENTS.md`, `../SPEC.md`, `README.md`
5. Generate Cursor skills

---

**Parent Hub**: [Agent Coordination Hub](../README.md)  
**Specification**: [../SPEC.md](../SPEC.md)  

*These rules are part of the METAINFORMANT Core layer. Stability: HIGH. Backwards-incompatible changes require deprecation period.*
