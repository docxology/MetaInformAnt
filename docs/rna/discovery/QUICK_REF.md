# Quick Reference: Ant Species Discovery

## One-Line Commands

```bash
# Basic discovery (all ants with ≥1 RNA-seq sample)
python3 scripts/rna/discover_ant_species_with_rnaseq.py

# Well-studied species only (≥50 samples)
python3 scripts/rna/discover_ant_species_with_rnaseq.py --min-samples 50

# Deploy all generated configs
cp output/ant_discovery/configs/*.yaml config/amalgkit/

# Deploy specific species
cp output/ant_discovery/configs/amalgkit_camponotus_floridanus.yaml config/amalgkit/

# Run workflow for discovered species
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_SPECIES.yaml
```

## Prerequisites

```bash
# Using uv (recommended)
uv pip install biopython ncbi-datasets-pylib
export NCBI_EMAIL="your.email@example.com"
```

## Typical Discovery Results

Expected to find 30-100 ant species with RNA-seq data, including:

**High Priority** (>100 samples, chromosome genomes):
- Camponotus floridanus (Florida carpenter ant) - ~300 samples
- Solenopsis invicta (Red fire ant) - ~350 samples

**Medium Priority** (10-100 samples, good genomes):
- Pogonomyrmex barbatus (Red harvester ant) - ~80 samples
- Monomorium pharaonis (Pharaoh ant) - ~100 samples
- Acromyrmex echinatior (Leafcutter ant)
- Atta cephalotes (Leafcutter ant)

**Emerging** (1-10 samples):
- Many additional species with initial studies

## Output Files

```
output/ant_discovery/
├── DISCOVERY_REPORT.md           # Read this first!
├── ant_species_rnaseq_data.json  # Machine-readable data
└── configs/
    └── amalgkit_*.yaml            # Ready-to-use configurations
```

## Next Steps After Discovery

1. **Review report**: `cat output/ant_discovery/DISCOVERY_REPORT.md`
2. **Check priorities**: Species with chromosome genomes + many samples
3. **Deploy configs**: `cp output/ant_discovery/configs/amalgkit_*.yaml config/amalgkit/`
4. **Start workflows**: Use existing scripts with new configs
5. **Monitor progress**: `python3 scripts/rna/orchestrate_workflows.py --status`

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No NCBI_EMAIL set | `export NCBI_EMAIL="your@email.com"` |
| ncbi-datasets-pylib missing | `uv pip install ncbi-datasets-pylib` |
| Slow SRA search | Normal - NCBI API, wait 5-10 minutes |
| FTP URL invalid | Manually verify assembly on NCBI, update YAML |

## Complete Documentation

See [GUIDE.md](GUIDE.md) for comprehensive guide.

## See Also

- **[GUIDE.md](GUIDE.md)**: Comprehensive discovery system documentation
- **[README.md](README.md)**: Discovery results and current status
- **[WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution
- **[ORCHESTRATION/README.md](../ORCHESTRATION/README.md)**: Orchestrator overview




