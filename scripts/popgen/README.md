# Population Genetics Scripts

Thin orchestrator scripts for the population genetics workflow (generate -> analyze -> report -> visualize). All analysis business logic lives in [`src/metainformant/popgen/`](../../src/metainformant/popgen/); simulation generators in [`src/metainformant/simulation/models/popgen.py`](../../src/metainformant/simulation/models/popgen.py).

## Scripts

| Script | Purpose |
|--------|---------|
| `analysis.py` | End-to-end workflow entry point (CLI) |
| `analyze.py` | Analysis orchestration helper (delegates to `metainformant.popgen`) |
| `generate_dataset.py` | Multi-scenario synthetic dataset generation |
| `report.py` | Human-readable markdown report writer |
| `visualize.py` | Plot generation for analysis results |

## Usage

```bash
uv run python scripts/popgen/analysis.py --output-dir output/popgen
```

## Related

- [Popgen Module](../../src/metainformant/popgen/) - reusable analysis methods
- [DNA Population](../../src/metainformant/dna/population/) - core statistics
- [Simulation Models](../../src/metainformant/simulation/models/) - generators
