# Population Genetics Scripts

Thin orchestrator script for the population genetics workflow (generate -> analyze -> report -> visualize). All analysis business logic lives in [`src/metainformant/popgen/`](../../src/metainformant/popgen/); simulation generators in [`src/metainformant/simulation/models/popgen.py`](../../src/metainformant/simulation/models/popgen.py).

## Scripts

| Script | Purpose |
|--------|---------|
| `analyze.py` | End-to-end workflow entry point (CLI); delegates to `metainformant.popgen.workflow` |

## Usage

```bash
uv run python scripts/popgen/analyze.py --output-dir output/popgen
```

## Related

- [Popgen Module](../../src/metainformant/popgen/) - reusable analysis methods
- [DNA Population](../../src/metainformant/dna/population/) - core statistics
- [Simulation Models](../../src/metainformant/simulation/models/) - generators
