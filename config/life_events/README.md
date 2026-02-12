# Life Events Configuration

Configuration templates for the life events analysis module (life course and event sequence analysis).

## Files

- `life_events_template.yaml` — Template with embedding, model, workflow, and output settings

## Usage

```bash
python3 scripts/life_events/run_life_events_analysis.py \
  --input data/life_events/sequences.json \
  --config config/life_events/life_events_template.yaml
```

## Related

- [Life Events Module](../../src/metainformant/life_events/)
- [Life Events Script](../../scripts/life_events/)
