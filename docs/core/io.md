# Core: io

Functions: `ensure_directory`, `open_text_auto`, `load_json`, `dump_json`, `read_jsonl`, `write_jsonl`, `read_delimited`, `write_delimited`.

```python
from metainformant.core import io

# Default to writing under output/
io.ensure_directory("output/example")
io.dump_json({"a": 1}, "output/example/a.json")
rows = list(io.read_delimited("output/example/a.csv"))
```
