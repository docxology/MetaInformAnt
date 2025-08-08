### Core: cache

Functions: `cache_json`, `load_cached_json`

```python
from pathlib import Path
from metainformant.core import cache

cache_dir = Path("./.cache")
cache.cache_json(cache_dir, "runs/demo", {"ok": True})
obj = cache.load_cached_json(cache_dir, "runs/demo", ttl_seconds=3600)
```


