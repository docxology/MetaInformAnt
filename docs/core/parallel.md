### Core: parallel

Function: `thread_map`

```python
from metainformant.core import parallel

def double(x: int) -> int:
    return 2 * x

out = parallel.thread_map(double, [1, 2, 3], max_workers=4)
```


