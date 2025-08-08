### Core: hash

Functions: `sha256_bytes`, `sha256_file`

```python
from metainformant.core import hash

digest = hash.sha256_file("README.md")
raw = hash.sha256_bytes(b"hello")
```


