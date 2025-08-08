### Core: text

Functions: `normalize_whitespace`, `slugify`, `safe_filename`

```python
from metainformant.core import text

s = text.normalize_whitespace("  Hello   world \n")
slug = text.slugify("My Title: Hello World!")
fname = text.safe_filename("My Report.pdf")
```


