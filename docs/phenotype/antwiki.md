### Phenotype: AntWiki

Function: `load_antwiki_json`

```python
from pathlib import Path
from metainformant.phenotype import antwiki

entries = antwiki.load_antwiki_json(Path("tests/data/phenotype/antwiki_dataset_sorted_final_01.json"))
```
