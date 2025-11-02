### Visualization: Plots

Functions: `lineplot`, `heatmap`, `pairplot_dataframe`.

```python
import pandas as pd
from metainformant.visualization import lineplot, pairplot_dataframe

ax = lineplot(None, [1, 3, 2, 5], label="series A")

df = pd.DataFrame([[1,2],[3,4]], columns=["x","y"])
grid = pairplot_dataframe(df)
```
