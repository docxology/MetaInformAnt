### Visualization: Plots

Functions: `lineplot`, `heatmap`, `pairplot_dataframe`.

```python
import pandas as pd
from metainformant.visualization import plots

ax = plots.lineplot(None, [1, 3, 2, 5], label="series A")

df = pd.DataFrame([[1,2],[3,4]], columns=["x","y"])
grid = plots.pairplot_dataframe(df)
```


