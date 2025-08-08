### Visualization: Animations

Function: `animate_time_series`

```python
from metainformant.visualization import animations

fig, anim = animations.animate_time_series([0,1,3,2,5], interval_ms=150)
# anim.save("series.mp4")  # optional
```


