from pathlib import Path
from metainformant.core.io import io

path = Path("test.json.gz")
print("Loading...")
data = io.load_json(path)
print(f"Data: {data}")
