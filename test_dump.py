from pathlib import Path
from metainformant.core.io import io

path = Path("test.json.gz")
io.dump_json({"a": 1}, path)
print(f"File created: {path.exists()}")
with open(path, "rb") as f:
    print(f"File content: {f.read()}")
