from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_antwiki_json(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return [data]
    return []


