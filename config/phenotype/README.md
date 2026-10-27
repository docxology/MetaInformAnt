# PHENOTYPE

## Overview
Configuration assets for Phenotype workflows. Keep templates and examples synchronized with the shared loaders in `metainformant.core.utils.config`.

## Contents
- [antwiki_scraper.yaml](antwiki_scraper.yaml)

## Usage
```python
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.phenotype.data.scraper import AntWikiScraperConfig

mapping = load_mapping_from_file("config/phenotype/antwiki_scraper.yaml")
config = AntWikiScraperConfig(**mapping)
```

The keys of `antwiki_scraper.yaml` are `AntWikiScraperConfig` keyword arguments. `load_scraper_config` in `metainformant.phenotype.data.scraper` reads JSON files only, and no environment variable overrides these settings, so build the config from the loaded mapping as above.
