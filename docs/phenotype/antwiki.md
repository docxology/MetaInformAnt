### Phenotype: AntWiki

#### Loading AntWiki JSON Data

Function: `load_antwiki_json`

```python-snippet
from pathlib import Path
from metainformant.phenotype.data import antwiki

entries = antwiki.load_antwiki_json(Path("tests/data/phenotype/antwiki_dataset_sorted_final_01.json"))
```

#### Web Scraping AntWiki

Web scraper that extracts species data from AntWiki pages. Pages are fetched with `requests`; the optional `cloudscraper` package is not used, and a Cloudflare or other bot-protection challenge is not bypassed (a blocked request surfaces as an error). The site's `robots.txt` is honoured by default (see the robots.txt section below).

**Basic Usage:**
```python-snippet
from metainformant.phenotype.data.scraper import AntWikiScraper, load_scraper_config

# Load configuration
config = load_scraper_config()

# Initialize scraper
scraper = AntWikiScraper(config)

# Scrape single species
data = scraper.scrape_species_page("Camponotus_pennsylvanicus")
print(f"Species: {data['species']}")
print(f"Traits: {data['traits']}")
print(f"Measurements: {data['measurements']}")
```

`scrape_species_page` takes a page name or a full URL. It returns the dictionary described under Extracted Data, or `None` when this scraper instance has already scraped that URL or the page could not be parsed. It raises:

- `RobotsDisallowedError` when robots.txt forbids the page, or cannot be retrieved;
- `NetworkError` when the page cannot be fetched (HTTP 404, or a 5xx response that persists through the retries);
- `requests.HTTPError` for HTTP 403 and other 4xx responses.

**Scraping Several Species:**
```python
from pathlib import Path

# Scrape named species; failures are reported instead of aborting the run
stats = scraper.scrape_all_species(
    output_dir=Path("output/phenotype/antwiki/"),
    species=["Camponotus_pennsylvanicus", "Formica_rufa"],  # page names or full URLs
    limit=10,  # Optional: scrape at most the first 10
)

print(f"Completed: {stats['completed']}, Failed: {stats['failed']}")

# Retry only the species that failed
if stats["failed_species"]:
    retry = scraper.scrape_all_species(
        output_dir=Path("output/phenotype/antwiki/"),
        species=stats["failed_species"],
    )
```

- `species` is required. AntWiki's species index is not retrieved, so omitting it raises `NotImplementedError`, as do `scrape_species_list()` and `get_species_list()`.
- The returned statistics hold `total` (species attempted, after removing duplicates and applying `limit`), `completed`, `failed`, and `failed_species`: the page name (or the full URL, for pages outside `base_url`) of every species that failed because robots.txt forbids it (or cannot be retrieved), the page could not be fetched, or the page could not be parsed.
- Each species is written to `<output_dir>/species/<name>.json`. `<output_dir>/all_species.json` holds the results of the latest call and is rewritten on every call. There is no checkpoint or resume: a repeated call fetches every named species again.
- `output_dir` defaults to `config.output_dir`; a `ValueError` is raised when neither is set.
- `bulk_scrape_species(species_urls, output_dir)` is the lighter alternative: it returns the list of scraped dictionaries, logs and skips failures (it reports no `failed_species`), and writes `<output_dir>/<name>.json`.

**Command-line Interface:**
```bash
# Scrape single species
python3 scripts/phenotype/scrape_antwiki.py --species Camponotus_pennsylvanicus

# Custom delay and output directory
python3 scripts/phenotype/scrape_antwiki.py --species Formica_rufa --delay 3.0 --output output/phenotype/antwiki/

# Scraper settings from a JSON file (keys are AntWikiScraperConfig keyword arguments)
python3 scripts/phenotype/scrape_antwiki.py --species Formica_rufa --config scraper_settings.json
```

`--species` writes `<output>/species/<name>.json` (default output directory `output/phenotype/antwiki/`) and exits with status 1 when robots.txt forbids the page (or cannot be retrieved), the page cannot be fetched, or it cannot be parsed. `--all` needs the AntWiki species index, which is not retrieved: it prints an error and exits with status 2, so name each species with `--species` (or call `scrape_all_species(species=[...])` from Python). `--resume` is rejected as unsupported. robots.txt is checked by default; `--check-robots` forces it on, and the opt-out is `"check_robots": false` in the `--config` file.

**Extracted Data:**

`scrape_species_page` returns a dictionary with the keys below. Values are matched by keyword in the page text under headings, so a field is empty or absent when the page has no matching heading.

- `url`, `data_source` (`"antwiki_scraper"`), `scraped_at` (Unix time) and `confidence_score` (a completeness score from 0 to 1)
- `species` and `species_name` (from the page title, `Genus species - AntWiki`, falling back to the page name when one was given) and `genus` (the first word of the title)
- `taxonomy`: `subfamily` and `tribe`, read from page elements whose class mentions taxonomy or classification
- `measurements` (the same dictionary as `morphology`): `body_length_mm`, `head_width_mm` and `color`, read from the text under headings that mention morphology, size, color, head or body. A size is a number followed by `mm` or `cm` and is stored as written, without unit conversion.
- `behavior`: `foraging_strategy`, `nest_type` and `colony_size_category`, from headings that mention behavior, foraging, nest or colony; `traits` is the list of these string values
- `ecology`: `habitat`, from headings that mention ecology, habitat or environment
- `distribution`: `regions`, from headings that mention distribution, range or location
- `description`: the first three paragraphs of the page body
- `images`: up to ten image URLs whose address contains `.jpg`, `.jpeg`, `.png` or `.gif`

`morphology`, `behavior` and `ecology` are present only when something was extracted.

**Configuration:**

`AntWikiScraperConfig` takes these keyword arguments:

| Setting | Default | Meaning |
|---------|---------|---------|
| `base_url` | `https://www.antwiki.org/wiki/` | Prefix that species page names are appended to |
| `delay_seconds` (alias `delay_between_requests`) | `2.0` | Minimum seconds between requests; values below 0.1 raise `ValidationError` |
| `timeout_seconds` (alias `timeout`) | `30` | Request timeout in seconds |
| `max_retries` | `3` | Retries of failed requests |
| `user_agent` | `MetaInformant/1.0 (research@metainformant.org)` | Sent with every request; robots.txt groups are matched against its leading product token (the text before the first `/`), so keep the token that identifies the scraper first |
| `check_robots` | `True` | Consult robots.txt (see below) |
| `output_dir` | `None` | Default output directory of `scrape_all_species` |

`load_scraper_config()` returns these defaults. Given a path, it reads a JSON file whose keys are the keyword arguments above:

```json
{"delay_seconds": 3, "check_robots": false, "output_dir": "output/phenotype/antwiki/"}
```

```python-snippet
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.phenotype.data.scraper import AntWikiScraperConfig, load_scraper_config

config = load_scraper_config("scraper_settings.json")

# config/phenotype/antwiki_scraper.yaml holds the same settings as a YAML template
config = AntWikiScraperConfig(**load_mapping_from_file("config/phenotype/antwiki_scraper.yaml"))
```

`load_scraper_config` reads JSON only: a YAML file raises `metainformant.core.utils.errors.IOError`, a missing file raises `ValidationError`, and no environment variable overrides a setting. `create_default_scraper_config(path)` writes the default settings as JSON that `load_scraper_config` reads back.

**robots.txt:**

With `check_robots=True` (the default), the scraper downloads `<origin>/robots.txt` once per origin (scheme, host and port) and per scraper instance, and applies it to every request, including each redirect hop:

- A URL that robots.txt disallows for the configured user agent raises `RobotsDisallowedError` (a `NetworkError` subclass) before any request for that URL is sent.
- The wait between requests is the larger of the configured delay and the site's `Crawl-delay`; fractional `Crawl-delay` values are rounded up to whole seconds.
- A robots.txt that cannot be retrieved (HTTP 5xx, 408 or 429, or a network error, after `max_retries` retries) makes the scraper refuse to proceed. HTTP 401 and 403 disallow everything. Any other 4xx status is read as no robots.txt, which allows everything.

`load_robots_txt(robots, site_url=None)` installs rules for an origin up front, from the file's text, its lines or a parsed `urllib.robotparser.RobotFileParser`, so that origin's robots.txt is never downloaded. `can_fetch(url)`, `crawl_delay(url)` and `effective_delay(url)` report the policy in force. `check_robots=False` is an explicit opt-out: robots.txt is neither requested nor applied. Parsing uses `urllib.robotparser`, so `*` and `$` wildcards inside paths are not understood, and rules are matched first to last instead of by longest match.

```python-snippet
from metainformant.phenotype.data.scraper import AntWikiScraper, AntWikiScraperConfig, RobotsDisallowedError

scraper = AntWikiScraper(AntWikiScraperConfig(delay_seconds=2.0))

# Supply robots.txt up front (a saved copy, for example); it is then never downloaded for this origin
scraper.load_robots_txt("User-agent: *\nDisallow: /wiki/Special:\nCrawl-delay: 5\n")

url = scraper.config.base_url + "Formica_rufa"
print(scraper.can_fetch(url))        # True
print(scraper.effective_delay(url))  # 5.0: the Crawl-delay exceeds the configured 2.0 seconds

try:
    scraper.scrape_species_page("Special:Random")
except RobotsDisallowedError as error:
    print(error)  # robots.txt of https://www.antwiki.org disallows https://www.antwiki.org/wiki/Special:Random ...
```

**Rate Limiting:**
After each fetched page, each redirect hop and the first robots.txt download for a site, the scraper waits the larger of the configured delay (default 2 seconds; `--delay` on the command line) and the site's `Crawl-delay`. Failed page requests are retried up to `max_retries` times with exponential backoff (1, 2, 4, ... seconds), and a backoff is never shorter than the `Crawl-delay`.
