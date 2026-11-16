# Projects

Domain projects and templates built on `metainformant`.

## Included Public Template

- [`template_bioinformatics_project/`](template_bioinformatics_project/) — reproducible bioinformatics project scaffold.

## Optional Private Extensions

Two research campaign repositories are optional private extensions and are omitted from `.gitmodules` so `git clone --recursive` of this public repository succeeds without credentials:

- `projects/hymenoptera_amalgkit` (`docxology/hymenoptera_amalgkit`) — optional private Hymenoptera Amalgkit RNA-seq campaign configs and orchestration scripts.
- `projects/apis_gwas` (`docxology/apis_gwas`) — optional private *Apis mellifera* BeeWAS GWAS campaign workflows and analyses.

Authorized users with repository access can check them out locally alongside the public tree (both paths are ignored by `.gitignore`):

```bash
git clone https://github.com/docxology/hymenoptera_amalgkit.git projects/hymenoptera_amalgkit
git clone https://github.com/docxology/apis_gwas.git projects/apis_gwas
```
