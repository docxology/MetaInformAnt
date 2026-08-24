# Current Amalgkit step reference

The pages below describe the commands and output contracts used by the
MetaInformAnt workflow. The installed CLI remains authoritative; verify help
output before adding a command-specific option.

| Order | Command | Purpose | Primary output |
|---:|---|---|---|
| 1 | [metadata](01_metadata.md) | Retrieve sample metadata | `work/metadata/metadata.tsv` |
| 2 | [dataset](02_dataset.md) | Initialize workspace and rule assets (preparation) | workspace/rule files |
| 3 | [select](03_select.md) | Apply `select_rules.tsv` | selected metadata |
| 4 | [getfastq](04_getfastq.md) | Retrieve reads | FASTQ files |
| 5 | [integrate](05_integrate.md) | Attach local read paths | integrated metadata |
| 6 | [quant](06_quant.md) | Quantify samples | `work/quant/<run>/abundance.tsv` |
| 7 | [merge](07_merge.md) | Merge abundances | merged expression tables |
| 8 | [wsfilter](09_wsfilter.md) | Within-species filtering | filtered expression tables |
| 9 | [finalize](10_finalize.md) | Finalize analysis inputs | final expression tables |
| 10 | [sanity](11_sanity.md) | Validate outputs | sanity report |

`dataset` is a workspace-preparation command, not an additional biological
analysis stage. The required per-species analysis chain is
`metadata → select → getfastq → integrate → quant → merge → wsfilter →
finalize → sanity`.

`csfilter` is documented separately in [10_csfilter](10_csfilter.md) because
it is an optional cross-species branch rather than a required per-species
step.

`cstmm` is also optional and belongs after independently validated
cross-species inputs; it is not part of the default per-species chain. The
`busco` and `rerun` commands are auxiliary current 0.16.60 commands and are
covered in the [command reference](../commands.md).

## Real-data rules

- A non-empty directory is not sufficient evidence of completion.
- Every table must be tied to its input metadata and configuration hash.
- Missing species and partial steps remain visible in the evidence manifest.
- Bulk data stays on `AMALGKIT_DATA_ROOT`; only compact review artifacts belong
  in the repository.
