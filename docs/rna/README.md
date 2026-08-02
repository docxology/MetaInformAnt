# RNA-seq and Amalgkit workflow

MetaInformAnt provides a typed Python workflow layer over the installed
Amalgkit command line interface. The active contract is discovered from the
installed Amalgkit version; the local verified environment uses Amalgkit
0.16.32 and the commands:

`metadata → select → getfastq → integrate → quant → merge → wsfilter → finalize → sanity`

For public-read archive discovery and download, the workflow uses the ENA
metadata/download path exposed by the installed Amalgkit release. A metadata
record or download attempt is not analysis evidence: sample selection, read
validation, quantification, merge, filtering, finalization, and sanity checks
must all produce their expected outputs.

`cstmm` and `csfilter` are optional cross-species branches. They are not
silently substituted for the per-species finalization chain.

## Start here

- [Getting started](GETTING_STARTED.md)
- [Configuration](CONFIGURATION.md)
- [Workflow execution](workflow.md)
- [Command reference](amalgkit/commands.md)
- [Step reference](amalgkit/steps/README.md)
- [Validation protocol](VALIDATION.md)
- [End-to-end validation](END_TO_END_VALIDATION.md)
- [Mounted-data setup](EXTERNAL_DRIVE_SETUP.md)
- [Hymenoptera project documentation](../../projects/hymenoptera_amalgkit/README.md)

## Execution contract

```text
configuration
    ↓
metadata → select → getfastq → integrate → quant → merge
                                                    ↓
                                  wsfilter → finalize → sanity
                                                    ↓
                            evidence manifest + report + analysis
```

The wrappers perform three separate jobs:

1. validate configuration and installed-command options;
2. execute a named Amalgkit command with logs and a JSONL manifest;
3. inspect outputs without treating a directory or zero exit code alone as
   proof of a valid result.

## Real-data boundary

Large outputs must live outside Git. Set the active root explicitly before
running a real cohort:

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
```

The root is an input/output boundary, not a claim that every configured
species has data. Always record configured species, materialized species,
selected metadata rows, valid abundance files, completed steps, warnings, and
the exact tool versions in the evidence manifest.

## Analysis boundary

Quantification, merge, filtering, and finalization are computational pipeline
outputs. They do not by themselves establish biological significance. A
publication-ready analysis additionally requires explicit sample inclusion
rules, metadata provenance, orthology strategy, normalization rationale,
effect sizes and uncertainty, multiple-testing control, sensitivity analyses,
and result-derived figures and captions. The Hymenoptera project records
those requirements in its manuscript and evidence documents.

## Public API

The stable workflow surface is in
`metainformant.rna.amalgkit.amalgkit`, `metainformant.rna.engine.workflow`,
and `metainformant.rna.steps`. The command registry is the source of truth
for current subcommands; use `amalgkit --help` and the project verification
scripts to confirm the installed environment before a real run.
