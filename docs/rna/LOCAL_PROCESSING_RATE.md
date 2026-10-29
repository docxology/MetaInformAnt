# Local processing rate and artifact inventory

Optimize the rate of **new restore-bound sample completions** while preserving
the frozen identity set, provenance and resumable inputs. Raw transfer growth,
quantifier activity and reused checkpoints measure different work. None alone
establishes a best achievable rate, scientific completion or release readiness.

## What unconfirmed means

The operational observer computes:

```text
current_lower_bound = producer_audited_at_start + new_restore_bound
supported_unconfirmed = eligible - current_lower_bound
frozen_total = eligible + held_prerequisites
```

The startup baseline is an aggregate producer audit. Newly counted successes
must have matching restore proof, provenance contract and abundance digest.
Reused checkpoints are reported separately and do not increase new-sample rate.
The aggregate remainder does not enumerate filesystem absence. A primary task
selection can contain both acquisition work and existing outputs needing
checkpoint verification; it is not the unconfirmed set.

The [2026-10-23 inventory reconciliation](../../output/hymenoptera_presence_audit_20261023/20261023T160213Z/reconciliation.md)
records 18,200 identities across 27 species. At that capture, 7,150 supported
identities had quant files matching the checked contract metadata. Only 61
supported runs received the bounded full portable quant validation; 7,089
remained incomplete or unvalidated in that unit. Another 71 had raw/partial
data without active quant, 10,781 had no expected active payload observed, five
remained live/unknown, and 193 retained prerequisite holds. These mutually
exclusive categories describe that timestamped capture, not current totals.
All identities and exact paths/reasons remain in its TSV/JSON receipts.

Presence is not quantification validity. Absence at expected active paths is
not biological zero, scientific exclusion or proof that a run never started.
Raw cleanup after validated quant can be expected: consult quant/provenance,
restore witnesses and reclamation journals together. Existing checkpoints
without active quant require a separate recovery assessment; never restore
them merely to complete an inventory.

## Processing decisions

| Evidence | Next local processing decision |
|---|---|
| Current-contract quant with content verification | Reuse through the normal journal/checkpoint path; count reuse separately |
| Quant metadata matches, content unchecked | Verify before relying on reuse; do not mark complete from presence |
| Quant mismatched, incomplete or failed | Preserve evidence and raw inputs; resolve the exact contract/error before retry |
| Checkpoint without active quant | Validate recovery custody and ownership before an authorized restore |
| Raw/partial without active quant | Preserve resumable parts; validate current mates and gzip/records before quant |
| No expected payload observed | Treat as an acquisition candidate after checking live ownership and prior attempts |
| Live/changed/unknown | Keep ambiguity; let the current owner progress and recapture later |
| Held prerequisite | Resolve the method/input/tool contract separately; retain the frozen accession |

At the inventory capture, the deployed classifier reported 193 holds: 144
submitted METAGENOMIC/PCR libraries titled 16S and 49 long-read
backend/reference/bootstrap prerequisites. A later parent-owned frozen-metadata
title screen identified 15 additional Apis GENOMIC/PCR libraries titled
`16Sseq of bee gut` that the METAGENOMIC-only check had missed. The corrected
operational skip policy is therefore **159 16S plus 49 long-read identities,
208 in total**, pending activation of the corrected classifier. The 193 count
remains the deployed runtime's recorded policy until a fresh source/runtime
binding verifies the change. All 18,200 frozen identities remain visible.
These are operational processing decisions, not scientific exclusions or
biological zeros. Increasing concurrency cannot resolve those method
boundaries; do not reinterpret them as host-transcriptome successes or silently
switch quantification backends.

## Measuring useful local rate

Use two completed observations with identical runtime, frozen manifest,
selection and journal bindings. Divide the difference in `new_restore_bound`
by capture-to-capture hours. Report reuse and failures beside it. If a completed
baseline is missing or the phase changed, the comparable completion delta is
unavailable; do not substitute a historical phase. Short file-growth windows
are retained growth lower bounds, not network traffic, capacity or completion
estimates.

The reviewed local profile at the inventory capture used a 48-identity rolling
frontier, size-band ordering, a 256 GiB compressed-size reservation estimate,
four active acquisitions/eight file requests, two record validators, four quant
slots/eight quant threads and one dedicated recovery slot. Frozen reservations
are planning estimates; physical free-space and extraction/staging reservations
remain separate checks. A larger frontier does not mean more active transfers.

Keep strict FASTQ record/mate validation, provenance checks and restore-bound
checkpointing enabled. Measure acquisition, validation, extraction, quant and
checkpoint occupancy before changing caps. Use comparable workloads and retain
negative controls; raising worker count alone does not establish improvement.
One short higher growth window cannot prove a causal speedup.

The local continuous CLI accepts `--quant-slots` from 1 through 8, also
configurable with `AMALGKIT_LOCAL_QUANT_SLOTS`; its default remains four.
The total Kallisto thread budget remains eight, so eight slots use one native
thread per sample. Invalid values, or more slots than task workers, refuse
before runtime controls are written. The selected slot count is recorded in
`active_runtime.json` and the effective resource profile in `active_scope.json`.
Quantification slots include fresh input validation and preparation as well as
native Kallisto execution; slot occupancy is not a native process count.

The bounded throughput trial uses eight whole-sample acquisitions, sixteen
file connections and eight quantification slots sharing the same eight-thread
budget. It retains two metadata connections, two record validators, one
dedicated recovery, the 48-task frontier and all existing validation,
reservation and deadline controls. Native fixture comparisons qualify this
profile's execution and replay; useful production improvement still requires
comparable new restore-bound completions per hour after activation.

## Change and activation boundary

Prepare throughput changes in an isolated GitVault checkout against the exact
deployed source manifest. Preserve the live producer and its source hashes.
Exercise real files, local services and native tools with bounded fixtures;
retain before/after evidence and exact errors. Obtain independent review of the
candidate rather than auditing work built by the same owner.

A candidate that passes isolated checks is not deployed. A later source/profile
cutover must satisfy the current ownership, frozen-input, preservation,
quiescence, writer-absence, source-review and launch gates. Preserve partials,
active controls, historical data and unrelated dirty work. Do not replay a
one-shot historical launch operator against a new candidate. Keep measuring
new completions after an admitted cutover before claiming a campaign gain.

See [generic acquisition](GENERIC_ACQUISITION.md) for worker contracts and
[RNA workflow validation](VALIDATION.md) for validation scope.
# Opt-in HD-first storage

The local launcher accepts `AMALGKIT_LOCAL_STORAGE_MODE=hd_first`. It disables
SSD FASTQ, quant-workspace and fasterq scratch overrides, selecting the existing
canonical bulk workspaces instead. Source, the operational progress database,
small runtime caches and the optional index cache remain on GitVault.
In continuous HD-first operation, temporary checkpoint restore verification also
uses the bulk volume's `.checkpoint-verify` directory. Hash/reference validation
and raw-reclamation ordering are unchanged; a symlinked restore workspace refuses.
The default `ssd_compute` retains the previous staging behavior. Invalid modes
refuse before cache directories or controls are created.

Inspect effective configuration without creating directories or invoking tools:

```bash
AMALGKIT_LOCAL_STORAGE_MODE=hd_first \
  bash output/hymenoptera_local_20261019/run.sh print-storage-profile
```

Supply the ordinary `AMALGKIT_DATA_ROOT` and `AMALGKIT_FROZEN_ROOT` bindings.
This inspection is configuration only; it does not verify mounts, native
qualification, available capacity, or deployment approval. No free-space floor,
validation rule, quantification deadline, retry history or thread budget changes.
An HD-first profile requires its own reviewed launch bindings before activation.
Compare full representative samples with the same total CPU budget; fixture
timing and high CPU usage do not establish higher validated completions/hour.
The continuous CLI rejects malformed worker/slot/frontier arguments before
loading the scientific dependency stack; valid invocations use the same engine.

## EHD-only runtime data

The opt-in `AMALGKIT_LOCAL_STORAGE_MODE=ehd_only` separates read-only software
from RNA runtime data. `METAINFORMANT_ROOT` names the existing software checkout;
`AMALGKIT_OPERATIONS_ROOT` names an existing, disjoint runtime directory below
the mounted `AMALGKIT_DATA_NAMESPACE`. That namespace must be the volume's
`Data` directory, with `AMALGKIT_DATA_VOLUME_UUID` matching the actual mounted,
unlocked volume. Both Python prefixes remain on the distinct source volume.

Supply existing Data paths for `AMALGKIT_DATA_ROOT`, `AMALGKIT_FROZEN_ROOT`,
`AMALGKIT_CONFIG_DIR`, `AMALGKIT_PROGRESS_DB_PATH`, the source and toolchain
manifests, and any size overlay. The runtime root contains the unchanged frozen
envelope at `frozen/manifest.jsonl`. Progress policy must be `DELETE`; selecting
this mode neither creates nor repairs a database. Qualified native software
continues to use `AMALGKIT_RUNTIME_DIR` and `AMALGKIT_TOOLCHAIN_BIN` on the
read-only source volume.

The launcher places HOME, TMP, uv/XDG/Matplotlib/Numba caches and the index cache
below the Data runtime root, clears scratch overrides, disables FASTQ staging
and source bytecode writes, and refuses mismatched or symlinked roots before
creating controls. Runtime receipts record both roots and the executing source
program. Configuration inspection uses the existing `print-storage-profile`
command and performs the EHD root checks without creating directories.

The source-side `sample_rate.py` accepts `--source-root` independently of
`--operations`. For EHD-only runtimes it requires `--stdout`; a separately
reviewed operator owns receipt publication. It rechecks source and manifest
bytes after sampling. The unchanged acquisition observer takes explicit Data
operations, frozen manifest, selection and journal paths. Pathname growth
remains a lower bound, without authenticated native-FD or network attribution.

Host process collection uses `host_sample.py` outside the worker's file-write
sandbox: macOS refuses the setuid process-query tool inside that sandbox.
Supply the exact sampler and `/bin/ps` SHA-256 bindings using
`--sampler-sha256` and `--ps-sha256`, together with `--operations`,
`--source-root`, `--manifest` and the bounded `--seconds` window. The collector
always requests stdout-only output, refuses Python file mutations and permits
only the pinned process-query command. Its Python audit hook is a guard on
reviewed code, not an OS sandbox for arbitrary native code. The operator owns
the inherited stdout destination. Worker and artifact-writer write restrictions
remain separate; moving host observation does not permit their writes outside
the selected Data workspace.

This mode preserves existing capacity guards, thread budgets, deadlines,
validation, checkpoint identity and attempt history. It supplies no migration
or launch authority. A deployment still requires source/profile/operator review,
quiescence, raw database/sidecar preservation and private integrity/backup checks.
Private native fixtures do not establish a representative full-sample thread
optimum, hardware power-loss durability or repair of prior storage faults.


### Acquisition startup dependencies

RNA subpackages load on first access. Importing the campaign acquisition worker
does not initialize analysis plots or platform font discovery. Existing public
RNA module names and explicit analysis imports remain available. This controls
startup overhead; it does not establish full-input throughput or optimal native
thread counts.


### Validated EHD capacity admission

The opt-in `ehd_only` mode validates the mounted Data UUID, ordinary control/bulk/config/progress roots, both external Python prefixes, and every managed cache/temp/home/index/log path before the engine opens logs or SQLite. SSD scratch and FASTQ staging are refused. This mode checks the actual EHD output filesystem and retains at least256GiB of external headroom; it has no RNA admission dependency on unused internal free space. The production launcher must also apply its reviewed native write sandbox. Legacy SSD/HD-first modes keep their existing system-space guard. The EHD root identities and path environment are rechecked before acquisition; path redirection, changed declarations or capacity telemetry failure refuse without accepting inputs. This concerns managed RNA writes; macOS may independently use internal storage for its own services.
