# Durable RNA quantification

The completion workflow freezes the configured species inventory, validates existing
quantifications, and processes verified missing accessions through bounded AWS jobs.
Each successful sample is preserved individually before replacement work is admitted.
Completion requires a full restored-output verification, rather than an instance exit
code, a journal count, or an end-of-run archive.

## Stored evidence

`metainformant.rna.engine.durable_quant` exposes `DirectoryStore`, `S3Store`,
`lock_quantification`, `restore_quantification`, and `validate_quantification`.
The AWS SDK is available through the `aws` dependency extra.

A locked sample requires the current Amalgkit runtime and contract, a matching
accession and species, the recorded abundance checksum, unique feature identifiers,
finite non-negative numeric values, positive finite count and TPM totals, positive
lengths and effective lengths, and positive integer Kallisto read counters.
Pseudoaligned reads cannot exceed processed reads; an optional target counter must
match the unique abundance features. Boolean and fractional read counters fail.
Fractional estimated counts and zero-expression rows remain valid; neither exact
estimated-count/read-count equality nor exact TPM normalization is required. Production callers additionally bind the
expected species configuration checksum.

Output files are stored by SHA-256 under `blobs/sha256/`. A sample receipt identifies
its files, sizes, hashes, contract, and cohort. Receipts are published only after
blob readback; local publication uses a same-filesystem hard link and S3 publication
uses conditional writes and server-side checksums. Conflicting output cannot replace
an existing receipt. S3 versioning adds recovery history; it is separate from legal
retention or S3 Object Lock.

Restoration checks every blob before publishing a fresh sample directory. The
portable validator retains original source paths as provenance and verifies output
content without requiring those historical machine paths to exist.

## Prepare one frozen cohort

Choose separate canonical and completion roots. Preparation reads canonical metadata,
indexes, exclusions, and progress without changing the canonical campaign. The
inventory appends newly discovered RNA-seq accessions within the configured taxa;
it preserves existing rows and their batch positions. It records ENA source tables,
configuration and input hashes, and explicit permanent exclusions.

```bash
export AMALGKIT_CANONICAL_DATA_ROOT=/path/to/canonical-data
export AMALGKIT_CAMPAIGN_ROOT=/path/to/completion-data
export AMALGKIT_DURABLE_COHORT=my-frozen-cohort
export AMALGKIT_BUCKET=my-authorized-bucket
export AWS_PROFILE=my-authorized-profile
export AWS_DEFAULT_REGION=us-east-2

uv run --extra aws --extra rna python - <<'PY'
import os
from pathlib import Path
from metainformant.rna.engine.completion_inventory import freeze_inventory, seal_existing_outputs

config = Path("projects/hymenoptera_amalgkit/config/amalgkit")
canonical = Path(os.environ["AMALGKIT_CANONICAL_DATA_ROOT"])
campaign = Path(os.environ["AMALGKIT_CAMPAIGN_ROOT"])
inventory = freeze_inventory(canonical, config, campaign)
seal_existing_outputs(
    inventory, canonical, config, campaign,
    bucket=os.environ["AMALGKIT_BUCKET"], cohort=os.environ["AMALGKIT_DURABLE_COHORT"],
    profile=os.environ["AWS_PROFILE"], region=os.environ["AWS_DEFAULT_REGION"],
)
PY
```

Keep the same root and cohort identifier when resuming. Files with missing provenance,
bad hashes, invalid values, or conflicting receipts remain unresolved. Public metadata
without a read object or a size reservation does not become a completed sample.

`legacy_quant_recovery.recover_indexed_archives` can recover complete quant members
from previously indexed truncated tar objects. Reads are bounded by member offsets
and bound to the source object's ETag. Recovery writes an isolated local store and
the same S3 receipt namespace; it does not promote an old live-source snapshot into
the canonical data root or manuscript evidence.

## Run bounded AWS processing

Launch only after recovery, input checks, and an authorized gross cost envelope.
Check the regional On-Demand vCPU quota and other running resources before raising
`--max-workers`; the fleet cap is an explicit operator setting, not a quota request.
Use an EC2 instance profile with access to the selected bucket. The controller checks
current regional compute and gp3 storage prices, reserves the complete job duration
plus a storage allowance, and keeps a durable local ledger. The hourly bound is the
maximum of the operator floor and compute plus provisioned gp3 storage divided by
672 hours (the shortest calendar month), one public IPv4 address at $0.005/hour,
and a $0.05/hour operating margin. Larger disks therefore increase reservations.
Each new job preserves its admitted hourly bound; earlier jobs retain the legacy
ledger rate. Termination requests continue accruing charges until EC2 termination
is observed. Credits never reduce this gross usage calculation.

```bash
uv run --extra aws --extra rna python scripts/rna/complete_hymenoptera.py \
  --campaign-root "$AMALGKIT_CAMPAIGN_ROOT" --repo "$PWD" \
  --bucket "$AMALGKIT_BUCKET" --cohort "$AMALGKIT_DURABLE_COHORT" \
  --profile "$AWS_PROFILE" --region "$AWS_DEFAULT_REGION" \
  --budget "$AUTHORIZED_TOTAL_GROSS_USD" --historical-gross "$VERIFIED_PRIOR_GROSS_USD" \
  --ami "$VERIFIED_AMI_ID" --instance-profile "$AUTHORIZED_INSTANCE_PROFILE" \
  --instance-type c7i.2xlarge --hourly-upper-bound 0.55
```

`--once` reconciles state and admits at most one job, then returns. Normal operation
continues through disjoint missing-task partitions. `--max-workers` defaults to
one and permits a bounded concurrent fleet when explicitly increased. Each active
worker retains task ownership until observed termination, so delayed receipts
cannot cause duplicate admissions. Every new admission reserves its full deadline
plus all outstanding worker deadlines under the same gross ceiling. Terminating
workers retain a shutdown allowance. If the remaining envelope cannot admit more
work, existing workers remain supervised until their results are reconciled.
Once all receipts exist, active workers are drained before final verification. The local controller should run
inside a persistent terminal session. Its ownership lock rejects a second controller
for the same root; idempotent launch tokens recover an uncertain API response.

Partitions have task and declared-byte limits. Unknown byte sizes are not admitted.
An individual oversized run is isolated, receives an appropriately sized disk and a
transfer-sized deadline, and must still fit the gross budget. The worker honors
explicit `max_in_flight`, reserves declared raw bytes, preserves the caller's free-space
floor, and bounds acquisition, extraction, and quantification by the job deadline.
The rendered startup script uses disk-backed temporary files and independent shutdown
deadlines. OS shutdown terminates the instance even if EC2 termination permission is
unavailable. Only verified quantification files are necessary for durable recovery;
failed raw reads are never packaged into the scientific result archive.

The controller stops on budget exhaustion or when remaining tasks exhaust their
bounded attempts or lack source-size evidence. These are unresolved outcomes, not
implicit scientific exclusions. Credits do not increase the authorized gross budget.

## Completion and scientific use

`aws_controller.json` records eligible, locked, and missing task counts, instance
ownership, deadlines, immutable source/input bindings, and the conservative spend
estimate. Per-sample failures remain in the worker journals and preparation ledgers.

Once all expected receipts exist, `verify_locked_campaign` restores and validates
every eligible sample under `completed_quant/`. It refuses empty or incomplete
inventories and configuration drift. Only then does it write
`quant_completion_certificate.json` with `all_quant_locked=true`.

That certificate establishes quantification storage and integrity. It does not replace
metadata harmonization, current downstream `merge → wsfilter → finalize → sanity`,
the finalized-matrix manifest, or the biological and manuscript release gates.

## Independent frozen-index binding

Production seals, worker reuse and the final certificate require
`expected_reference_index_sha256` from the frozen inventory. Locking verifies
that the provenance-bound complete reference manifest selects actual index
bytes with that digest. The verified manifest is stored as a content-addressed
blob; the bound receipt records both manifest and index hashes.

Bound receipts use `cohort/reference-bound-receipts/species/accession.json`.
Earlier `cohort/receipts/` objects remain immutable recovery artifacts and do not
satisfy the completion gate. Available original output/manifest/index files can
be resealed into the bound namespace without overwriting those objects. Missing
legacy reference evidence requires reacquisition/reprocessing or a separately
reviewed migration; no inferred index identity substitutes for content hashes.
Archive recovery now requires the frozen index hash explicitly and refuses
binding when the original manifest/index evidence is unavailable. Restore checks
independent expected configuration and index hashes before exposing outputs.

## Module boundaries

`quant_storage` owns append-only object storage, `quant_validation` owns output
integrity, and `durable_quant` owns receipt publication and restoration. Existing
imports from `durable_quant` remain supported. `aws_inputs` owns immutable input
archives and startup rendering; `aws_resources` owns price validation and elapsed
usage accounting. The controller retains the single-writer admission/reconciliation
state machine so resource reservation and persisted launch identity stay atomic.
Worker startup installs the RNA and AWS extras from the committed frozen lock.
