# Generic idempotent Amalgkit acquisition

The public API is `metainformant.rna.amalgkit.acquisition`. Shared methods live
in `metainformant.rna.engine.acquisition_*`; parent and project scripts are thin
adapters. The same manifest worker runs locally and on AWS. It stops after
per-sample quantification: downstream analyses require their own complete-cohort
and scientific evidence gates.

## Commands

```bash
uv run python scripts/rna/acquisition.py --help
uv run python scripts/rna/acquisition.py local --help
uv run python scripts/rna/acquisition.py aws --help
```

The command family provides `freeze`, `plan`, `estimate`, `quote-aws`, `local`
(also `worker`), and `aws`. The nested Hymenoptera project delegates through
`scripts/acquisition.py`; its existing `scripts/cloud/cloud_worker.py` remains
a compatibility adapter to the shared worker.

AWS collection, pricing, durable S3 access and controller execution require
the optional `aws` extra. Local SQLite/file status helpers and command discovery
do not import the AWS SDK. Local quantification additionally requires the `rna`
extra and working external tools.

## Freeze and allocate

Use existing Amalgkit-selected metadata, reference indexes and progress DB as
the source. `freeze` appends newly discovered ENA runs into an isolated campaign
snapshot, preserving canonical inputs. The generic command accepts any nonempty
configured species set. The legacy Hymenoptera freeze API retains its 27-species
default. Reusing a frozen directory with different configuration bytes fails.

```bash
uv run python scripts/rna/acquisition.py freeze \
  --data-root "$AMALGKIT_DATA_ROOT" --config-dir "$AMALGKIT_CONFIG_DIR" \
  --output-dir "$AMALGKIT_CAMPAIGN_ROOT"

uv run python scripts/rna/acquisition.py plan \
  --campaign-root "$AMALGKIT_CAMPAIGN_ROOT" --backend hybrid \
  --local-fraction 0.25 --cloud-status "$AMALGKIT_STATUS_JSON" \
  --output-dir output/acquisition/plan
```

`plan` also accepts an existing `--manifest`. A cloud-status snapshot is optional
but must be supplied when reconciling an active cloud campaign. Its receipt and
live-assignment IDs exclude completed and reserved work from new assignments.
Observation age matters: the AWS controller rechecks current ownership before
launching. Plans are deterministic count allocations, not sample-size balancing.

The immutable plan contains `allocation.json` and nonempty lane selections
`local_partition.json` / `aws_partition.json`. Repeating the same plan is
idempotent; changing assignments in the same directory fails. Each task belongs
exactly once to completed, reserved, local-pending or AWS-pending. Existing AWS
reservations retain AWS retry ownership. Manifest/snapshot hashes bind lane
selections; the AWS controller additionally checks the original inventory hash.
Stop both lanes and use a new coordinated allocation generation before
repartitioning; do not run independently generated overlapping local/AWS plans.
For a hybrid plan, start the allocated AWS coordinator first. The local worker
refuses acquisition until the persisted controller ledger acknowledges the same
allocation hash. The coordinator also refuses changed or omitted allocations on
resume; a plan proposal alone does not reserve tasks against a running controller.

## Local execution

```bash
uv run --extra aws --extra rna python scripts/rna/acquisition.py local \
  --manifest "$AMALGKIT_CAMPAIGN_ROOT/manifest.jsonl" \
  --task-selection output/acquisition/plan/local_partition.json \
  --data-root "$AMALGKIT_LOCAL_WORK_ROOT" --config-dir "$AMALGKIT_CONFIG_DIR" \
  --stage-inputs --workers 4 --threads 8 --quant-slots 2 --fastq-slots 1 \
  --fastq-threads 2 --compression-threads 2 --validation-slots 2 \
  --max-in-flight 8 --max-raw-bytes 64424509440 \
  --durable-bucket "$AMALGKIT_DURABLE_BUCKET" \
  --durable-cohort "$AMALGKIT_DURABLE_COHORT" --profile dev-agent --region us-east-2
```

Choose a writable, adequately sized work root. `--stage-inputs` copies frozen
metadata/index files without replacing different existing files. Configurations
are supplied explicitly and checked against generic-envelope checksums. Unknown
raw sizes cannot be admitted under a positive raw-byte bound; resolve source
metadata first. The shared AWS controller already consumes hash-bound cached
NCBI source resolutions for such tasks.

One process owns each data root exclusively; concurrent library invocations in
one process are refused because existing RNA paths use process-wide environment
configuration. Use separate worker processes and disjoint selections for
parallel lanes. Internal sample execution and quantification, FASTQ extraction,
compression, validation and submission-window controls use the existing shared
resource-profile builder, including its effective CPU limits.

Current local quantifications are reused. With durable storage configured,
reference-bound receipts are restored and validated before acquisition; corrupt
existing blobs fail rather than triggering a new download. Portable config/index
witnesses preserve original provenance. Each invocation retains a separate task
journal and result under `acquisition_runs/`; `quantified` is the sum of `reused`
and `newly_quantified`. Calibration must exclude reused samples. Local raw
reclamation is opt-in (`--reclaim-raw-after-quant`); AWS bootstrap opts into the
existing provenance-gated cleanup.

## AWS execution and limits

```bash
uv run --extra aws python scripts/rna/acquisition.py aws \
  --campaign-root "$AMALGKIT_CAMPAIGN_ROOT" --repo "$METAINFORMANT_REPO" \
  --config-dir "$AMALGKIT_CONFIG_DIR" \
  --task-allocation output/acquisition/plan/allocation.json \
  --bucket "$AMALGKIT_DURABLE_BUCKET" --cohort "$AMALGKIT_DURABLE_COHORT" \
  --profile dev-agent --region us-east-2 --ami "$AMALGKIT_WORKER_AMI" \
  --instance-profile "$AMALGKIT_INSTANCE_PROFILE" --instance-type c7i.2xlarge \
  --budget 750 --historical-gross 0 --max-workers 6 \
  --worker-workers 16 --worker-threads 8 --worker-quant-slots 4 \
  --worker-fastq-slots 1 --worker-fastq-threads 2 \
  --worker-compression-threads 2 --worker-validation-slots 4 \
  --worker-max-in-flight 12 --priority-species ""
```

This is an operator example, not authorization to create another campaign or
reset historical spend. Use the actual historical gross usage for the selected
ledger. IAM role, bucket permissions, regional quota, AMI and network access must
already support the workload. Existing live workers keep their admitted source,
request, deadline and price; new settings apply to subsequently admitted jobs.

Fallback FASTQ compression defaults to pigz level 6. Select level 1 for temporary
FASTQ files when CPU time matters more than compressed scratch size: pass
`--compression-level 1` to local/worker execution or
`--worker-compression-level 1` to the generic AWS controller. Both accept levels
1–9 and reject invalid values. Direct streaming callers can set
`AMALGKIT_PIPELINE_COMPRESSION_LEVEL=1`. Compression remains lossless; the existing
FASTQ validation and provenance checks still apply. Level 1 can require more
scratch space, so retain raw-byte reservations and disk headroom and measure
the workload before changing concurrency. Worker results and AWS job records
retain the selected compression level for cost and throughput comparisons.

The default bootstrap supports Amazon Linux 2023 on x86_64. Other architectures
need an explicit custom startup template; licensed AMIs are excluded from this
Linux pricing model. Burstable instance families use standard CPU credits to
avoid surplus charges. Regional Linux compute and baseline gp3 catalog prices,
one public IPv4 address, the existing hourly margin and operator floor bound each
job. The gross ceiling includes historical charge, elapsed job charges, active
reservations, the full proposed deadline and storage reserve before admission.
Credits are not deducted. Completion is not guaranteed within a configured cap.

Partition bytes/task counts, minimum/maximum disk, expansion factor and disk
reserve are configurable. Defaults remain 60 GiB raw reservations, 120 tasks
after initial admission, and 600–2,000 GiB disks. Larger configured disks are
repriced; the controller imposes its own 16,384 GiB maximum.

`--disk-throughput-mibps` sets gp3 throughput for new AWS jobs and price quotes
(default 125 MiB/s). Values 125–750 retain baseline 3,000 IOPS. The live AWS
catalog's GiB/s-month throughput rate is converted to MiB/s-month; throughput
above 125 is charged separately and included in the complete job reservation
before admission. Missing or invalid rates fail closed. Check the selected
instance's sustained EBS bandwidth before increasing the volume setting.
Existing admissions keep their original request and price.

Full-cohort receipt coverage and fleet drainage are required before restoring all outputs and
writing a completion certificate. An exhausted or unresolved lane does not
produce a completion claim.

## Cost and time scenarios

Read current prices without launching workers:

```bash
uv run --extra aws python scripts/rna/acquisition.py quote-aws \
  --profile dev-agent --region us-east-2 --instance-type c7i.2xlarge --disk-gib 600
```

Provide an explicit rates JSON. This example is illustrative, not a measured
throughput promise:

```json
{
  "local": {
    "evidence": {"observed_units": 4, "samples_per_hour_low": 8,
      "samples_per_hour_high": 12, "source": "replace with benchmark receipt"},
    "costs": {"hourly_usd_per_unit": 0}
  },
  "aws": {
    "evidence": {"observed_units": 6, "samples_per_hour_low": 30,
      "samples_per_hour_high": 40, "fleet_rate_cap": 50,
      "source": "replace with cohort receipt observation window"},
    "costs": {"hourly_usd_per_unit": 0.55, "fixed_usd": 10,
      "setup_hours": 0, "retries_per_task": 0}
  }
}
```

```bash
uv run python scripts/rna/acquisition.py estimate \
  --allocation output/acquisition/plan/allocation.json \
  --rates "$AMALGKIT_RATES_JSON" --local-units 4 --aws-units 6 \
  --spent-usd "$AMALGKIT_OBSERVED_GROSS" --reserved-usd "$AMALGKIT_ACTIVE_RESERVE" \
  --ceiling-usd 750 --output output/acquisition/estimate.json
```

Units are local sample-worker slots or AWS instances, respectively. Rate bounds
are scenario evidence, not confidence intervals. Similar sample/source mix is
assumed; different capacity is marked extrapolated. A fleet-rate cap models
saturation. Setup, retries and saturation can increase cost when parallelism
rises. Concurrent lane duration is the maximum lane duration; costs add. Local
zero cost is an explicit marginal-cost assumption and omits electricity and
equipment cost. S3 storage/requests/network and other charges outside the hourly
quote belong in `fixed_usd` or a separate reserve.

The estimate covers pending lane work. Existing reserved tasks remain an
explicit completion dependency; their reserved cost is counted separately and
must not also enter the pending-lane estimate. A conservative budget-fit flag is
a model result, not permission to launch. Benchmark scaling changes before
adopting them, and keep the runtime admission ceiling authoritative.
