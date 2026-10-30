# RNA recovery observations and preservation

RNA recovery distinguishes immutable software, actual native dependencies,
managed data, and retained historical evidence. A passing private fixture is
component evidence; it does not establish production database integrity,
permission to launch, scientific completion, or optimal threading.

## Process exit races

process_is_confirmed_inert() checks a captured PID and original psutil creation
time against current kernel state. It returns true only for current
NoSuchProcess, or two observations of a zombie with the same creation time.
Live, inaccessible and unconfirmed processes remain unresolved. A changed
creation time raises ProcessIdentityChangedError; callers must not adopt it.

```python
import psutil

from metainformant.rna.engine.process_observation import process_is_confirmed_inert

owner = psutil.Process()
captured_time = owner.create_time()
assert process_is_confirmed_inert(owner.pid, captured_time) is False
```

This primitive does not enumerate producers, descendants, controllers, working
directories, descriptors or mapped writers. Operational recovery still needs
the complete native scope. AccessDenied itself never proves absence.

## Raw preservation before SQLite

With exclusive admission, ordinary source/bulk/runtime locks and complete
writer absence, preserve raw current and original databases plus every
present -journal, -wal and -shm member. Retain full hashes and source metadata.
Unexpected, missing, redirected or changing members refuse. Never open an
original database merely to see whether its journal is recoverable.

Use isolated copies for full integrity checks, SQLite API backup, exact schema
and row/history comparison, and a read-only check of the actual backup
locator and hash. Recheck the complete live family before private SQL,
publication and component replay. Operational database rows are not the
frozen biological cohort denominator.

## Qualification and activation

- Preserve failed tests and consumed invocations. Component replay does not
  grant permission to replay a production one-shot.
- Verify mounted physical roots, UUIDs, boot, source, native executable/library
  closure, interpreter prefixes and actual environment at use time.
- Reboot invalidates old process and passive transport identities. A new
  transport requires explicit qualification and operational review; preserve
  an exact owned stale socket through a non-overwriting rename.
- Later STOP/cancellation wins, including directories and broken symlinks.
  Never archive an unowned or later pause.
- Verify frozen inputs/IDs, protected outputs, resumable parts, real file and
  directory fsync, publication and actual blob-directory readback.
- Admit at most one producer. Bind a unique new selection and journal against
  that launch's inventories. Separate startup, reuse and new restore-bound
  outputs.

Validated EHD-only operation checks its actual external write roots and
capacity. Managed RNA files must not stage on another disk or depend on an
unrelated internal-space requirement. Legacy modes retain their own guards.
Moving retained evidence or changing storage policy needs separate review;
these observation primitives perform neither.

## Native fixture deadlines

The private runner in tests/rna/recovery_fixture.py accepts explicitly owned
fixture commands. Timeout remains failure even when cleanup succeeds. Before
each SIGTERM it records and rechecks the exact descendant PID and creation
time. It does not use global or forced kills, and preserves native output and
terminal/timeout receipts. It is test infrastructure, not a production
termination operator.

Real tests exercise file changes, symlink and identity refusals,
live/zombie/reaped states, native inspection denial and nested subprocess
timeouts. POSIX process and macOS sandbox cases declare their platform
prerequisites before execution. A portable CI result cannot substitute for
native operator qualification. No check here establishes full-input
throughput or resolves historical storage, quantification, metadata or
FASTQ failures.

See [external data roots](EXTERNAL_DRIVE_SETUP.md),
[local observations](LOCAL_PROCESSING_RATE.md), and
[filesystem verification](../core/verification.md).
