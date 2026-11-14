# Durable RNA quantification

Validated per-sample quantification outputs can be locked into a
content-addressed local store and restored later with every byte re-verified.
A locked receipt is storage and integrity evidence for one sample; it is not a
completed cohort and does not replace the downstream release gates.

## Stored evidence

Three modules divide the work:

- `metainformant.rna.engine.quant_storage` owns append-only object storage:
  `DirectoryStore` (a local store that publishes by same-filesystem hard link),
  the `ObjectStore` protocol, and `receipt_key()`.
- `metainformant.rna.engine.quant_validation` owns output integrity:
  `validate_quantification()`.
- `metainformant.rna.engine.durable_quant` owns receipt publication and
  restoration: `lock_quantification()`, `restore_quantification()` and
  `bound_receipt_key()`.

A locked sample requires the current Amalgkit runtime and contract, a matching
accession and species, the recorded abundance checksum, unique feature identifiers,
finite non-negative numeric values, positive finite count and TPM totals, positive
lengths and effective lengths, and positive integer Kallisto read counters.
Pseudoaligned reads cannot exceed processed reads; an optional target counter must
match the unique abundance features. Boolean and fractional read counters fail.
Fractional estimated counts and zero-expression rows remain valid; neither exact
estimated-count/read-count equality nor exact TPM normalization is required.
Production callers additionally bind the expected species configuration checksum.

Output files are stored by SHA-256 under `blobs/sha256/`. A sample receipt identifies
its files, sizes, hashes, contract, and cohort. Receipts are published only after
blob readback. Conflicting output cannot replace an existing receipt.

Restoration checks every blob before publishing a fresh sample directory. The
portable validator retains original source paths as provenance and verifies output
content without requiring those historical machine paths to exist.

```python
from pathlib import Path

from metainformant.rna.engine.durable_quant import lock_quantification, restore_quantification
from metainformant.rna.engine.quant_storage import DirectoryStore

store = DirectoryStore(Path("/path/to/receipt-store"))
receipt = lock_quantification(
    store, "my-frozen-cohort", "apis_mellifera", "SRR12345678",
    Path("/path/to/data/apis_mellifera/work/quant/SRR12345678"),
    expected_config_sha256=config_sha256,
    expected_reference_index_sha256=frozen_index_sha256,
)
restore_quantification(
    store, "my-frozen-cohort", "apis_mellifera", "SRR12345678",
    Path("/path/to/restored/SRR12345678"),
    expected_config_sha256=config_sha256,
    expected_reference_index_sha256=frozen_index_sha256,
)
```

## Frozen-index binding

With `expected_reference_index_sha256`, locking verifies that the
provenance-bound complete reference manifest selects actual index bytes with
that digest. The verified manifest is stored as a content-addressed blob; the
bound receipt records both manifest and index hashes.

Bound receipts use `cohort/reference-bound-receipts/species/accession.json`
(`bound_receipt_key()`); receipts locked without an index digest use
`cohort/receipts/species/accession.json` and never satisfy a bound lookup.
No inferred index identity substitutes for content hashes. Restore checks
independent expected configuration and index hashes before exposing outputs.

`scripts/rna/import_quantifications.py` can read a reference-bound receipt
store as import evidence (`quant_import.ReceiptStore`); see the
[engine README](../../src/metainformant/rna/engine/README.md).

## Portable downstream restoration

`restore_quantification(..., verified_config_path=path,
expected_reference_index_sha256=frozen_index_hash)` verifies the supplied frozen
configuration bytes and stored reference manifest, then installs separate local
input witnesses below `.metainformant_restored_inputs`. Original quant sidecars
remain byte-identical. The current-method classifier verifies the witness schema,
contract identity, local path containment and configuration/reference checksums;
corrupt or mismatched witnesses fail closed. Without an explicit witness it keeps
the original input-path checks. This supports downstream work after a source
checkout moves, without changing the original command or pretending
quantification was rerun.

An isolated diagnostic subset must identify its selected run IDs and full frozen
species denominators. A complete subset is not a complete species. Keep its data
root separate from the active acquisition root while validating merge, wsfilter,
finalize, sanity and comparative methods.
