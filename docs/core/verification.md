# Verification within one filesystem observation

```python
import os
from pathlib import Path

from metainformant.core.io.verification import BoundReads

source = Path(os.environ["VERIFICATION_INPUT"]).resolve(strict=True)
with BoundReads() as reads:
    expected = reads(source)
    assert reads(source) == expected
    bytes_read = reads.bytes_hashed
```

BoundReads computes a full SHA-256 digest once per absolute physical pathname
in one operation. Repeated uses check device, inode, mode, ownership, size,
mtime and ctime. Successful context exit rechecks every path, including paths
read only once. A new instance always reads bytes again.

VerificationChangedError refuses missing paths, directories, symlinks,
relative paths and detected changes. Descriptor reads use O_NOFOLLOW.
Opening or reading failures propagate. Existing caller exceptions remain
visible instead of being replaced by an exit-check error.

This observation does not replace ordinary filesystem locks, prevent changes
after its final check, prove database integrity, or grant deployment
permission. Recovery still needs current source/native/root bindings, writer
absence, raw database and sidecar preservation before SQLite, and independent
operational review.

See [I/O](io.md), [hashing](hash.md), and [RNA recovery](../rna/recovery.md).
