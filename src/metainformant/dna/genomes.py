from __future__ import annotations

import re

_ASSEMBLY_REGEX = re.compile(r"^GC[FA]_[0-9]{9}(?:\.[0-9]+)?$")


def is_valid_assembly_accession(accession: str) -> bool:
    """Return True if the string looks like a valid NCBI assembly accession.

    Examples: GCF_000001405.39, GCA_000001405
    """
    return bool(_ASSEMBLY_REGEX.match(accession))
