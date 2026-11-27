# Public interface architecture

`metainformant.core.io` is the public façade for table and JSON helpers. The benchmark uses `write_delimited` and `read_delimited` with an explicit tab delimiter. The reader yields dictionaries keyed by the input headings and represents ordinary field values as strings. The writer defaults to an atomic write path. We exercise the successful write path here, with crash recovery and concurrent writes outside the experiment.

`metainformant.core.data.validation` checks typed values. Its schema helper checks required keys and types. The benchmark checks ordered columns separately, validates a converted integer count, enforces a nonnegative range, and tracks record identity within each table. The policy belongs to the caller, so the same I/O layer can serve tables with other field types or constraints.

The JSON output is written using `dump_json` and checked using `load_json`. The chart consumes the checked per-case summaries. This ordering lets readers audit the numeric summary and trace the plotted values to the realized acceptance decisions.
