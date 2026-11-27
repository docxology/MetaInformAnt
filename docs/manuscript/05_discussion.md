# Discussion

Testing both valid and invalid inputs makes the boundary between syntax and semantics visible. A duplicate identifier can be a valid delimited record and still violate a project’s identity constraint. Conversely, a correctly shaped table can carry values outside a domain range. Here these checks compose without requiring any provider, network access, or study-specific inputs.

This benchmark covers a regular family of small local tables with explicit types. It does not characterize performance, robustness under concurrent writes, or behavior on arbitrary formats and missing values. The checks establish software behavior on controlled inputs. They leave sampling quality, normalization decisions, biological interpretation, and broader package verification outside the evidence.
