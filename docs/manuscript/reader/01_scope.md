# Scope

This reader accompanies MetaInformAnt’s local table I/O [@metainformant]. It demonstrates the round trip from a deliberately synthetic table to Python records and back to a file. The table uses invented record labels and arbitrary nonnegative counts. It contains no sequences, observed measurements, sample identifiers, or accession numbers.

Tabular files carry syntax, such as column headings and delimiters, separately from meaning, such as the type or interpretation of a cell [@w3c_tables]. MetaInformAnt’s delimited reader uses the headings as dictionary keys and returns field values as strings, so the application must explicitly convert counts before calculating totals.
