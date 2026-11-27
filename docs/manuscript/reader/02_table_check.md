# Local table check

The fixture records counts for invented groups A, B, and C. The export procedure checks the ordered columns, record labels, and nonnegative integer counts. It writes a second TSV via MetaInformAnt’s table writer and checks that the parsed records are equal on rereading.

![Synthetic table counts. Group A has 4, group B has 7, and group C has 2.](../../output/figures/synthetic_table_counts.png)

The chart summarizes the fixture rather than representing a biological observation. A successful round trip establishes preservation of this table’s field values. Questions about sampling, normalization, and interpretation require independent evidence.
