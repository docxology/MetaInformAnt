# Results

Across three input scales, all three valid tables preserved their ordered fields and matched the independent total. All twelve malformed tables were rejected for the intended reason. The smallest tables contained 32 rows, the middle tables 256 rows, and the largest tables 2048 rows. These are deterministic fixture counts, so confidence intervals and significance tests do not apply.

![Paired observed and expected count totals agree at all three input scales. A heatmap marks all four malformed controls rejected for each scale.](../../output/figures/table_contract.png)

The JSON receipt records the observed totals alongside the oracle totals and the observed verdict alongside the expected verdict. The TSV hashes identify the exact files examined. A change to the renderer cannot silently change those underlying records.
