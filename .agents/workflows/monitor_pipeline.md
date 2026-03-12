---
description: how to monitor amalgkit pipeline progress safely
---
# Monitoring amalgkit Pipeline

For large-scale runs, direct database queries can hang. Use these steps to safely check progress.

1. **Check High-Level Status**
   Use the specialized monitoring script which uses safe SQLite patterns:
   ```bash
   python3 scripts/rna/check_pipeline_status.py
   ```

2. **Detailed Progress (Read-Only)**
   If you need to query the database directly, use a read-only URI:
   ```bash
   sqlite3 "file:output/amalgkit/pipeline_progress.db?mode=ro" "SELECT status, COUNT(*) FROM progress GROUP BY status;"
   ```

3. **Intensive Analysis (Snapshot)**
   For complex reporting (e.g., pivot tables), copy the DB to `/tmp` first:
   ```bash
   cp output/amalgkit/pipeline_progress.db /tmp/pipeline_status.db
   python3 -c "import sqlite3; conn = sqlite3.connect('/tmp/pipeline_status.db'); ..."
   ```

4. **Monitoring Logs**
   Check for recent failures in the orchestrator log:
   ```bash
   tail -n 100 output/amalgkit/run_all_species_incremental.log
   ```
