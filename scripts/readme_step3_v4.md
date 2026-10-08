### Step 3. Overlap with the earlier WF22 catalog

A random subsample of WF23 and WF24 proteins was searched with MMseqs2 against the
earlier WF22-only catalog, using query-side coverage. A positive control of WF22
proteins was recovered at 98.25%.

**Result.** 45.89% of WF23 and WF24 proteins match the earlier catalog at gene
level and 66.49% at family level.

| File | Purpose |
|---|---|
| `scripts/step2c_diagnostic_epyc.sh` | Subsample search, gene level |
| `scripts/step2d_control_and_tiers_n88.sh` | Positive control and family level |
| `scripts/step2e_covmode2_n88.sh` | Query-side coverage search |

