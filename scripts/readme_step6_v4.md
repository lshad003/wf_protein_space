### Step 6. Sequencing batch

Sequencing run, centre, platform and read depth were tabulated for all 88
metagenomes to check cohort comparisons against batch.

**Result.** WF23 and WF24 share one run and yield similar protein numbers per
metagenome (422,857 and 431,504), so cohort alone does not change yield. WF24
spans two runs that differ 1.75-fold (431,504 against 756,692) on 1.24-fold more
reads, driven by assembly size; gene density per Mb is constant. WF22 was
sequenced on separate runs, so WF22 against WF23 and WF24 cannot be separated
from batch.

| File | Purpose |
|---|---|
| `scripts/step6b_batch_table_n88.sh` | Batch table and run contrasts |

Output: `results/batch_table_88.tsv`

