### Step 2. Clustering

Proteins were clustered with MMseqs2 at 95%, 90% and 50% amino-acid identity.

**Result.** 6,171,602 clusters at 95% identity, 5,352,641 at 90% and 2,918,630
families at 50%.

| File | Purpose |
|---|---|
| `scripts/step2_cluster.sh` | Database build and clustering at three thresholds |
| `scripts/step25b_membership_n88.py` | Cluster tables at 88 metagenomes |

Output: `catalog/db/LsPS_AA_{95,90,50}_cluster`, `results/clusters_{95,90,50}_n88.tsv`

