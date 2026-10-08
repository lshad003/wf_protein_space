### Step 7. Rarefaction

Clusters were accumulated against proteins sampled rather than metagenomes, with
WF24 split by sequencing run so that cohort and run are not conflated.

**Result.** At 3,600,000 proteins, WF22 and WF23 give 1,169,677 and 1,159,262
clusters (0.9% apart). WF24 gives 1,295,364 and 1,687,560 in its two runs (30.3%
apart), while WF23 and WF24 within one run differ by 11.7%. Sequencing run, not
cohort, drives richness. The pooled curve reaches 6,171,602 clusters at
40,257,962 proteins and still adds 436,425 clusters between 35 and 40 million
proteins, so the catalog is not saturated.

| File | Purpose |
|---|---|
| `scripts/step7_rarefaction_n88.sh` | Accumulation against proteins sampled |
| `scripts/step7b_wf24_byrun_n88.sh` | Richness at matched depth, WF24 by run |

Output: `results/rarefaction_95_n88.tsv`

