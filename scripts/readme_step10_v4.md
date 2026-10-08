### Step 10. The unannotated fraction (eggNOG)

eggNOG annotation of the 6,171,602 representatives was summarized by cluster size
and completeness, then carried to the 50% family level with support filters.

**Result.** 48.29% of representatives are annotated, rising from 36.83% for
singletons to 90.53% for clusters of twenty-one or more. Proteins spanning a
whole short contig annotate better than complete open reading frames (65.1%
against 40.4% among singletons), so partial calls are not low-quality sequence.
At the family level, 76.03% of all families have no eggNOG-annotated member,
falling to 55.40% of the 1,020,958 supported families (565,611 families); the
value changes by less than 0.1 percentage points with a stricter support filter.

| File | Purpose |
|---|---|
| `scripts/step10_annotation_summary_n88.sh` | Annotated fraction by cluster size and completeness |
| `scripts/step11_interaction_and_supported_n88.sh` | Family-level unannotated fraction under support filters |

Output: `results/annotation_summary_n88.tsv`, `results/family_dark_fraction_n88.tsv`

