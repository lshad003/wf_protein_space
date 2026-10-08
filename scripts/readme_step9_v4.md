### Step 9. Open reading frame completeness

Prodigal flags whether each protein has both a start and a stop codon or is
truncated at a contig edge. Completeness was summarized overall and by cluster size.

**Result.** 46.89% of the 40,257,962 proteins are complete. Among 95%
representatives, completeness rises with cluster size: 11.58% for singletons,
32.98% for pairs, 61.16% for three to five, 81.99% for six to twenty and 86.11%
for twenty-one or more. Singletons are mostly contig-edge fragments, which
supports using the 2,067,011 supported clusters as the catalog.

| File | Purpose |
|---|---|
| `scripts/step9_orf_completeness_n88.sh` | Completeness by cluster size |

Output: `results/orf_completeness_n88.tsv`

