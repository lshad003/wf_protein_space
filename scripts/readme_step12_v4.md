### Step 12. Combined annotation and protein classes

eggNOG and Pfam (hmmsearch, gathering thresholds) were combined: a protein is
known if either source annotates it.

**Result.** The union annotates 3,222,677 of 6,171,602 representatives (52.22%).
Of the 2,067,011 supported representatives, 585,347 (28.32%) have no annotation.
At the 50% family level, the supported unannotated fraction is 51.47%, against
55.40% with eggNOG alone. Supported representatives fall into three classes:
K, with a Pfam domain (1,343,020, of which 76,515 carry only DUF or UPF domains);
KWP, eggNOG only (138,644); and U, unannotated (585,347).

| File | Purpose |
|---|---|
| `scripts/step12d_union_replevel_n88.sh` | eggNOG and Pfam union per representative |
| `scripts/step12e_family_union_n88.py` | Union at the 50% family level |
| `scripts/step12f_rep_classes_n88.sh` | K, KWP and U labels |

Output: `results/union_rep_level_n88.tsv`, `results/family_union_dark_n88.tsv`,
`results/rep_classes_n88.tsv`

