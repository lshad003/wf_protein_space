### Step 22. Community-level tests

Whole-community gene abundance (Bray-Curtis on Hellinger-transformed counts) and
gene presence (Jaccard, at least five reads) were tested with PERMANOVA. In WF22,
between-animal terms (treatment, egg mass) were tested by permuting whole animals
(42 samples, 14 animals) and on animal centroids (15 animals); month was permuted
within animals. A Pfam-defined carbohydrate subset (59,297 representatives, 7,891
tested) was analysed the same way and per gene.

**Result.**
- Treatment: R2 0.051 (p = 0.928) by whole-animal permutation and 0.069
  (p = 0.657) on centroids; presence R2 0.062 (p = 0.665) and 0.106 (p = 0.465).
  In WF24 (34 animals) R2 is 0.173 (p = 0.679) and 0.170 for presence
  (p = 0.770), at chance level for 6 of 33 degrees of freedom. Genes detected per
  sample do not differ (p = 0.3676).
- Egg mass: R2 0.327 by whole-animal permutation and 0.513 on centroids
  (both p = 0.001), with equal dispersions (p = 0.214). Month: R2 0.127 (p = 0.001).
- Carbohydrate subset: zero genes respond to treatment, against 1,541 and 2,046 for
  egg mass and 1,761 and 2,586 for month; composition gives treatment R2 0.053
  (p = 0.454) and egg mass 0.673 (p = 0.001).
- Treatment groups differ in dispersion (p = 0.003) but not in centroid.

| File | Purpose |
|---|---|
| `scripts/step22b_permanova_fix_n88.R` | Abundance composition, WF22 |
| `scripts/step22c_jaccard_n88.R` | Presence composition and richness, WF22 |
| `scripts/step22e_permanova_wf24_n34_n88.R` | Abundance and presence composition, WF24 |
| `scripts/step22d_cazy_subset_n88.R` | Carbohydrate subset |

Output: logs `step22b_n88` to `step22e_n88`, `results/limma_cazy_*_n88.csv`

