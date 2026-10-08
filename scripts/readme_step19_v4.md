### Step 19. Properties of the unknown fraction

Length, completeness, prevalence, persistence across years and abundance were
compared across four classes: K (Pfam domain), KWP (eggNOG only), GU (genomic
unknown) and EU (environmental unknown).

**Result.**
- Length: mean coding length is 1,045 bp for K, 795 for KWP, 498 for GU and 367
  for EU. Complete reading frames are 79.5%, 69.8%, 76.0% and 63.7%, so unknown
  genes are short but mostly complete.
- Prevalence: 82.77% of unannotated genes occur in more than 20 of 88 samples,
  3,496 in all 88, and 0.65% in five or fewer.
- Persistence: 496,266 unannotated genes occur in all three years and carry
  91.15% of unannotated reads; 18,896 single-year genes carry 1.23%, and 260
  genes receive no reads.
- Abundance: unannotated genes are 28.32% of supported genes, 7.04% of raw mapped
  reads and 19.32% after length normalization; the raw share is 7.11% in WF22,
  6.28% in WF23 and 7.16% in WF24 (per sample 4.38 to 16.12%).

| File | Purpose |
|---|---|
| `scripts/step16y_complete_by_class_n88.py` | Length and completeness by class |
| `scripts/step16v_prevalence_n88.py` | Prevalence by class |
| `scripts/step19a_core_unknown_n88.py` | Persistence across years |
| `scripts/step16w_dark_abundance_n88.py` | Abundance share, raw and length-normalized |
| `scripts/step18b_fig12.py` | Prevalence and abundance figures |

Output: `results/completeness_by_class_n88.tsv`, `results/prevalence_primary_n88.tsv.gz`,
`results/core_unknown_by_years_n88.tsv.gz`, `results/dark_abundance_by_sample_n88.tsv`

