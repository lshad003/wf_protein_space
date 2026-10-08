### Step 17. Treatment, egg mass and month: per-gene tests

Gene abundance was tested per gene with limma-voom (TMM normalization), each cohort
separately. WF22 samples are repeated measures of 15 animals, so animal was
blocked with duplicateCorrelation.

**Result.** Within-animal correlation is 0.2864; without blocking, WF22 treatment
contrasts return 5,005 and 7,263 genes, and with blocking 3 and 155. In the same
blocked model egg mass returns 57,120 and 51,808 genes and month up to 141,580.
Within the balanced egg mass 3 subset, treatment returns zero genes. WF23 returns
zero; WF24 (34 animals) returns zero for five of six strains and 3 genes for one.
The unannotated share of gene abundance does not differ by treatment in any cohort
(all p > 0.1).

| File | Purpose |
|---|---|
| `scripts/step17c_wf22_meta.sh` | WF22 design table |
| `scripts/step17g_limma_wf22_n88.R` | WF22, blocked and unblocked |
| `scripts/step17h_limma_wf2324_n88.R` | WF23 |
| `scripts/step17i_limma_wf24_n34_n88.R` | WF24, 34 animals |
| `scripts/step17b_class_shift_n88.py` | Unannotated abundance share by treatment |

Output: `results/limma_*_n88.csv`, `metadata/wf22_design.tsv`

