### Step 21. Gene neighbourhood of unknown families

For the 8,464 unknown families with three or more members, annotated domains on
adjacent genes were counted. A domain beside at least half of a family's members
is taken as a conserved neighbour, a hypothesis about function rather than an
assignment.

**Result.** 3,614 families (42.7%) have a conserved annotated neighbour. The
strongest cases sit beside oxidoreduction and respiration domains (GSDH, Rieske,
COX1, adh_short, Aldo_ket_red), transporters (MFS_1, ATP_bind_1), regulators
(HTH_1, TetR-like, PhoU), surface proteins (Fimbrial, AsmA) and peptidases.
Members of these families are 6.29% of supported unannotated genes.

| File | Purpose |
|---|---|
| `scripts/step21b_neighbours_n88.py` | Neighbour domains of unknown families |

Output: `results/unk_family_neighbours_n88.tsv`

