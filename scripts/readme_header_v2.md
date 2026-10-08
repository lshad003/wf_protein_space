# Wood frog gut protein space: analysis record

This repository documents the construction of a non-redundant protein catalog
from wood frog (*Lithobates sylvaticus*) faecal metagenomes collected over three
years, and the analyses that characterize it: size and saturation, cohort
contributions, technical drivers of richness, and functional annotation.

## Sample set

88 faecal metagenomes from 59 post-metamorphic froglets. Faeces were collected
weekly and pooled by month; each metagenome is one monthly pool from one animal
(stems named `<animal>.<timepoint>`).

| cohort | metagenomes | animals | design |
|---|---|---|---|
| WF22 | 44 | 15 | three monthly pools per animal (one animal with two) |
| WF23 | 9 | 9 | one pool per animal |
| WF24 | 35 | 35 | one pool per animal |

**WF22.** Three egg masses; two *Basidiobolus* strains and control; three monthly
timepoints. Six of 50 metagenomes were excluded: one control with *Basidiobolus*
reads, one on quality review, and four animals that died during the experiment.

**WF23.** Egg mass 3; 5 control and 4 UHM520.7734.

**WF24.** 32 animals from egg mass 4 and 3 from egg mass 2, across seven treatment
codes. Treated animals received at least two inoculations within the month. Tank
is confounded with egg mass.

**Exclusions.** UHM586.41010, a WF24 control with elevated *Basidiobolus* reads,
is excluded from the catalog, applying the same criterion as WF22 sample
UHM56.10839. Clustering was run on 89 metagenomes; UHM586.41010 members were then
removed, the support rule (at least three members from at least two metagenomes)
was reapplied, and original representatives were retained. UHM590.41012 remains in
the catalog but, with UHM586.41010, is excluded from treatment tests, which use 34
WF24 animals (4 controls, 5 per strain). Old and new values of every affected
number are in `results/n88_number_table.tsv`.

| File | Purpose |
|---|---|
| `scripts/step25a_filter_n88.py` | Supported set, count matrix and sample table at 88 metagenomes |
| `scripts/step25b_membership_n88.py` | Cluster tables at 88 metagenomes |
| `scripts/step25c_number_table.py` | Old and new values of all reported numbers |

Output: `results/supported_reps_95_n88.txt`, `results/count_matrix_primary_n88.tsv.gz`,
`metadata/samples_88_treatment.tsv`, `results/clusters_{95,90,50}_n88.tsv`,
`results/n88_number_table.tsv`

