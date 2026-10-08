### Step 1. Gene prediction and catalog input

Proteins were predicted from per-sample assembled contigs with Prodigal in
metagenomic mode, with headers carrying the source metagenome. For 31 of the 45
WF23 and WF24 metagenomes, existing predictions made with the same version and
parameters were reused.

**Main result.** The 45 WF23 and WF24 metagenomes contribute 24,151,134 proteins
(WF23 3,642,871; WF24 20,508,263) and the 44 WF22 metagenomes 16,399,461, for a
clustering input of 40,550,595 from 89 metagenomes; removing UHM586.41010
(292,633 proteins) gives 40,257,962 at 88. Gene density is constant across
prediction batches (median 4,580 and 4,602 proteins per Mb, range 4,047 to
4,969), so protein yield per sample tracks assembly size.

| File | Purpose |
|---|---|
| `scripts/step0_check_existing_predictions.sh`, `step0b_check_by_assembly.sh`, `step0c_locate_predict_script.sh`, `step1b_preflight.sh` | Inventory of existing predictions and input contigs |
| `scripts/step1c_prodigal_array.sh` | Prodigal on the 14 remaining metagenomes |
| `scripts/step1d_rename_45.sh` | Headers prefixed with source metagenome |
| `scripts/step1e_protein_counts.sh` | Protein counts against assembly size |

Output: `results/step1_protein_counts_45.tsv`, `results/step0_by_assembly_check.tsv`

