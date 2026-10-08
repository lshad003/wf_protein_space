### Step 1. Gene prediction

Proteins were predicted from assembled contigs with Prodigal in metagenomic mode.

**Result.** 40,550,595 proteins from 89 metagenomes; 40,257,962 after excluding
one WF24 sample.

| File | Purpose |
|---|---|
| `scripts/step1c_prodigal_array.sh` | Gene prediction |
| `scripts/step1d_rename_45.sh` | Headers prefixed with source metagenome |
| `scripts/step1e_protein_counts.sh` | Protein counts per metagenome |

Output: `results/step1_protein_counts_45.tsv`

