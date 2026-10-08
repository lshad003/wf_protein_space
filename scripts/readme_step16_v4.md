### Step 16. Abundance

Reads were mapped with bwa-mem2 2.3 to the 2,069,453 supported representatives
(nucleotide coding sequences). The matrix used in all analyses is restricted to
the 2,067,011 clusters supported at 88 metagenomes.

**Result.** Mean mapped fraction is 0.9587 in WF22, 0.9634 in WF23 and 0.9667 in
WF24, so the catalog represents all three years equally.

| File | Purpose |
|---|---|
| `scripts/step16h_supported_cds.py` | Mapping reference |
| `scripts/step16j_index_reference.sh` | Reference index |
| `scripts/step16q_map_array.sh` | Read mapping and counting |
| `scripts/step16u_matrix_n88.py` | Count matrix at 88 metagenomes |

Output: `results/count_matrix_primary_n88.tsv.gz`, `results/mapping_qc_88.tsv`

