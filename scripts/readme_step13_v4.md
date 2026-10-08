### Step 13. Genomic and environmental unknowns

The 585,347 unannotated supported representatives were searched against NCBI
ClusteredNR (release 20260128, 470,748,714 sequences) with DIAMOND 2.1.24
--very-sensitive. A hit requires E <= 1e-10 with at least 50% coverage of both
query and subject.

**Result.** 430,299 representatives (20.82% of the supported catalog) have no hit
in ClusteredNR (environmental unknowns). 155,048 (7.50%) match only sequences
annotated as hypothetical (genomic unknowns). Across looser and stricter criteria
the environmental unknown share ranges from 63.92% to 76.63% of unannotated
representatives.

| File | Purpose |
|---|---|
| `scripts/step13d_extract_u.py` | Unannotated representatives to fasta |
| `scripts/step13f_nrclust_search.sh` | DIAMOND search against ClusteredNR |
| `scripts/step13h_gu_eu_split_n88.sh` | Genomic and environmental split under several criteria |
| `scripts/step13i_fourway_n88.py` | Final four-class labels |

Output: `results/step13_nrclust20260128_hits.tsv`, `results/fourway_classes_n88.tsv`

