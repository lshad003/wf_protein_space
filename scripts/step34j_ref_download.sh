#!/bin/bash
# Step 34j: reference files for structural-hit names (curl only). Each into its own new directory; refuses to overwrite.
# PDB entry titles: wwPDB entries.idx (IDCODE, HEADER, ACCESSION DATE, COMPOUND, SOURCE, AUTHOR LIST, RESOLUTION,
# EXPERIMENT TYPE). CATH: cath-names.txt (node names) and cath-domain-list.txt (domain -> C.A.T.H), latest release.
W=/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j
for d in $W/ref_pdb $W/ref_cath; do if [ -e $d ]; then echo "REFUSING: $d exists"; exit 1; fi; done
mkdir -p $W/ref_pdb $W/ref_cath
curl -sSfL -o $W/ref_pdb/entries.idx https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx
curl -sSfL -o $W/ref_cath/cath-names.txt http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-names.txt
curl -sSfL -o $W/ref_cath/cath-domain-list.txt http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt
{ date -Iseconds; for f in $W/ref_pdb/entries.idx $W/ref_cath/cath-names.txt $W/ref_cath/cath-domain-list.txt; do echo "$f $(wc -l < $f) lines md5 $(md5sum < $f | cut -c1-32)"; done; grep -m2 -i "release\|version" $W/ref_cath/cath-names.txt; } | tee $W/ref_versions.txt
