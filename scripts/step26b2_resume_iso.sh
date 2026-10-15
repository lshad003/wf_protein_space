#!/bin/bash
#SBATCH -p stajichlab -c 1 --mem 2G --time 12:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26b2.log
D=/bigdata/stajichlab/lshad003/wf_protein_space/catalog/nmpfams2
U=$(/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 -c "
import json; r=json.load(open('$D/zenodo_record.json'))
print([f['links']['self'] for f in r['files'] if f['key']=='iso_clusters25_names.tsv.bz2'][0])")
echo "$U"
for i in 1 2 3 4 5; do curl -L -C - -o $D/iso_clusters25_names.tsv.bz2 "$U" && break; sleep 60; done
md5sum $D/iso_clusters25_names.tsv.bz2
echo "want d6ae4851bdb76667b233bb4a3bb2dc19"
