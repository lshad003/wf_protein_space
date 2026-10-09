#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 1 -c 2 --mem 8gb --time=24:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step26b.%j.log
# Step 26b: download NMPFams2 (Zenodo 10.5281/zenodo.17225887, "Quadrupling the
# protein family space with global metagenomics") with curl, verify md5.
# Download and inspection only; no search.
D=/bigdata/stajichlab/lshad003/wf_protein_space/catalog/nmpfams2
mkdir -p $D
cd $D
curl -s https://zenodo.org/api/records/17225887 > zenodo_record.json
while read f m; do
  if [ -s $f ] && [ "$(md5sum $f | cut -d' ' -f1)" = "$m" ]; then echo "already ok $f"; continue; fi
  echo "$(date) downloading $f"
  curl -sS -L --retry 5 -C - -o $f https://zenodo.org/api/records/17225887/files/$f/content
  got=$(md5sum $f | cut -d' ' -f1)
  if [ "$got" = "$m" ]; then echo "md5 OK $f"; else echo "MD5 MISMATCH $f got $got want $m"; fi
done << 'LIST'
NMPFAMSDB2_MODELS_SCORES.txt cde30f38f161004e57672b28fbc0d315
foldseek_results.tar.gz 254b7bf55518eada4b657278ebbff57d
structures.tar.gz 93a8288cd6fde9b290ccbcfff09380e5
metag_clusters25_names.tsv.bz2 677b976daebebacc289667d9d1fd48b6
iso_clusters25_names.tsv.bz2 d6ae4851bdb76667b233bb4a3bb2dc19
LIST
ls -la $D
echo "DONE step26b"
