#!/bin/bash
#SBATCH -p short -w r02 -c 32 --mem 64G -t 2:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step28f.%j.log
# Step 28f: ProstT5 (CPU) + Foldseek search of the calibration sets, as step28e
# (search vs PDB and CATH50, -e 1e-1); convertalis fields as step28e2 (no prob, TM or lddt).
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
M=/bigdata/stajichlab/lshad003/software/foldseek_gpu/weights/prostt5
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
O=$W/results/prostt5_calib; mkdir -p $O/tmp
for q in known740 decoy740; do
  if [ -e $O/${q}_vs_PDB.tsv ]; then echo "ABORT: outputs exist for $q"; exit 1; fi
done
declare -A IN=([known740]=$W/results/prostt5_calib_known_lenmatched.faa [decoy740]=$W/results/prostt5_decoy_eu740rev.faa)
for q in known740 decoy740; do
  date
  $F createdb ${IN[$q]} $O/$q --prostt5-model $M --threads 32
  for t in PDB CATH50; do
    $F search $O/$q $DB/$t $O/${q}_res_$t $O/tmp -e 1e-1 --threads 32
    $F convertalis $O/$q $DB/$t $O/${q}_res_$t $O/${q}_vs_$t.tsv \
      --format-output "query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov"
    echo "$q $t: hit rows $(wc -l < $O/${q}_vs_$t.tsv); queries with a hit $(cut -f1 $O/${q}_vs_$t.tsv | sort -u | wc -l)"
  done
done
date
echo "DONE step28f"
