#!/bin/bash
#SBATCH -J step34j_dbs
#SBATCH -p epyc
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -t 3-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b2_foldseek_dbs.%j.log
# Step 34j B2: download the Foldseek databases Alphafold/UniProt50 (AFDB50) and ESMAtlas30 with `foldseek databases`
# (Foldseek GPU build 463739e0..., as Steps 28/29e/34d) into /bigdata/stajichlab/lshad003/software/foldseek_dbs.
# Names: AFDB50, ESMAtlas30. Versions and dates go to results/step34j/foldseek_db_versions.tsv.
# Refuses to overwrite an existing database of either name. PDB and CATH50 are not touched.
F=/bigdata/stajichlab/lshad003/software/foldseek_gpu/foldseek/bin/foldseek
DB=/bigdata/stajichlab/lshad003/software/foldseek_dbs
W=/bigdata/stajichlab/lshad003/wf_protein_space
V=$W/results/step34j/foldseek_db_versions.tsv
for n in AFDB50 ESMAtlas30; do
  if [ -e $DB/$n ] || [ -e $DB/$n.dbtype ]; then echo "REFUSING: $DB/$n exists"; exit 1; fi
done
if [ -e $V ]; then echo "REFUSING: $V exists"; exit 1; fi
hostname; date; $F version
echo -e "db\tfoldseek_name\tversion_file\tdownload_started\tdownload_finished\tn_entries\texit" > $V
for pair in "AFDB50:Alphafold/UniProt50" "ESMAtlas30:ESMAtlas30"; do
  IFS=: read n src <<< "$pair"
  echo "== $n ($src)"; s=$(date -Iseconds)
  $F databases $src $DB/$n $DB/tmp_step34j_$n --threads $SLURM_CPUS_PER_TASK
  rc=$?
  e=$(date -Iseconds)
  ver=$(cat $DB/$n.version 2>/dev/null | tr '\n\t' '  ')
  ne=$(wc -l < $DB/$n.index 2>/dev/null)
  echo -e "$n\t$src\t${ver:-NA}\t$s\t$e\t${ne:-NA}\t$rc" >> $V
  echo "$n exit $rc entries ${ne:-NA} version ${ver:-NA}"; ls -la $DB | grep " $n" | head -30
  rm -rf $DB/tmp_step34j_$n
done
date
