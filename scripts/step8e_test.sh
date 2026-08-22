#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -c 16 --mem 64gb --time=2:00:00
#SBATCH --out /bigdata/stajichlab/lshad003/wf_protein_space/logs/step8e_test.%j.log
WPS=/bigdata/stajichlab/lshad003/wf_protein_space
SRC=/bigdata/operations/pkgadmin/srv/projects/db/eggNOG/5.0.2
LOCAL=/scratch/lshad003_eggnog_db
mkdir -p $LOCAL $WPS/results/eggnog_test
for F in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
  [ -s "$LOCAL/$F" ] || cp "$SRC/$F" "$LOCAL/$F"
done
echo "staged on $(hostname)"; ls -lh $LOCAL

echo; echo "=== TEST A: 2.1.9 full run, chunk 600 ==="
module load eggnog-mapper/2.1.9
emapper.py -i $WPS/catalog/db/LsPS_AA_95_rep__split/LsPS.600 \
  -o A600 --output_dir $WPS/results/eggnog_test --cpu 16 \
  -m diamond --itype proteins --evalue 0.00001 --sensmode very-sensitive \
  --data_dir $LOCAL 2>&1 | tail -8
echo "A rows: $(grep -vc '^#' $WPS/results/eggnog_test/A600.emapper.annotations 2>/dev/null || echo FAIL)"

echo; echo "=== TEST B: 2.1.9 annotate from existing hits, chunk 1 ==="
emapper.py -m no_search --annotate_hits_table $WPS/results/eggnog/LsPS.1.emapper.hits \
  -o B1 --output_dir $WPS/results/eggnog_test --cpu 16 --data_dir $LOCAL --override 2>&1 | tail -8
echo "B rows: $(grep -vc '^#' $WPS/results/eggnog_test/B1.emapper.annotations 2>/dev/null || echo FAIL)"

echo; echo "=== TEST C: 2.1.7 annotate from existing hits, chunk 1 ==="
module unload eggnog-mapper/2.1.9
module load eggnog-mapper/2.1.7
emapper.py -m no_search --annotate_hits_table $WPS/results/eggnog/LsPS.1.emapper.hits \
  -o C1 --output_dir $WPS/results/eggnog_test --cpu 16 --data_dir $LOCAL --override 2>&1 | tail -8
echo "C rows: $(grep -vc '^#' $WPS/results/eggnog_test/C1.emapper.annotations 2>/dev/null || echo FAIL)"
