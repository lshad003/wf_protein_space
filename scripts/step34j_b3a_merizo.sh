#!/bin/bash
#SBATCH -J step34j_merizo
#SBATCH -p epyc
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -t 2-00:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/wf_protein_space/logs/step34j_b3a_merizo.%j.log
# Step 34j B3a: domain segmentation with Merizo (software/Merizo commit 41d12fb84e6e8fdb586c2c859d12161dc7bb5bfd, env
# condaenvs/merizo, CPU) of every ESMFold model >= 200 aa (results/step34j/esmfold_input_labels.tsv lengths;
# models in results/step34j/esmfold_pdb). --length_conditional_iterate (iterative mode for >= 512 aa, as Merizo
# advises for long models). Models are copied to results/step34j/merizo/in (Merizo writes next to its input).
# Domain PDBs are written from the Merizo residue ranges (result column) by this script, not by Merizo:
# results/step34j/domains/<model>__D<k>.pdb, residues of domain k (discontinuous segments kept together).
# Table: results/step34j/merizo/merizo_results.tsv. Refuses to overwrite.
W=/bigdata/stajichlab/lshad003/wf_protein_space
D=$W/results/step34j
for p in $D/merizo $D/domains; do if [ -e $p ]; then echo "REFUSING: $p exists"; exit 1; fi; done
mkdir -p $D/merizo/in $D/domains
hostname; date
awk -F'\t' 'NR>1 && $4>=200 {print $1}' $D/esmfold_input_labels.tsv > $D/merizo/ids_ge200.txt
miss=0
while read id; do
  if [ -s $D/esmfold_pdb/$id.pdb ]; then cp $D/esmfold_pdb/$id.pdb $D/merizo/in/; else miss=$((miss+1)); echo "no model: $id"; fi
done < $D/merizo/ids_ge200.txt
echo "models >= 200 aa listed $(wc -l < $D/merizo/ids_ge200.txt), copied $(ls $D/merizo/in | grep -c '\.pdb$'), missing $miss"
cd $D/merizo
/bigdata/stajichlab/lshad003/condaenvs/merizo/bin/python /bigdata/stajichlab/lshad003/software/Merizo/predict.py \
  -d cpu -t $SLURM_CPUS_PER_TASK --length_conditional_iterate --output_headers -i $D/merizo/in/*.pdb > $D/merizo/merizo_results.tsv
rc=$?
echo "merizo exit $rc; rows $(grep -vc '^input' $D/merizo/merizo_results.tsv)"
if [ $rc -ne 0 ]; then echo "ABORT: merizo exit $rc"; exit 1; fi
# 2026-10-19: the writer below fails with KeyError 'filename' (Merizo names the column 'input' and repeats the header
# per model); job 30050461 hit this and the domains were written by scripts/step34j_b3a2_domains.py instead.
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - <<'PYEOF'
import os
D = "/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j/"
rows = [l.rstrip("\n").split("\t") for l in open(D + "merizo/merizo_results.tsv") if l.strip()]
h = rows[0]; print("header:", h); rows = [dict(zip(h, r)) for r in rows[1:]]
nd = 0; nm = 0
for r in rows:
    name = os.path.basename(r["filename"]).replace(".pdb", "")
    res = r["result"].strip()
    if not res or res in ("NA", "-"): continue
    doms = res.split(",")
    lines = [l for l in open(D + "esmfold_pdb/%s.pdb" % name) if l.startswith("ATOM")]
    nm += 1
    for k, d in enumerate(doms, 1):
        keep = set()
        for seg in d.split("_"):
            a, b = seg.split("-"); keep |= set(range(int(a), int(b) + 1))
        with open(D + "domains/%s__D%d.pdb" % (name, k), "w") as o:
            for l in lines:
                if int(l[22:26]) in keep: o.write(l)
            o.write("END\n")
        nd += 1
print("models with domains %d of %d, domain PDBs written %d" % (nm, len(rows), nd))
PYEOF
rc=$?
echo "domain writer exit $rc"; date
# 2026-10-19: exit with the main command's status so afterok dependencies see failures (CLAUDE.md rule)
exit $rc
