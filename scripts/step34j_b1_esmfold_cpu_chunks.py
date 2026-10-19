# Step 34j B1: re-chunk results/step34j/esmfold_input.faa (8,759 sequences) into N CPU chunks, length-sorted round
# robin, for scripts/step34j_b1_esmfold_cpu.sh. The GPU array (job 30049857) was cancelled before any task started:
# scheduler gave no start time (queued gpu AlphaFold job 30048074 estimated 2026-10-31). Refuses to overwrite.
import os, sys
D = "/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j/"; O = D + "esmfold_chunks_cpu/"
N = int(sys.argv[1])
if os.path.exists(O): raise SystemExit("REFUSING: " + O + " exists")
seq, h = {}, None
for l in open(D + "esmfold_input.faa"):
    if l[0] == ">": h = l[1:].strip(); seq[h] = ""
    else: seq[h] += l.strip()
ids = sorted(seq, key=lambda k: (-len(seq[k]), k)); os.makedirs(O)
for c in range(N):
    with open(O + "C.%d.faa" % c, "w") as o:
        for k in sorted(ids[c::N], key=lambda k: len(seq[k])): o.write(">%s\n%s\n" % (k, seq[k]))
print("sequences %d in %d chunks; residues per chunk %s" % (len(ids), N, sorted(sum(len(seq[k]) for k in ids[c::N]) for c in range(N))[::max(1, N // 6)]))
