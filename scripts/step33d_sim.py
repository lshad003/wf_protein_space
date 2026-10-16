# Step 33d: simulated host positive control for step33c.
# 5,000 fragments from the wood frog genome (catalog/host/GCA_028564925.1_aRanSyl1.merge_genomic.fna.gz, uppercased).
# Lengths: drawn with replacement from the lengths (Tiara len=) of the 630,389 set-A contigs (Tiara eukarya + organelle
# >= 1 kb, 88 samples); positions uniform over the genome (sequence chosen in proportion to length); a fragment that
# runs off its sequence or contains any N is redrawn (count reported). Then 1% substitutions: each position
# independently with p = 0.01, replaced by one of the three other bases. random.Random(33) throughout.
# Output (refuses to overwrite): results/step33d_sim/sim5000.fa, headers sim_<i>|<seqid>|<start0>|<len>|<n_subs>
import os, re, gzip, bisect, random
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OD = R + "step33d_sim"; OUT = OD + "/sim5000.fa"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i]
lens = []
for s in stems:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2)
            if f[1] in ("eukarya", "organelle"): lens.append(int(re.search(r"len=(\d+)", f[0]).group(1)))
if len(lens) != 630389: raise SystemExit("ABORT: %d set-A lengths, expected 630,389" % len(lens))
ids, seqs, cur = [], [], []
with gzip.open(W + "/catalog/host/GCA_028564925.1_aRanSyl1.merge_genomic.fna.gz", "rt") as fh:
    for l in fh:
        if l[0] == ">":
            if cur: seqs.append("".join(cur).upper()); cur = []
            ids.append(l[1:].split()[0])
        else: cur.append(l.strip())
    seqs.append("".join(cur).upper())
cum, t = [], 0
for s in seqs: t += len(s); cum.append(t)
print("genome: %d sequences, %d bp" % (len(seqs), t), flush=True)
if t != 5150352190: raise SystemExit("ABORT: genome length differs from NCBI")
rng = random.Random(33); ALT = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
os.makedirs(OD, exist_ok=True); redraw = 0; L = []
with open(OUT, "w") as o:
    for i in range(5000):
        n = rng.choice(lens)
        while True:
            g = rng.randrange(t); j = bisect.bisect_right(cum, g); st = g - (cum[j - 1] if j else 0)
            frag = seqs[j][st:st + n]
            if len(frag) == n and "N" not in frag and set(frag) <= set("ACGT"): break
            redraw += 1
        fr = list(frag); k = 0
        for p in range(n):
            if rng.random() < 0.01: fr[p] = rng.choice(ALT[fr[p]]); k += 1
        o.write(">sim_%d|%s|%d|%d|%d\n%s\n" % (i, ids[j], st, n, k, "".join(fr))); L.append((n, k))
L.sort()
print("fragments 5000 | redrawn (off end or N) %d | length min %d median %d max %d | substitution rate %.4f"
      % (redraw, L[0][0], L[2500][0], L[-1][0], sum(k for _, k in L) / sum(n for n, _ in L)))
print("set-A lengths: median %d, >= 3 kb %.2f%% | simulated >= 3 kb %.2f%%" % (sorted(lens)[len(lens) // 2],
      100.0 * sum(x >= 3000 for x in lens) / len(lens), 100.0 * sum(n >= 3000 for n, _ in L) / 5000))
print("wrote", OUT)
