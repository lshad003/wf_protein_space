# Step 33d summary: host-detection positive controls for step33c.
# 1. Simulated fragments (step33d_sim.py) mapped with the step33c settings (aRanSyl1.mmi, minimap2 2.30 -x asm10
#    --secondary=no); host rule as step33c_summary.py: merged query coverage >= 0.50 and identity (sum matches / sum
#    block length) >= 0.95, sensitivity >= 0.90. Also: fraction whose longest alignment lands on the source sequence
#    overlapping the source interval. By length bin (1-3 kb, >= 3 kb) and overall.
# 2. Mitochondrion: aRanSyl1 (GCA_028564925.1) has no sequence labelled mitochondrial (13 chromosomes + HiC_scaffold_*),
#    so organelle contigs cannot be tested against a frog mitogenome. Reported instead: Tiara-organelle contigs
#    (first stage, 88 samples; second-stage split) and their step33c host calls (results/step33c_contig_host.tsv).
# Output (refuses to overwrite): results/step33d_sim/sim5000_host.tsv
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
D = R + "step33d_sim/"; OUT = D + "sim5000_host.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
src = {}
for l in open(D + "sim5000.fa"):
    if l[0] == ">": q = l[1:].strip(); f = q.split("|"); src[q] = (f[1], int(f[2]), int(f[3]))
if len(src) != 5000: raise SystemExit("ABORT: %d simulated fragments" % len(src))
iv, mt, bl, best = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
for l in open(D + "sim5000.host.paf"):
    f = l.split("\t"); q = f[0]
    if q not in src: raise SystemExit("ABORT: unknown query " + q)
    iv[q].append((int(f[2]), int(f[3]))); mt[q] += int(f[9]); bl[q] += int(f[10])
    if q not in best or int(f[10]) > best[q][0]: best[q] = (int(f[10]), f[5], int(f[7]), int(f[8]))
def cov(v):
    v.sort(); c, cs, ce = 0, -1, -1
    for a, b in v:
        if a > ce: c += max(0, ce - cs); cs, ce = a, b
        else: ce = max(ce, b)
    return c + ce - cs
rows = []
for q, (sid, st, n) in src.items():
    c = cov(iv[q]) / n if q in iv else 0.0; idt = mt[q] / bl[q] if q in iv else 0.0
    b = best.get(q); loc = bool(b and b[1] == sid and b[2] < st + n and b[3] > st)
    rows.append((q, n, c, idt, c >= .5 and idt >= .95, c >= .5 and idt >= .90, q in iv, loc))
with open(OUT, "w") as o:
    o.write("fragment\tlength\tcov\tidentity\thost95\thost90\tany_alignment\tat_source\n")
    for r in rows: o.write("%s\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\n" % r)
print("1. simulated frog fragments (1% substitutions), host rule as step33c")
print("   %-8s %6s %16s %16s %16s %16s" % ("bin", "n", "host95", "host90", "any alignment", "at source locus"))
for lab, lo, hi in (("1-3 kb", 0, 3000), (">= 3 kb", 3000, 10 ** 12), ("all", 0, 10 ** 12)):
    s = [r for r in rows if lo <= r[1] < hi]; n = len(s)
    print("   %-8s %6d %s" % (lab, n, " ".join("%7d (%5.2f%%)" % (k, 100.0 * k / n) for k in
          (sum(r[4] for r in s), sum(r[5] for r in s), sum(r[6] for r in s), sum(r[7] for r in s)))))
ids = sorted(r[3] for r in rows if r[6])
print("   identity of fragments with an alignment: median %.4f, 5th percentile %.4f" % (ids[len(ids) // 2], ids[len(ids) // 20]))
fail = [r for r in rows if not r[4]]
print("   not host95: %d | of these no alignment %d, coverage < 0.5 %d, identity < 0.95 with coverage >= 0.5 %d"
      % (len(fail), sum(1 for r in fail if not r[6]), sum(1 for r in fail if r[6] and r[2] < .5), sum(1 for r in fail if r[2] >= .5 and r[3] < .95)))
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i]
org = {}
for s in stems:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 3)
            if f[1] == "organelle": org[(s, f[0].split()[0])] = f[2]
H = {}
for i, l in enumerate(open(R + "step33c_contig_host.tsv")):
    if i:
        f = l.rstrip("\n").split("\t")
        if (f[0], f[1]) in org: H[(f[0], f[1])] = (f[5] != "0", f[7] == "1", f[8] == "1")
if len(H) != len(org): raise SystemExit("ABORT: %d organelle contigs missing from step33c_contig_host.tsv" % (len(org) - len(H)))
print("\n2. aRanSyl1 has no sequence labelled mitochondrial; frog-mitogenome test not possible with this assembly.")
print("   Tiara-organelle contigs (88 samples) and step33c host calls (vs the nuclear assembly):")
for sub in sorted(set(org.values())) + ["all"]:
    ks = [k for k in org if sub == "all" or org[k] == sub]
    print("   %-14s n=%6d | any alignment %d | host95 %d | host90 %d" % (sub, len(ks), sum(H[k][0] for k in ks), sum(H[k][1] for k in ks), sum(H[k][2] for k in ks)))
print("wrote", OUT)
