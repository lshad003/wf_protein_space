# Step 33d v2 summary: simulated host positive control rerun with minimap2 -c (supersedes the step33d host calls, whose
# identity used seed-chain matches). Same 5,000 fragments (results/step33d_sim/sim5000.fa, step33d_sim.py: frog genome,
# set-A length distribution, 1% substitutions, seed 33). Host rule as step33c2: merged query coverage >= 0.50 and
# identity >= 0.95 (sensitivity 0.90). Identity: sum(base-level matches) / sum(block length) from the CIGAR alignment
# (main); 1 - de (de:f, gap-compressed per-base divergence from the -c alignment; minimap2 writes de instead of dv with -c), weighted by query span as second estimate. Every row must carry cg:Z.
# GATE: exits 4 if fewer than 90% of fragments are host at >= 95% (CIGAR identity).
# Output (refuses to overwrite): results/step33d_sim_v2/sim5000_host.tsv
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
D = R + "step33d_sim_v2/"; OUT = D + "sim5000_host.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
src = {}
for l in open(R + "step33d_sim/sim5000.fa"):
    if l[0] == ">": q = l[1:].strip(); f = q.split("|"); src[q] = (f[1], int(f[2]), int(f[3]))
if len(src) != 5000: raise SystemExit("ABORT: %d simulated fragments" % len(src))
iv, mt, bl, best = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
dew, sp = collections.Counter(), collections.Counter()
for l in open(D + "sim5000.host.paf"):
    f = l.split("\t"); q = f[0]
    if q not in src: raise SystemExit("ABORT: unknown query " + q)
    if "\tcg:Z:" not in l: raise SystemExit("ABORT: PAF row without cg:Z")
    dv = float(re.search(r"\tde:f:([\d.eE+-]+)", l).group(1)); a = int(f[3]) - int(f[2])
    iv[q].append((int(f[2]), int(f[3]))); mt[q] += int(f[9]); bl[q] += int(f[10]); dew[q] += dv * a; sp[q] += a
    if q not in best or int(f[10]) > best[q][0]: best[q] = (int(f[10]), f[5], int(f[7]), int(f[8]))
def cov(v):
    v.sort(); c, cs, ce = 0, -1, -1
    for a, b in v:
        if a > ce: c += max(0, ce - cs); cs, ce = a, b
        else: ce = max(ce, b)
    return c + ce - cs
rows = []
for q, (sid, st, n) in src.items():
    m = q in iv; c = cov(iv[q]) / n if m else 0.0; ic = mt[q] / bl[q] if m else 0.0; idv = 1 - dew[q] / sp[q] if m else 0.0
    b = best.get(q); loc = bool(b and b[1] == sid and b[2] < st + n and b[3] > st)
    rows.append((q, n, c, ic, idv, c >= .5 and ic >= .95, c >= .5 and ic >= .90, c >= .5 and idv >= .95, c >= .5 and idv >= .90, m, loc))
with open(OUT, "w") as o:
    o.write("fragment\tlength\tcov\tid_cigar\tid_1mde\thost95\thost90\thost95_1mde\thost90_1mde\tany_alignment\tat_source\n")
    for r in rows: o.write("%s\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\n" % r)
print("simulated frog fragments (1% substitutions) mapped with -c; host rule as step33c2")
print("   %-8s %6s %16s %16s %16s %16s %16s %16s" % ("bin", "n", "CIGAR host95", "CIGAR host90", "1-de host95", "1-de host90", "any alignment", "at source"))
for lab, lo, hi in (("1-3 kb", 0, 3000), (">= 3 kb", 3000, 10 ** 12), ("all", 0, 10 ** 12)):
    s = [r for r in rows if lo <= r[1] < hi]; n = len(s)
    print("   %-8s %6d %s" % (lab, n, " ".join("%7d (%5.2f%%)" % (k, 100.0 * k / n) for k in [sum(r[i] for r in s) for i in (5, 6, 7, 8, 9, 10)])))
for i, lab in ((3, "CIGAR"), (4, "1-de")):
    x = sorted(r[i] for r in rows if r[9]); print("   %s identity, fragments with an alignment: median %.4f, 5th percentile %.4f" % (lab, x[len(x) // 2], x[len(x) // 20]))
fail = [r for r in rows if not r[5]]
print("   not CIGAR host95: %d | no alignment %d, coverage < 0.5 %d, identity < 0.95 with coverage >= 0.5 %d"
      % (len(fail), sum(1 for r in fail if not r[9]), sum(1 for r in fail if r[9] and r[2] < .5), sum(1 for r in fail if r[2] >= .5 and r[3] < .95)))
print("wrote", OUT)
h = sum(r[5] for r in rows)
if h < 0.9 * 5000: print("GATE FAILED: %d of 5000 (%.2f%%) host at >= 95%%, below 90%%. Stopping." % (h, h / 50)); raise SystemExit(4)
print("GATE: %d of 5000 (%.2f%%) host at >= 95%%, >= 90%%, passed" % (h, h / 50))
