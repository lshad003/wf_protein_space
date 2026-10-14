import glob, os, re, time, statistics as st
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; P = W + "/results/esmfold_pdb"
lab = {}
with open(W + "/results/esmfold_input_labels.tsv") as fh:
    next(fh)
    for l in fh:
        r, g, n, a = l.rstrip("\n").split("\t"); lab[r] = (g, int(n), int(a))
if len(lab) != 618: raise SystemExit("ABORT: expected 618 labels, got %d" % len(lab))
log = sorted(glob.glob(W + "/logs/step29b_esmfold.*.log"), key=os.path.getmtime)[-1]
rx = re.compile(r"Predicted structure for (\S+) with length (\d+), pLDDT ([\d.]+), pTM ([\d.]+)")
done = {m.group(1): (float(m.group(3)), float(m.group(4))) for m in rx.finditer(open(log).read())}
rows = []
for r, (pl, pt) in done.items():
    f = P + "/" + r + ".pdb"
    if r not in lab or not os.path.exists(f) or time.time() - os.path.getmtime(f) < 120: continue
    X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(f) if l.startswith("ATOM") and l[12:16].strip() == "CA"])
    n = len(X); d2 = np.linalg.norm(X[2:] - X[:-2], axis=1); d3 = np.linalg.norm(X[3:] - X[:-3], axis=1)
    hel = float(np.mean((d3 >= 4.5) & (d3 <= 6.0))); strd = float(np.mean((d2[:-1] >= 6.2) & (d3 >= 9.0)))
    rgr = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean())) / (2.2 * n ** 0.38)
    cls = ("helix_only" if hel >= 0.8 and strd < 0.05 else "folded") if (pl >= 70 and pt >= 0.5) else "low_confidence"
    rows.append((r, lab[r][0], n, lab[r][2], pl, pt, hel, strd, rgr, cls))
print("log %s | completed in log %d of 618 | read %d" % (os.path.basename(log), len(done), len(rows)))
def summ(sub, tag):
    print("== %s" % tag)
    for g in ("HIT", "NOHIT", "NOHIT_LONG"):
        s = [x for x in sub if x[1] == g]; tot = sum(1 for v in lab.values() if v[0] == g)
        if not s: print("  %-10s 0/%d" % (g, tot)); continue
        pc = lambda c: 100.0 * sum(1 for x in s if c(x)) / len(s)
        print("  %-10s %3d/%3d | len %4d | pLDDT %.1f, >=70 %4.1f%% | pTM %.2f, >=0.5 %4.1f%%, >=0.7 %4.1f%% | folded %4.1f%%, helix-only %d" % (
            g, len(s), tot, st.median(x[2] for x in s), st.median(x[4] for x in s), pc(lambda x: x[4] >= 70), st.median(x[5] for x in s),
            pc(lambda x: x[5] >= 0.5), pc(lambda x: x[5] >= 0.7), pc(lambda x: x[9] == "folded"), sum(1 for x in s if x[9] == "helix_only")))
summ(rows, "all finished")
Lmax = max(x[2] for x in rows); summ([x for x in rows if x[2] >= 100], "length-matched, 100-%d aa" % Lmax)
print("top NOHIT/NOHIT_LONG 'folded' by pTM (len, n_animals, pLDDT, pTM, helix, strand, Rg ratio):")
for x in sorted([x for x in rows if x[1] != "HIT" and x[9] == "folded"], key=lambda x: -x[5])[:12]:
    print("  %-34s %-10s %4d %3d %5.1f %.2f %.2f %.2f %.2f" % (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]))
with open(W + "/results/esmfold_peek_v2.tsv", "w") as o:
    o.write("rep\tgroup\tlength\tn_animals\tplddt\tptm\thelix_frac\tstrand_frac\trg_ratio\tclass\n")
    for x in rows: o.write("%s\t%s\t%d\t%d\t%.1f\t%.3f\t%.2f\t%.2f\t%.2f\t%s\n" % x)
