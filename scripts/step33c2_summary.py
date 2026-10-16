# Step 33c v2 summary (supersedes step33c_summary.py, whose identity used seed-chain matches from PAFs mapped without
# -c; see step33d): host (wood frog) status of contigs from the step33c2 PAFs (results/step33c_paf_v2, minimap2 -c).
# Every PAF row must carry a cg:Z CIGAR; aborts otherwise. Identity (main) = sum(col 10, base-level matches) /
# sum(col 11, alignment block length incl. gaps); second estimate 1 - de (de:f, gap-compressed per-base divergence from the -c alignment; minimap2 writes de instead of dv with -c), weighted by query span.
# Per contig over all PAF rows (minimap2 --secondary=no; primary and supplementary): coverage = merged query interval
# length / contig length; identity = sum(matches, col 10) / sum(block length, col 11).
# Host = coverage >= 0.50 and identity >= 0.95 (main); identity >= 0.90 sensitivity. Same rule for Basidiobolus.
# Gate: if the set C negative control (Tiara-prokaryotic >= 3 kb, 2%) has host fraction > 1% at the main rule, prints
# A/B/C and exits 3 before the remaining reports.
# Tiara group and length bin per contig as step32b (prok/euk/unknown x ge3kb/1to3kb; lt1kb = absent from Tiara);
# UHM586.41010 Tiara output only for the EU reps whose own contig is there.
# Outputs (refuse to overwrite): results/step33c_contig_host_v2.tsv, results/step33c_eu_host_v2.tsv,
#   results/step33c_family_host_v2.tsv
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; L = R + "step33c_lists/"; P = R + "step33c_paf_v2/"
OUT = [R + "step33c_contig_host_v2.tsv", R + "step33c_eu_host_v2.tsv", R + "step33c_family_host_v2.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
CSE2 = "UHM636.41030__k141_332084_5"; X = "UHM586.41010"
GRP = {"bacteria": "prok", "archaea": "prok", "prokarya": "prok", "eukarya": "euk", "organelle": "euk", "unknown": "unknown"}
CATS = ["prok_ge3kb", "euk_ge3kb", "unknown_ge3kb", "prok_1to3kb", "euk_1to3kb", "unknown_1to3kb", "lt1kb"]
stems = [l.strip() for l in open(L + "stems.txt") if l.strip()]
if len(stems) != 89 or stems[-1] != X: raise SystemExit("ABORT: stems.txt is not the 88 + " + X)
sets, tia = {}, {}
for s in stems:
    for l in open(L + s + ".tsv"):
        c, v = l.rstrip("\n").split("\t"); sets[(s, c)] = v
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2); ln = int(re.search(r"len=(\d+)", f[0]).group(1))
            tia[(s, f[0].split()[0])] = GRP[f[1]] + ("_ge3kb" if ln >= 3000 else "_1to3kb")
def stats(kind):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    dew, sp = collections.Counter(), collections.Counter()
    for s in stems:
        p = P + "%s.%s.paf" % (s, kind)
        if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
        for l in open(p):
            f = l.split("\t"); k = (s, f[0])
            if k not in sets: raise SystemExit("ABORT: unlisted contig in " + p)
            if "\tcg:Z:" not in l: raise SystemExit("ABORT: PAF row without cg:Z (not mapped with -c) in " + p)
            dv = re.search(r"\tde:f:([\d.eE+-]+)", l)
            if not dv: raise SystemExit("ABORT: PAF row without de:f in " + p)
            ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
            dew[k] += float(dv.group(1)) * (int(f[3]) - int(f[2])); sp[k] += int(f[3]) - int(f[2])
    out = {}
    for k, v in iv.items():
        v.sort(); cov, cs, ce = 0, -1, -1
        for a, b in v:
            if a > ce: cov += max(0, ce - cs); cs, ce = a, b
            else: ce = max(ce, b)
        cov += ce - cs
        out[k] = (cov / ql[k], mt[k] / bl[k], ql[k], 1 - dew[k] / sp[k])
    return out
H = stats("host")
def call(st, k, idt, e=1):
    v = st.get(k); return bool(v and v[0] >= 0.5 and v[e] >= idt)
def frac(keys, st, idt, e=1):
    n = len(keys); h = sum(call(st, k, idt, e) for k in keys); return n, h, "%.2f" % (100.0 * h / n) if n else "NA"
by = {x: [k for k, v in sets.items() if x in v] for x in "ABCF"}
print("host = coverage >= 50% and identity >= 95% (sensitivity: identity >= 90%); identity from the CIGAR alignment (main) and 1 - de")
by["B_excl_" + X] = [k for k in by["B"] if k[0] != X]
for x, lab in (("A", "Tiara-eukaryotic >= 1 kb"), ("B", "contigs with an EU supported rep"), ("C", "negative control: 2% Tiara-prokaryotic >= 3 kb"), ("F", "contigs with a member of the 15,580 families"), ("B_excl_" + X, "B excluding contigs on " + X)):
    print("  %s %-48s n=%d | CIGAR host95 %d (%s%%) host90 %d (%s%%) | 1-de host95 %d (%s%%) host90 %d (%s%%) | any alignment %d" % ((x, lab) + frac(by[x], H, .95)[:1] + frac(by[x], H, .95)[1:] + frac(by[x], H, .90)[1:] + frac(by[x], H, .95, 3)[1:] + frac(by[x], H, .90, 3)[1:] + (sum(1 for k in by[x] if k in H),)))
nC, hC, _ = frac(by["C"], H, .95)
if hC > 0.01 * nC: print("GATE FAILED: negative control host fraction above 1%; thresholds too loose. Stopping."); raise SystemExit(3)
print("GATE: negative control host fraction <= 1%, passed")
with open(OUT[0], "w") as o:
    o.write("stem\tcontig\tsets\ttiara\tlength\thost_cov\thost_id_cigar\thost_id_1mde\thost95\thost90\thost95_1mde\thost90_1mde\n")
    for k, v in sorted(sets.items()):
        h = H.get(k); o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (k[0], k[1], v, tia.get(k, "lt1kb"), h[2] if h else "NA",
            "%.4f" % h[0] if h else "0", "%.4f" % h[1] if h else "NA", "%.4f" % h[3] if h else "NA", call(H, k, .95), call(H, k, .90), call(H, k, .95, 3), call(H, k, .90, 3)))
def contig(g):
    s, rest = g.split("__", 1); return s, rest[:rest.rfind("_")]
eu = collections.defaultdict(lambda: [0, 0, 0]); eux = [0, 0, 0]
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]
        if c != "EU": continue
        k = contig(r); t = tia.get(k, "lt1kb"); e = eu[t]; e[0] += 1; e[1] += call(H, k, .95); e[2] += call(H, k, .90)
        if k[0] != X: eux[0] += 1; eux[1] += call(H, k, .95); eux[2] += call(H, k, .90)
rows = [("tiara_group_bin", "n_EU_reps", "n_host95", "pct_host95", "n_host90", "pct_host90")]
for t in CATS + ["all"]:
    e = [sum(eu[x][i] for x in CATS) for i in range(3)] if t == "all" else eu[t]
    rows.append((t, e[0], e[1], "%.2f" % (100.0 * e[1] / e[0]) if e[0] else "NA", e[2], "%.2f" % (100.0 * e[2] / e[0]) if e[0] else "NA"))
with open(OUT[1], "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("\nEU supported reps (430,299, n88 catalog incl. %s) by Tiara group and contig length (CIGAR identity)" % X)
for r in rows: print("\t".join(map(str, r)))
print("EU reps excluding %s: n %d | host95 %d (%.3f%%) | host90 %d (%.3f%%)" % (X, eux[0], eux[1], 100.0 * eux[1] / eux[0], eux[2], 100.0 * eux[2] / eux[0]))
fam = set()
with open(R + "unk_direct30_eu_recurrence.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ia, ic = h.index("n_animals"), h.index("n_cohorts")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if int(f[ia]) >= 10 and int(f[ic]) >= 2: fam.add(f[0])
fm = {f: [0, 0, 0] for f in fam}
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in fm:
        k = contig(m); v = fm[r]; v[0] += 1; v[1] += call(H, k, .95); v[2] += call(H, k, .90)
mh95 = sum(1 for v in fm.values() if 2 * v[1] >= v[0]); mh90 = sum(1 for v in fm.values() if 2 * v[2] >= v[0])
anyh = sum(1 for v in fm.values() if v[1])
with open(OUT[2], "w") as o:
    o.write("family_rep\tn_members\tn_host95\tn_host90\tmajority_host95\tmajority_host90\n")
    for f in sorted(fm):
        v = fm[f]; o.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (f, v[0], v[1], v[2], 2 * v[1] >= v[0], 2 * v[2] >= v[0]))
print("\n15,580 families (majority host = host members >= 50%% of members): host95 %d (%.2f%%), host90 %d (%.2f%%); "
      "families with >= 1 host95 member %d" % (mh95, 100.0 * mh95 / 15580, mh90, 100.0 * mh90 / 15580, anyh))
v = fm[CSE2]
print("Cse2 family %s: %d members, host95 %d, host90 %d (expected 0): %s" % (CSE2, v[0], v[1], v[2], "PASS" if v[1] == 0 and v[2] == 0 else "FAIL"))
B = stats("basid")
nh = [k for k in by["A"] if not call(H, k, .95)]
print("\nBasidiobolus (six genomes), Tiara-eukaryotic >= 1 kb contigs that are not host95: n=%d | CIGAR basid95 %d (%s%%) basid90 %d (%s%%) | 1-de basid95 %d (%s%%) basid90 %d (%s%%)"
      % ((len(nh),) + frac(nh, B, .95)[1:] + frac(nh, B, .90)[1:] + frac(nh, B, .95, 3)[1:] + frac(nh, B, .90, 3)[1:]))
org = {}
for s in stems[:-1]:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 3)
            if f[1] == "organelle": org[(s, f[0].split()[0])] = f[2]
for sub in sorted(set(org.values())) + ["all"]:
    ks = [k for k in org if sub == "all" or org[k] == sub]
    print("Tiara-organelle %-13s n=%d | host95 %d host90 %d | basid95 %d" % ((sub, len(ks)) + (frac(ks, H, .95)[1], frac(ks, H, .90)[1], frac(ks, B, .95)[1])))
for p in OUT: print("wrote", p)
