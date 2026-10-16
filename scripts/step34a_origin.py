# Step 34a origin: one origin label per supported 95% rep (results/fourway_classes_n88.tsv, 2,067,011, n88 catalog),
# from the rep's own contig (protein -> contig as step29d).
# Fly / frog calls (minimap2 2.30 -x asm10 --secondary=no -c; coverage = merged query interval / contig length >= 0.50
# and identity = sum(matches) / sum(block length) from the CIGAR alignment >= 0.95):
#   contigs in results/step33c_lists: frog from results/step33c_contig_host_v2.tsv (host95), fly from
#   results/step33g_contig_fly.tsv (fly95); all other rep contigs: computed here from results/step34a_paf/<stem>.
#   {frog,fly}.paf (step34a_map_array.sh, same rule). Every rep contig must have a call; aborts otherwise.
# Origin (first match): fly | frog | Tiara-eukaryotic other (eukarya + organelle) | Tiara-prokaryotic (bacteria +
#   archaea + prokarya) | Tiara-unknown | <1 kb unclassified (contig absent from Tiara -m 1000). Tiara as step32b
#   (UHM586.41010 Tiara output used only for reps on its contigs). Length bin: ge3kb / 1to3kb (Tiara len=) / lt1kb.
# Candidate microbial sets: strict = Tiara-prokaryotic (contigs >= 1 kb, not fly or frog); broad = strict +
#   Tiara-unknown + <1 kb unclassified (all not fly or frog).
# Outputs (refuse to overwrite): results/step34a_rep_origin.tsv, results/step34a_origin_by_class.tsv
import os, re, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = [R + "step34a_rep_origin.tsv", R + "step34a_origin_by_class.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
X = "UHM586.41010"
GRP = {"bacteria": "Tiara-prokaryotic", "archaea": "Tiara-prokaryotic", "prokarya": "Tiara-prokaryotic",
       "eukarya": "Tiara-eukaryotic other", "organelle": "Tiara-eukaryotic other", "unknown": "Tiara-unknown"}
stems = [l.strip() for l in open(R + "step34a_lists/stems.txt") if l.strip()]
if len(stems) != 89 or stems[-1] != X: raise SystemExit("ABORT: stems.txt is not the 88 + " + X)
frog, fly = {}, {}
for i, l in enumerate(open(R + "step33c_contig_host_v2.tsv")):
    if i: f = l.rstrip("\n").split("\t"); frog[(f[0], f[1])] = f[8] == "1"
for i, l in enumerate(open(R + "step33g_contig_fly.tsv")):
    if i: f = l.rstrip("\n").split("\t"); fly[(f[0], f[1])] = f[6] == "1"
if set(frog) != set(fly): raise SystemExit("ABORT: step33c2 and step33g contig sets differ")
def calls(p, s):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    for l in open(p):
        f = l.split("\t"); k = (s, f[0])
        if "\tcg:Z:" not in l: raise SystemExit("ABORT: PAF row without cg:Z in " + p)
        ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
    out = set()
    for k, v in iv.items():
        v.sort(); c, cs, ce = 0, -1, -1
        for a, b in v:
            if a > ce: c += max(0, ce - cs); cs, ce = a, b
            else: ce = max(ce, b)
        if (c + ce - cs) / ql[k] >= 0.5 and mt[k] / bl[k] >= 0.95: out.add(k)
    return out
nnew = 0
for s in stems:
    new = [(s, l.strip()) for l in open(R + "step34a_lists/%s.txt" % s) if l.strip()]
    if not new: continue
    fr = calls(R + "step34a_paf/%s.frog.paf" % s, s); fl = calls(R + "step34a_paf/%s.fly.paf" % s, s)
    for k in new:
        if k in frog: raise SystemExit("ABORT: contig mapped twice " + str(k))
        frog[k] = k in fr; fly[k] = k in fl; nnew += 1
print("contig calls: %d from step33c2/step33g, %d from step34a" % (len(frog) - nnew, nnew), flush=True)
tia = {}
for s in stems:
    with open(R + "tiara/%s.tiara.tsv" % s) as fh:
        next(fh)
        for l in fh:
            f = l.split("\t", 2); ln = int(re.search(r"len=(\d+)", f[0]).group(1))
            tia[(s, f[0].split()[0])] = (GRP[f[1]], "ge3kb" if ln >= 3000 else "1to3kb")
ORI = ["fly", "frog", "Tiara-eukaryotic other", "Tiara-prokaryotic", "Tiara-unknown", "<1 kb unclassified"]
CL = ["K", "KWP", "GU", "EU"]; BINS = ["ge3kb", "1to3kb", "lt1kb"]
tab = collections.Counter(); both = 0; cont = collections.defaultdict(set)
with open(R + "fourway_classes_n88.tsv") as fh, open(OUT[0], "w") as o:
    next(fh); o.write("rep\tclass4\torigin\tlen_bin\tfly95\tfrog95\n")
    for l in fh:
        r, c = l.split("\t", 2)[:2]; s, rest = r.split("__", 1); k = (s, rest[:rest.rfind("_")])
        if k not in fly: raise SystemExit("ABORT: no fly/frog call for the contig of " + r)
        t = tia.get(k); b = t[1] if t else "lt1kb"
        ori = "fly" if fly[k] else "frog" if frog[k] else t[0] if t else "<1 kb unclassified"
        both += fly[k] and frog[k]; tab[(ori, c, b)] += 1; cont[ori].add(k)
        o.write("%s\t%s\t%s\t%s\t%d\t%d\n" % (r, c, ori, b, fly[k], frog[k]))
N = sum(tab.values())
if N != 2067011: raise SystemExit("ABORT: %d reps labelled" % N)
rows = [("len_bin", "origin") + tuple("n_" + c for c in CL) + ("n_all",) + tuple("pct_of_" + c for c in CL) + ("pct_of_all",)]
for b in BINS + ["all"]:
    bb = BINS if b == "all" else [b]
    col = {c: sum(tab[(o_, c, x)] for o_ in ORI for x in bb) for c in CL}; tot = sum(col.values())
    print("\n== %s (reps: %d) | counts and %% of each class" % (b, tot))
    print("  %-24s " % "origin" + " ".join("%18s" % c for c in CL + ["all"]))
    for o_ in ORI:
        v = {c: sum(tab[(o_, c, x)] for x in bb) for c in CL}; a = sum(v.values())
        rows.append((b, o_) + tuple(v[c] for c in CL) + (a,) + tuple("%.2f" % (100.0 * v[c] / col[c]) if col[c] else "NA" for c in CL) + ("%.2f" % (100.0 * a / tot) if tot else "NA",))
        print("  %-24s " % o_ + " ".join("%18s" % ("%d (%.2f%%)" % (x, 100.0 * x / y if y else 0)) for x, y in [(v[c], col[c]) for c in CL] + [(a, tot)]))
with open(OUT[1], "w") as o:
    for r_ in rows: o.write("\t".join(map(str, r_)) + "\n")
print("\nreps on contigs called both fly and frog (labelled fly): %d" % both)
def cand(oris):
    v = {c: sum(tab[(o_, c, x)] for o_ in oris for x in BINS) for c in CL}; a = sum(v.values())
    return a, v, len(set().union(*(cont[o_] for o_ in oris)))
for lab, oris in (("strict (Tiara-prokaryotic >= 1 kb, not fly/frog)", ["Tiara-prokaryotic"]),
                  ("broad (strict + Tiara-unknown + <1 kb unclassified, not fly/frog)", ["Tiara-prokaryotic", "Tiara-unknown", "<1 kb unclassified"])):
    a, v, nc = cand(oris)
    print("%s: %d reps (%.2f%% of 2,067,011) on %d contigs | %s" % (lab, a, 100.0 * a / N, nc, " ".join("%s %d" % (c, v[c]) for c in CL)))
for p in OUT: print("wrote", p)
