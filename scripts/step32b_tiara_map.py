# Step 32b: genes and families on Tiara contig classes (step32a, results/tiara/<stem>.tiara.tsv, tiara 1.0.3 -m 1000).
# 0. All 88 stems of metadata/samples_88_treatment.tsv have a non-empty Tiara output and their array tasks
#    (jobs 29979396, 29979411) are COMPLETED 0:0 in sacct; aborts otherwise.
#    Excluded control UHM586.41010 (step32a2, job 29979778, same settings): used ONLY to label supported 95% reps
#    whose own contig is in UHM586.41010 (step 2, and the NMPFams2 queries in step 4, which are supported reps).
#    Family members (step 3) must come from the 88; a member from UHM586.41010 aborts.
# Join: stem + first word of Tiara's sequence_id (k141_X); contig length from its len= field. Protein IDs are
#   <stem>__k141_X_<gene>; contig parsed as step29d (split on the first "__", drop the last "_<gene>").
#   Contigs absent from Tiara are < 1 kb (checked on UHM102.10840: Tiara holds exactly the 22,719 contigs >= 1000 bp).
# Groups: prok = bacteria + archaea + prokarya; euk = eukarya + organelle; unknown; lt1kb = not classified (< 1 kb).
#   Length bins: ge3kb (main) and 1to3kb (sensitivity), reported separately.
# 2. Supported 95% reps with class4 K/KWP/GU/EU from results/fourway_classes_n88.tsv (step13i_fourway_n88.py;
#    K and KWP as step12f, GU/EU split of U by the ClusteredNR criterion).
# 3. The 15,580 recurrent EU-only 30% families (results/unk_direct30_eu_recurrence.tsv, n_animals >= 10 and
#    n_cohorts >= 2, as step28g; checked against results/prostt5_main_eu15580.faa). Members from
#    results/unk_direct30_clusters.tsv. Majority within a bin: the group (prok/euk/unknown) holding >= 50% of the
#    members on contigs in that bin, if exactly one group does; "mixed" otherwise; "unclassifiable" if no member is
#    in that bin.
# 4. NMPFams2: EU queries of results/unk_filtered_aplakidou.faa; matched = any row in results/nmpfams2_metag_rep_hits.tsv
#    or results/nmpfams2_iso_rep_hits.tsv (as step26d4, published 99.33% unmatched).
# 5. Positive control: Cse2 family UHM636.41030__k141_332084_5 should be majority prok at ge3kb.
# Outputs (refuse to overwrite): results/step32b_class_tiara.tsv, results/step32b_family_tiara.tsv,
#   results/step32b_nmpfams2_tiara.tsv
import os, re, subprocess, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
OUT = [R + "step32b_class_tiara.tsv", R + "step32b_family_tiara.tsv", R + "step32b_nmpfams2_tiara.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
CSE2 = "UHM636.41030__k141_332084_5"; X = "UHM586.41010"
GRP = {"bacteria": "prok", "archaea": "prok", "prokarya": "prok", "eukarya": "euk", "organelle": "euk", "unknown": "unknown"}
CATS = ["prok_ge3kb", "euk_ge3kb", "unknown_ge3kb", "prok_1to3kb", "euk_1to3kb", "unknown_1to3kb", "lt1kb"]
# 0
stems = [l.split("\t")[0] for i, l in enumerate(open(W + "/metadata/samples_88_treatment.tsv")) if i]
if len(stems) != 88: raise SystemExit("ABORT: %d stems, expected 88" % len(stems))
sa = subprocess.run(["sacct", "-j", "29979396,29979411,29979778", "-X", "-n", "-P", "-o", "JobID,State,ExitCode"],
                    capture_output=True, text=True).stdout.split()
task = {x.split("|")[0]: x.split("|")[1:] for x in sa}
idx = {"29979396_1": 1}; idx.update({"29979411_%d" % i: i for i in range(2, 89)})
idx["29979778"] = X
bad = [j for j in idx if task.get(j) != ["COMPLETED", "0:0"]]
if bad: raise SystemExit("ABORT: array tasks not COMPLETED 0:0: %s" % bad[:5])
tia, nc, raw, ncx = {}, collections.Counter(), collections.Counter(), collections.Counter()
for s in stems + [X]:
    p = R + "tiara/%s.tiara.tsv" % s
    if not os.path.isfile(p) or os.path.getsize(p) == 0: raise SystemExit("ABORT: missing or empty " + p)
    with open(p) as fh:
        if not next(fh).startswith("sequence_id\tclass_fst_stage"): raise SystemExit("ABORT: unexpected header in " + p)
        for l in fh:
            f = l.rstrip("\n").split("\t")
            c, ln = f[0].split()[0], int(re.search(r"len=(\d+)", f[0]).group(1))
            if f[1] not in GRP: raise SystemExit("ABORT: unknown Tiara class %s in %s" % (f[1], p))
            if ln < 1000: raise SystemExit("ABORT: contig < 1 kb in Tiara output " + p)
            tia[(s, c)] = GRP[f[1]] + ("_ge3kb" if ln >= 3000 else "_1to3kb")
            if s == X: ncx[tia[(s, c)]] += 1; continue
            raw[f[1]] += 1; nc[tia[(s, c)]] += 1
print("0. Tiara outputs: 88 of 88 present, non-empty; array tasks 88 of 88 COMPLETED 0:0; %s (job 29979778) COMPLETED 0:0" % X)
print("   contigs >= 1 kb classified: %d | first stage %s" % (len(tia), dict(raw.most_common())))
print("   by group and bin: %s" % {k: nc[k] for k in CATS[:6]})
print("   %s (rep labelling only): %d contigs >= 1 kb, %s" % (X, sum(ncx.values()), {k: ncx[k] for k in CATS[:6]}))
S88 = set(stems); nx = collections.Counter()
def cat(g, rep=False):
    s, rest = g.split("__", 1)
    if s == X and rep: nx["rep"] += 1
    elif s not in S88: raise SystemExit("ABORT: %s %s from a stem outside the 88" % ("rep" if rep else "member", g))
    return tia.get((s, rest[:rest.rfind("_")]), "lt1kb")
def pct(a, b): return "%.2f" % (100.0 * a / b) if b else "NA"
# 2
print("\n2. classes from results/fourway_classes_n88.tsv (written by scripts/step13i_fourway_n88.py)")
cc = collections.defaultdict(collections.Counter)
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]; cc[c][cat(r, True)] += 1
CL = ["K", "KWP", "GU", "EU"]
if set(cc) != set(CL): raise SystemExit("ABORT: classes %s" % sorted(cc))
if sum(sum(v.values()) for v in cc.values()) != 2067011: raise SystemExit("ABORT: supported reps are not 2,067,011")
hdr = ["class", "n_reps"] + ["n_" + k for k in CATS] + ["pct_all_" + k for k in CATS] \
      + ["pct_of_ge3kb_" + g for g in ("prok", "euk", "unknown")] + ["pct_of_1to3kb_" + g for g in ("prok", "euk", "unknown")]
rows = [hdr]
for c in CL + ["all"]:
    v = sum(cc.values(), collections.Counter()) if c == "all" else cc[c]; n = sum(v.values())
    b3 = sum(v[g + "_ge3kb"] for g in ("prok", "euk", "unknown")); b1 = sum(v[g + "_1to3kb"] for g in ("prok", "euk", "unknown"))
    rows.append([c, n] + [v[k] for k in CATS] + [pct(v[k], n) for k in CATS]
                + [pct(v[g + "_ge3kb"], b3) for g in ("prok", "euk", "unknown")] + [pct(v[g + "_1to3kb"], b1) for g in ("prok", "euk", "unknown")])
print("   supported reps whose own contig is in %s (labelled from its Tiara output): %d" % (X, nx["rep"]))
rx = collections.Counter(l.split("\t", 2)[1] for i, l in enumerate(open(R + "fourway_classes_n88.tsv")) if i and l.startswith(X + "__"))
print("   by class: %s" % dict(rx)); nx.clear()
with open(OUT[0], "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
# 3
fam = []
with open(R + "unk_direct30_eu_recurrence.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ia, ic = h.index("n_animals"), h.index("n_cohorts")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if int(f[ia]) >= 10 and int(f[ic]) >= 2: fam.append(f[0])
fa = {l[1:].split()[0] for l in open(R + "prostt5_main_eu15580.faa") if l[0] == ">"}
if len(fam) != 15580 or set(fam) != fa: raise SystemExit("ABORT: family set is not the 15,580 of prostt5_main_eu15580.faa")
want = set(fam); mc = {f: collections.Counter() for f in fam}
for l in open(R + "unk_direct30_clusters.tsv"):
    r, m = l.rstrip("\n").split("\t")
    if r in want: mc[r][cat(m)] += 1
def major(v, b):
    t = {g: v[g + "_" + b] for g in ("prok", "euk", "unknown")}; n = sum(t.values())
    if n == 0: return "unclassifiable"
    w = [g for g, k in t.items() if 2 * k >= n]
    return w[0] if len(w) == 1 else "mixed"
fs = {"ge3kb": collections.Counter(), "1to3kb": collections.Counter()}
with open(OUT[1], "w") as o:
    o.write("family_rep\tn_members\t" + "\t".join("n_" + k for k in CATS)
            + "\tmajority_ge3kb\tmajority_1to3kb\tfrac_members_ge1kb\tfrac_members_ge3kb\n")
    for f in fam:
        v = mc[f]; n = sum(v.values())
        if n == 0: raise SystemExit("ABORT: family %s has no members" % f)
        m3, m1 = major(v, "ge3kb"), major(v, "1to3kb"); fs["ge3kb"][m3] += 1; fs["1to3kb"][m1] += 1
        g3 = sum(v[k] for k in CATS if k.endswith("ge3kb")); g1 = n - v["lt1kb"]
        o.write("%s\t%d\t%s\t%s\t%s\t%.4f\t%.4f\n" % (f, n, "\t".join(str(v[k]) for k in CATS), m3, m1, g1 / n, g3 / n))
tm = sum(mc.values(), collections.Counter()); N = sum(tm.values())
print("\n3. 15,580 recurrent EU-only families: %d member proteins; %s" % (N, {k: "%d (%s%%)" % (tm[k], pct(tm[k], N)) for k in CATS}))
for b in ("ge3kb", "1to3kb"):
    print("   majority at %s: %s" % (b, {k: "%d (%s%%)" % (fs[b][k], pct(fs[b][k], 15580)) for k in ("prok", "euk", "unknown", "mixed", "unclassifiable")}))
# 4
q = [l[1:].split()[0] for l in open(R + "unk_filtered_aplakidou.faa") if l[0] == ">"]
c4 = {}
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]; c4[r] = c
hit = {l.split("\t", 1)[0] for f in ("metag", "iso") for l in open(R + "nmpfams2_%s_rep_hits.tsv" % f)}
eu = [x for x in q if c4.get(x) == "EU"]
if len(eu) != 266848: raise SystemExit("ABORT: EU queries %d, expected 266,848" % len(eu))
nq, nu = collections.Counter(), collections.Counter()
for x in eu:
    k = cat(x, True); nq[k] += 1; nu[k] += x not in hit
grp = [("all EU (published 99.33%)", CATS), ("prok_ge3kb", ["prok_ge3kb"]), ("prok_1to3kb", ["prok_1to3kb"]),
       ("euk (ge3kb + 1to3kb)", ["euk_ge3kb", "euk_1to3kb"]), ("euk_ge3kb", ["euk_ge3kb"]), ("euk_1to3kb", ["euk_1to3kb"]),
       ("unknown (ge3kb + 1to3kb)", ["unknown_ge3kb", "unknown_1to3kb"]), ("lt1kb (not classified)", ["lt1kb"]),
       ("unclassified = unknown + lt1kb", ["unknown_ge3kb", "unknown_1to3kb", "lt1kb"])]
rows = [("contig_group", "n_EU_queries", "n_matched_metag_or_iso", "n_unmatched", "pct_unmatched")]
for lab, ks in grp:
    a = sum(nq[k] for k in ks); b = sum(nu[k] for k in ks); rows.append((lab, a, a - b, b, pct(b, a)))
if rows[1][3] != 265061: raise SystemExit("ABORT: total unmatched %d, published 265,061" % rows[1][3])
with open(OUT[2], "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("\n4. NMPFams2 (metag or iso rep hits), EU queries; %d of them on %s contigs (labelled from its Tiara output)" % (nx["rep"], X))
for r in rows: print("\t".join(map(str, r)))
# 5
v = mc[CSE2]
print("\n5. Cse2 family %s: members %s | majority ge3kb %s, 1to3kb %s | rep's own contig: %s | %s"
      % (CSE2, {k: v[k] for k in CATS if v[k]}, major(v, "ge3kb"), major(v, "1to3kb"), cat(CSE2, True),
         "PASS" if major(v, "ge3kb") == "prok" else "FAIL"))
for p in OUT: print("wrote", p)
