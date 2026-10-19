# Step 34k final: Figure 1 inputs with member presence restricted to strict-eligible contigs.
# Part A: contig labels from results/step34k/contig_pre.tsv (scripts/step34k_prep.py); pending contigs resolved with
#   fly/frog calls from results/step34g/paf, results/step34j/ctl_paf, results/step34k/paf (minimap2 -c asm10; merged
#   query cover >= 0.50 with sum(matches)/sum(block) >= 0.90 for fly, >= 0.95 for frog) or results/step33c_contig_host_v2.tsv
#   (host95) + results/step33g_contig_fly.tsv (fly90). Contig absent from a PAF = no hit.
# Part B: presence of each strict_final rep (results/step34a_rep_origin_v2_final.tsv) in each of the 88 metagenomes
#   (metadata/samples_88_treatment.tsv) from its 95% cluster members (results/clusters_95_n88.tsv): strict_members =
#   members on strict-eligible contigs only (primary); all_members = any member (sensitivity, = scripts/fig_prep_accum_strict.py).
#   Accumulation as fig_prep_accum_strict.py: 200 random orders, numpy default_rng(1) re-seeded per rule so both rules use
#   the same orders; sets pooled (88), WF22, WF23, WF24.
# Part C: cohort combination of every strict_final rep under each rule, counted by class4.
# Part D: per-metagenome reads in results/count_matrix_primary_n88.tsv.gz (rows = supported reps; counts = primary
#   alignment records, each mate counted, scripts/step16q_map_array.sh) for all matrix rows and strict_final rows.
# Outputs (refuse to overwrite): results/step34k/contig_labels.tsv, results/fig_prep_accum_v2.tsv,
#   results/fig_prep_cohort_sharing.tsv, results/fig_prep_read_retention.tsv, results/step34k_summary.tsv
import os, gzip, collections
import numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step34k/"
OUT = [D + "contig_labels.tsv", R + "fig_prep_accum_v2.tsv", R + "fig_prep_cohort_sharing.tsv", R + "fig_prep_read_retention.tsv", R + "step34k_summary.tsv"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
SUM = []
def add(sec, metric, n, src):
    SUM.append((sec, metric, str(n), src)); print("%-10s %-78s %s" % (sec, metric, n), flush=True)
def paf_calls(p):
    iv, mt, bl, ql = collections.defaultdict(list), collections.Counter(), collections.Counter(), {}
    for l in open(p):
        f = l.split("\t"); k = f[0]
        ql[k] = int(f[1]); iv[k].append((int(f[2]), int(f[3]))); mt[k] += int(f[9]); bl[k] += int(f[10])
    out = {}
    for k, v in iv.items():
        v.sort(); c, cs, ce = 0, -1, -1
        for a, b in v:
            if a > ce: c += max(0, ce - cs); cs, ce = a, b
            else: ce = max(ce, b)
        out[k] = ((c + ce - cs) / ql[k], mt[k] / bl[k])
    return out
# ---- Part A
rows = [l.rstrip("\n").split("\t") for l in open(D + "contig_pre.tsv")]; H = rows[0]; rows = rows[1:]
pend_src = collections.defaultdict(set)
for r in rows:
    if r[8] == "1": pend_src[r[7]].add(r[0])
hit = {}   # (stem, contig) -> True if fly or frog
PD = {"step34g_paf": (R + "step34g/", "stems_to_map.txt", "paf"), "step34j_paf": (R + "step34j/", "ctl_stems_to_map.txt", "ctl_paf"),
      "step34k_paf": (D, "stems_to_map.txt", "paf")}
for tag, (dd, sl, pdir) in PD.items():
    if not pend_src[tag]: continue
    for s in (l.strip() for l in open(dd + sl) if l.strip()):
        if s not in pend_src[tag]: continue
        for kind, idt in (("fly", 0.90), ("frog", 0.95)):
            p = dd + "%s/%s.%s.paf" % (pdir, s, kind)
            if not os.path.exists(p): raise SystemExit("ABORT: missing " + p)
            for c, (cv, i) in paf_calls(p).items():
                if cv >= 0.5 and i >= idt: hit[(s, c, tag)] = True
t33 = {}
if pend_src["step33c_step33g_tables"]:
    for i, l in enumerate(open(R + "step33c_contig_host_v2.tsv")):
        if i: f = l.rstrip("\n").split("\t"); t33[(f[0], f[1])] = f[8] == "1"
    for i, l in enumerate(open(R + "step33g_contig_fly.tsv")):
        if i:
            f = l.rstrip("\n").split("\t"); k = (f[0], f[1])
            if k in t33: t33[k] = t33[k] or f[7] == "1"
elig = set(); A = collections.Counter(); AM = collections.Counter(); nme = 0
with open(OUT[0], "w") as o:
    o.write("stem\tcontig\tn_members\tlabel_source\torigin\tfly90\tlca_domain_v2\tcall_source\tfly_or_frog\tstrict_eligible\n")
    for s, c, n, ls, ori, f90, dom, cs, pend in rows:
        ff = "NA"
        if ls == "step34a_rep_label":
            e = ori == "Tiara-prokaryotic" and f90 == "0" and dom != "Eukaryota"; why = "rep_label_" + ("eligible" if e else "not")
        elif pend == "1":
            if cs == "step33c_step33g_tables":
                if (s, c) not in t33: raise SystemExit("ABORT: no step33 call for %s %s" % (s, c))
                h_ = t33[(s, c)]
            else: h_ = hit.get((s, c, cs), False)
            ff = str(int(h_)); e = not h_; why = "computed_" + cs + ("_eligible" if e else "_fly_or_frog")
        else:
            e = False; why = "computed_not_eligible_" + ("lt1kb" if ori == "<1 kb unclassified" else ("euk_lca" if ori == "Tiara-prokaryotic" else ori.replace(" ", "_")))
        if e: elig.add(s + "\t" + c); nme += int(n)
        A[why] += 1; AM[why] += int(n)
        o.write("\t".join([s, c, n, ls, ori, f90, dom, cs, ff, str(int(e))]) + "\n")
src = "results/step34k/contig_labels.tsv"
add("A", "member contigs", len(rows), src); add("A", "member proteins", sum(AM.values()), src)
add("A", "strict-eligible contigs", len(elig), src); add("A", "strict-eligible member proteins", nme, src)
for k in sorted(A): add("A", "contigs | " + k, A[k], src); add("A", "member proteins | " + k, AM[k], src)
del rows, hit, t33
# ---- Part B
meta = [l.rstrip("\n").split("\t") for l in open(W + "/metadata/samples_88_treatment.tsv")][1:]
stems = [r[0] for r in meta]; coh = np.array([r[1] for r in meta]); si = {s: i for i, s in enumerate(stems)}
if len(stems) != 88: raise SystemExit("ABORT: %d samples" % len(stems))
strict = {}; cls = []
with open(R + "step34a_rep_origin_v2_final.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t"); ir, ic, isf = h.index("rep"), h.index("class4"), h.index("strict_final")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[isf] == "1": strict[f[ir]] = len(strict); cls.append(f[ic])
NS = len(strict)
if NS != 1359273: raise SystemExit("ABORT: strict_final %d" % NS)
pa = [[] for _ in stems]; ps = [[] for _ in stems]; nmem = nel = nout = 0
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t"); j = strict.get(r)
        if j is None: continue
        s, rest = m.split("__", 1); k = si.get(s); nmem += 1
        if k is None: nout += 1; continue
        pa[k].append(j)
        if s + "\t" + rest.rsplit("_", 1)[0] in elig: ps[k].append(j); nel += 1
add("B", "strict_final member proteins in clusters_95_n88", nmem, "results/clusters_95_n88.tsv")
add("B", "  from stems outside the 88", nout, "results/clusters_95_n88.tsv")
add("B", "  on strict-eligible contigs", nel, "results/clusters_95_n88.tsv + results/step34k/contig_labels.tsv")
P = {"all_members": [np.unique(np.array(x, dtype=np.int32)) for x in pa], "strict_members": [np.unique(np.array(x, dtype=np.int32)) for x in ps]}
del pa, ps
CO = ["WF22", "WF23", "WF24"]; NP = 200
bits = {}
for rule in ("strict_members", "all_members"):
    b = np.zeros(NS, np.int8)
    for k, x in enumerate(P[rule]): b[x] |= 1 << CO.index(coh[k])
    bits[rule] = b
    add("B", "%s: strict_final reps present in >= 1 of 88" % rule, int((b > 0).sum()), "results/fig_prep_accum_v2.tsv")
with open(OUT[1], "w") as out:
    out.write("rule\tset\tn_metagenomes\tmean\tsd\tq025\tq975\n")
    for rule in ("strict_members", "all_members"):
        rng = np.random.default_rng(1)
        for name, idx in [("pooled", np.arange(88))] + [(c, np.where(coh == c)[0]) for c in CO]:
            cur = np.zeros((NP, len(idx)), dtype=np.int64)
            for b in range(NP):
                sn = np.zeros(NS, bool); tot = 0
                for t, k in enumerate(rng.permutation(idx)):
                    x = P[rule][k]; tot += int((~sn[x]).sum()); sn[x] = True; cur[b, t] = tot
            for t in range(len(idx)):
                v = cur[:, t]; out.write("%s\t%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\n" % (rule, name, t + 1, v.mean(), v.std(), np.percentile(v, 2.5), np.percentile(v, 97.5)))
            add("B", "%s %s: final accumulation (%d metagenomes)" % (rule, name, len(idx)), int(cur[0, -1]), "results/fig_prep_accum_v2.tsv")
            if len(idx) > 1:
                add("B", "%s %s: mean gain from last metagenome" % (rule, name), "%.1f" % (cur[:, -1] - cur[:, -2]).mean(), "results/fig_prep_accum_v2.tsv")
if os.path.exists(R + "fig_prep_accum_strict.tsv"):
    old = {(f[0], f[1]): f[2] for f in (l.rstrip("\n").split("\t") for l in open(R + "fig_prep_accum_strict.tsv")) if f[0] != "set"}
    new = {(f[1], f[2]): f[3] for f in (l.rstrip("\n").split("\t") for l in open(OUT[1])) if f[0] == "all_members"}
    add("B", "all_members rows differing in mean from results/fig_prep_accum_strict.tsv", sum(1 for k in old if new.get(k) != old[k]), "results/fig_prep_accum_strict.tsv")
# ---- Part C
LAB = {1: "WF22", 2: "WF23", 4: "WF24", 3: "WF22+WF23", 5: "WF22+WF24", 6: "WF23+WF24", 7: "all3", 0: "none"}
ORD = [1, 2, 4, 3, 5, 6, 7, 0]; CL = ["K", "KWP", "GU", "EU"]
ca = np.array(cls)
with open(OUT[2], "w") as o:
    o.write("rule\tcombination\tclass\tn_reps\n")
    for rule in ("strict_members", "all_members"):
        for b in ORD:
            sel = bits[rule] == b
            for c in CL + ["all"]:
                n = int(sel.sum()) if c == "all" else int((sel & (ca == c)).sum())
                o.write("%s\t%s\t%s\t%d\n" % (rule, LAB[b], c, n))
            add("C", "%s %s (all classes)" % (rule, LAB[b]), int(sel.sum()), "results/fig_prep_cohort_sharing.tsv")
# ---- Part D
with gzip.open(R + "count_matrix_primary_n88.tsv.gz", "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t")[1:]
    if sorted(h) != sorted(stems): raise SystemExit("ABORT: matrix columns are not the 88")
    cat = np.zeros(len(h), np.int64); st = np.zeros(len(h), np.int64); nr = ns = 0
    for l in fh:
        f = l.rstrip("\n").split("\t"); v = np.array(f[1:], dtype=np.int64); cat += v; nr += 1
        if f[0] in strict: st += v; ns += 1
add("D", "matrix rows (supported reps)", nr, "results/count_matrix_primary_n88.tsv.gz")
add("D", "matrix rows that are strict_final reps", ns, "results/count_matrix_primary_n88.tsv.gz")
add("D", "strict_final reps not in the matrix", NS - ns, "results/count_matrix_primary_n88.tsv.gz")
hc = {s: i for i, s in enumerate(h)}
with open(OUT[3], "w") as o:
    o.write("stem\tcohort\tcatalog_reads\tstrict_reads\tpct_retained\n")
    for k, s in enumerate(stems):
        i = hc[s]; o.write("%s\t%s\t%d\t%d\t%.2f\n" % (s, coh[k], cat[i], st[i], 100.0 * st[i] / cat[i]))
src = "results/fig_prep_read_retention.tsv"
for name, idx in [("pooled", np.arange(88))] + [(c, np.where(coh == c)[0]) for c in CO]:
    ci = [hc[stems[k]] for k in idx]; pct = 100.0 * st[ci] / cat[ci]
    add("D", "%s: catalog reads (sum)" % name, int(cat[ci].sum()), src); add("D", "%s: strict reads (sum)" % name, int(st[ci].sum()), src)
    add("D", "%s: pooled retention sum(strict)/sum(catalog), %%" % name, "%.2f" % (100.0 * st[ci].sum() / cat[ci].sum()), src)
    add("D", "%s: per-metagenome retention median / min / max, %%" % name, "%.2f / %.2f / %.2f" % (np.median(pct), pct.min(), pct.max()), src)
with open(OUT[4], "w") as o:
    o.write("section\tmetric\tvalue\tsource_file\n")
    for r in SUM: o.write("\t".join(r) + "\n")
print("wrote", ", ".join(OUT))
