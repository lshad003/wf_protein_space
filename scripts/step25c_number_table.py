# Step 25c: results/n88_number_table.tsv. For every README number from the
# steps rerun at 88 metagenomes, the old value is parsed from the 89-sample
# log or output and the new value from the _n88 log or output. No value is
# typed in by hand; a value that cannot be parsed is written as NOT_FOUND.
import re, glob, gzip, collections, csv
W = "/bigdata/stajichlab/lshad003/wf_protein_space/"
L = W + "logs/"; R = W + "results/"
def one(pat):
    f = sorted(glob.glob(pat))
    return f[-1] if f else None
def txt(p):
    try: return open(p).read() if p else ""
    except OSError: return ""
def rx(p, pat, g=1, flags=re.M):
    m = re.search(pat, txt(p), flags)
    return m.group(g) if m else "NOT_FOUND"
def sect(p, start, end=None):
    t = txt(p); i = t.find(start)
    if i < 0: return ""
    t = t[i:]
    if end:
        j = t.find(end, len(start))
        if j > 0: t = t[:j]
    return t
def rxs(p, start, pat, end=None, g=1):
    m = re.search(pat, sect(p, start, end), re.M)
    return m.group(g) if m else "NOT_FOUND"
def tsv(p):
    try: return list(csv.reader(open(p), delimiter="\t"))
    except OSError: return []
def row(p, key):
    for r in tsv(p):
        if r and r[0] == key: return r
    return None
def pct(a, b): return "%.2f" % (100.0 * a / b)
rows = []
def add(script, q, old, new, path): rows.append((script, q, str(old), str(new), path))

# ---- step7 / step7b
o7, n7 = one(L + "step7.*.log"), one(L + "step7_n88.*.log")
add("step7_rarefaction.sh", "pooled proteins (catalog input)", rx(o7, r"total proteins: ([\d,]+)"), rx(n7, r"total proteins: ([\d,]+)"), n7)
add("step7_rarefaction.sh", "pooled clusters at full depth", rx(o7, r"40,550,595 proteins ->\s+([\d,]+)"),
    rx(n7, r"([\d,]+) clusters\s+\(\+[\d,]+ over last [\d,]+\)\s*\n\s*\nwrote"), n7)
add("step7_rarefaction.sh", "clusters added 35M to 40M proteins", rx(o7, r"40,000,000 proteins ->.*\(\+([\d,]+) over"), rx(n7, r"40,000,000 proteins ->.*\(\+([\d,]+) over"), n7)
o7b, n7b = one(L + "step7b.*.log"), one(L + "step7b_n88.*.log")
for lab in ("WF22 \\(all 44\\)", "WF23 \\(all 9, run 0426\\)", "WF24 run 0426", "WF24 run 0730"):
    for T in ("3,600,000", "8,000,000"):
        f = lambda p: rxs(p, "richness at %s" % T, lab + r"\s+n=\s*(\d+) pool=\s*([\d,]+)\s+([\d,]+|insufficient)", "===  " if 0 else ("richness at 8" if T == "3,600,000" else None), 3)
        nn = lambda p: rxs(p, "richness at %s" % T, lab + r"\s+n=\s*(\d+)", ("richness at 8" if T == "3,600,000" else None), 1)
        add("step7b_wf24_byrun.sh", "%s clusters at %s proteins (n=%s -> %s)" % (lab.replace("\\", ""), T, nn(o7b), nn(n7b)), f(o7b), f(n7b), n7b)
def num(s): return float(s.replace(",", "")) if s not in ("NOT_FOUND", "insufficient") else None
for p_, tag in ((o7b, "old"), (n7b, "new")): pass
def ratios(p):
    g = lambda lab: num(rxs(p, "richness at 3,600,000", lab + r"\s+n=\s*\d+ pool=\s*[\d,]+\s+([\d,]+)", "richness at 8"))
    a, b, c, d = g("WF22 \\(all 44\\)"), g("WF23 \\(all 9, run 0426\\)"), g("WF24 run 0426"), g("WF24 run 0730")
    if None in (a, b, c, d): return ("NOT_FOUND",) * 3
    return ("%.1f" % (100 * (a - b) / b), "%.1f" % (100 * (d / c - 1)), "%.1f" % (100 * (c / b - 1)))
ro, rn = ratios(o7b), ratios(n7b)
add("step7b_wf24_byrun.sh", "WF22 vs WF23 difference at 3.6M, %", ro[0], rn[0], n7b)
add("step7b_wf24_byrun.sh", "WF24 run 0730 over run 0426 at 3.6M, %", ro[1], rn[1], n7b)
add("step7b_wf24_byrun.sh", "WF24 run 0426 over WF23 at 3.6M, %", ro[2], rn[2], n7b)

# ---- step9
for lab, key in (("all proteins", "all_proteins"), ("representatives", "representatives"),
                 ("singleton", "cluster_size_singleton"), ("size 2", "cluster_size_2"), ("size 3-5", "cluster_size_3-5"),
                 ("size 6-20", "cluster_size_6-20"), ("size 21+", "cluster_size_21+")):
    o, n = row(R + "orf_completeness.tsv", key), row(R + "orf_completeness_n88.tsv", key)
    add("step9_orf_completeness.sh", "complete ORF %% (%s), criterion partial=00" % lab, o[3] if o else "NOT_FOUND", n[3] if n else "NOT_FOUND", R + "orf_completeness_n88.tsv")
    if key in ("all_proteins", "representatives"):
        add("step9_orf_completeness.sh", "n (%s)" % lab, o[1] if o else "NOT_FOUND", n[1] if n else "NOT_FOUND", R + "orf_completeness_n88.tsv")

# ---- step10
for lab, key in (("all representatives", "all_representatives"), ("singleton", "size_singleton"), ("size 21+", "size_21+")):
    o, n = row(R + "annotation_summary.tsv", key), row(R + "annotation_summary_n88.tsv", key)
    add("step10_annotation_summary.sh", "eggNOG annotated %% (%s)" % lab, o[3] if o else "NOT_FOUND", n[3] if n else "NOT_FOUND", R + "annotation_summary_n88.tsv")
o, n = row(R + "annotation_summary.tsv", "all_representatives"), row(R + "annotation_summary_n88.tsv", "all_representatives")
add("step10_annotation_summary.sh", "representatives (n)", o[1] if o else "NOT_FOUND", n[1] if n else "NOT_FOUND", R + "annotation_summary_n88.tsv")
add("step10_annotation_summary.sh", "eggNOG annotated representatives (n)", o[2] if o else "NOT_FOUND", n[2] if n else "NOT_FOUND", R + "annotation_summary_n88.tsv")

# ---- step11
o11, n11 = one(L + "step11.*.log"), one(L + "step11_n88.*.log")
for c, lab in ((1, "complete"), (4, "spans contig")):
    pat = r"^singleton" + r"\s+([\d.]+)%" * 4
    add("step11_interaction_and_supported.sh", "singleton annotation %% among %s ORFs" % lab, rxs(o11, "WITHIN cluster size", pat, None, c), rxs(n11, "WITHIN cluster size", pat, None, c), n11)
for lab in ("all families", ">=3 members and >=2 metagenomes", ">=3 members and >=3 metagenomes"):
    o, n = row(R + "family_dark_fraction.tsv", lab), row(R + "family_dark_fraction_n88.tsv", lab)
    add("step11_interaction_and_supported.sh", "50%% families unannotated %% (%s)" % lab, o[4] if o else "NOT_FOUND", n[4] if n else "NOT_FOUND", R + "family_dark_fraction_n88.tsv")
    if lab.startswith(">=3 members and >=2"):
        add("step11_interaction_and_supported.sh", "supported families (n, >=3 members >=2 metagenomes)", o[1] if o else "NOT_FOUND", n[1] if n else "NOT_FOUND", R + "family_dark_fraction_n88.tsv")
        add("step11_interaction_and_supported.sh", "supported unannotated families (n)", int(o[1]) - int(o[2]) if o else "NOT_FOUND", int(n[1]) - int(n[2]) if n else "NOT_FOUND", R + "family_dark_fraction_n88.tsv")

# ---- step13 taxonomy
o13, n13 = one(L + "step13.*.log"), one(L + "step13_n88.*.log")
add("step13_taxonomy.sh", "annotated representatives (n)", rx(o13, r"annotated representatives \(n=([\d,]+)\)"), rx(n13, r"annotated representatives \(n=([\d,]+)\)"), n13)
for k in ("Bacteria", "Eukaryota", "Archaea", "Viruses"):
    add("step13_taxonomy.sh", "%% %s of annotated" % k, rx(o13, r"^\s+%s\s+[\d,]+\s+([\d.]+)%%" % k), rx(n13, r"^\s+%s\s+[\d,]+\s+([\d.]+)%%" % k), n13)
for k in ("Alphaproteobacteria", "Bacteroidetes", "Actinobacteria", "Gammaproteobacteria", "Betaproteobacteria", "Metazoa", "Fungi"):
    add("step13_taxonomy.sh", "max_annot_lvl %% %s" % k, rx(o13, r"\|%s\s+[\d,]+\s+([\d.]+)%%" % k), rx(n13, r"\|%s\s+[\d,]+\s+([\d.]+)%%" % k), n13)
add("step13_taxonomy.sh", "Fungi representatives (n)", rx(o13, r"\|Fungi\s+([\d,]+)"), rx(n13, r"\|Fungi\s+([\d,]+)"), n13)
for b in ("singleton", "21+"):
    pat = r"^%s\s+[\d.]+%%\s+([\d.]+)%%" % re.escape(b)
    add("step13_taxonomy.sh", "%% Eukaryota among annotated, cluster size %s" % b, rxs(o13, "clade against cluster size", pat), rxs(n13, "clade against cluster size", pat), n13)

# ---- step12d / step12f
u_o, u_n = tsv(R + "union_rep_level.tsv"), tsv(R + "union_rep_level_n88.tsv")
def U(t, scope, col):
    h = [r for r in t if r and r[0] == "scope"]; r = [r for r in t if r and r[0] == scope]
    if not h or not r: return None
    return int(r[0][h[0].index(col)])
for scope in ("all_reps", "supported"):
    for col in ("total", "eggnog", "pfam", "union", "neither"):
        a, b = U(u_o, scope, col), U(u_n, scope, col)
        add("step12d_union_replevel.sh", "%s %s (n)" % (scope, col), a, b, R + "union_rep_level_n88.tsv")
        if col != "total":
            add("step12d_union_replevel.sh", "%s %s %%" % (scope, col), pct(a, U(u_o, scope, "total")) if a is not None else "NOT_FOUND",
                pct(b, U(u_n, scope, "total")) if b is not None else "NOT_FOUND", R + "union_rep_level_n88.tsv")
o12, n12 = L + "step12f.log", L + "step12f_n88.log"
for i, k in enumerate(("K", "KWP", "U")):
    add("step12f_rep_classes.sh", "supported %s (n)" % k, rx(o12, r"sup K/KWP/U: (\d+) (\d+) (\d+)", i + 1), rx(n12, r"sup K/KWP/U: (\d+) (\d+) (\d+)", i + 1), R + "rep_classes_n88.tsv")
add("step12f_rep_classes.sh", "supported DUF/UPF-only K (n)", rx(n12, r"89 supported .* DUF/UPF-only K: (\d+)"), rx(n12, r"^sup DUF/UPF-only K: (\d+)"), n12)
add("step12f_rep_classes.sh", "supported strict K (n)", rx(n12, r"89 supported .* strict K: (\d+)"), rx(n12, r"^sup DUF/UPF-only K: \d+ +strict K: (\d+)"), n12)

# ---- step13i four-way (regenerated intermediate)
o13i, n13i = L + "step13i.log", L + "step13i_n88.log"
for k in ("K", "KWP", "GU", "EU"):
    pat = r"^%s\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)" % k
    for g, lab in ((1, "n"), (2, "% of supported"), (3, "mean CDS bp"), (4, "complete ORF %"), (6, "in all 3 years %"), (7, "% of mapped reads")):
        add("step13i_fourway.py (regenerated)", "%s %s" % (k, lab), rx(o13i, pat, g), rx(n13i, pat, g), R + "fourway_classes_n88.tsv")

# ---- step16v
o16, n16 = L + "step16v.log", L + "step16v_n88.log"
for k in ("K", "KWP", "U"):
    add("step16v_prevalence.py", "%s genes in all samples (89 -> 88)" % k, rxs(o16, "in all 89", r"^%s (\d+)" % k), rxs(n16, "in all 88", r"^%s (\d+)" % k), R + "prevalence_primary_n88.tsv.gz")
def prev(p, top):
    m = re.search(r"^U total (\d+) 1:(\d+)\S+ 2-5:(\d+)\S+ 6-20:(\d+)\S+ 21-60:(\d+)\S+ %s:(\d+)\S+ %s:(\d+)" % (top[0], top[1]), txt(p), re.M)
    if not m: return ("NOT_FOUND",) * 2
    v = [int(x) for x in m.groups()]
    return ("%.2f" % (100.0 * (v[4] + v[5] + v[6]) / v[0]), "%.2f" % (100.0 * (v[1] + v[2]) / v[0]))
po, pn = prev(o16, ("61-88", "89")), prev(n16, ("61-87", "88"))
add("step16v_prevalence.py", "U genes in more than 20 samples %", po[0], pn[0], R + "prevalence_primary_n88.tsv.gz")
add("step16v_prevalence.py", "U genes in 5 or fewer samples %", po[1], pn[1], R + "prevalence_primary_n88.tsv.gz")

# ---- step16w, means over samples from the per-sample tables
def dark(p):
    t = tsv(p)
    if len(t) < 2: return None
    h = t[0]; d = t[1:]
    raw = [float(r[h.index("pct_reads_U")]) for r in d]; rpk = [float(r[h.index("pct_rpk_U")]) for r in d]
    by = collections.defaultdict(list)
    for r in d: by[r[1]].append(float(r[h.index("pct_reads_U")]))
    return {"n samples": len(d), "U % raw reads, mean of samples": "%.2f" % (sum(raw) / len(raw)),
            "U % length-normalized, mean of samples": "%.2f" % (sum(rpk) / len(rpk)),
            "U % raw, WF22 mean": "%.2f" % (sum(by["WF22"]) / len(by["WF22"])),
            "U % raw, WF23 mean": "%.2f" % (sum(by["WF23"]) / len(by["WF23"])),
            "U % raw, WF24 mean": "%.2f" % (sum(by["WF24"]) / len(by["WF24"])),
            "U % raw, per-sample min": "%.2f" % min(raw), "U % raw, per-sample max": "%.2f" % max(raw)}
do, dn = dark(R + "dark_abundance_by_sample.tsv"), dark(R + "dark_abundance_by_sample_n88.tsv")
for k in (do or dn or {}):
    add("step16w_dark_abundance.py", k, do[k] if do else "NOT_FOUND", dn[k] if dn else "NOT_FOUND", R + "dark_abundance_by_sample_n88.tsv")
add("step16w_dark_abundance.py", "U % of supported genes", rx(L + "step16w.log", r"U is ([\d.]+)%"), rx(L + "step16w_n88.log", r"U is ([\d.]+)%"), L + "step16w_n88.log")

# ---- step17g / 17h / 17i
o, n = L + "step17g.log", L + "step17g_n88.log"
add("step17g_limma_wf22.R", "genes in >=35 of 44", rx(o, r"of 44: (\d+)"), rx(n, r"of 44: (\d+)"), n)
add("step17g_limma_wf22.R", "within-animal correlation (refined)", rx(o, r"refined correlation: ([\d.]+)"), rx(n, r"refined correlation: ([\d.]+)"), n)
for blk, s, e in (("B_naive", "B_naive", "A_blocked"), ("A_blocked", "A_blocked", "C_EM3"), ("C_EM3", "C_EM3", None)):
    for cf in ("treatmentSTP1710.7", "treatmentSTP1717.1", "egg_mass2", "egg_mass3", "month2", "month3"):
        pat = r"^\s+%s\s+padj<0.05:\s+(\d+)" % re.escape(cf)
        a, b = rxs(o, "=== " + s, pat, e and "=== " + e), rxs(n, "=== " + s, pat, e and "=== " + e)
        if a == b == "NOT_FOUND": continue
        add("step17g_limma_wf22.R", "%s %s genes padj<0.05" % (blk, cf), a, b, R + "limma_%s_%s_n88.csv" % (blk, cf.replace(".", "_")))
o, n = L + "step17h.log", L + "step17h_n88.log"
add("step17h_limma_wf2324.R", "WF23 genes in >=8 of 9", rxs(o, "=== WF23", r"genes present in >= 8 : (\d+)"), rxs(n, "=== WF23", r"genes present in >= 8 : (\d+)"), n)
add("step17h_limma_wf2324.R", "WF23 UHM520.7734 vs Control genes padj<0.05", rxs(o, "=== WF23", r"treatmentUHM520.7734\s+padj<0.05:\s+(\d+)"),
    rxs(n, "=== WF23", r"treatmentUHM520.7734\s+padj<0.05:\s+(\d+)"), R + "limma_WF23_treatmentUHM520_7734_n88.csv")
o, n = L + "step17i.log", L + "step17i_n88.log"
add("step17i_limma_wf24_n34.R", "WF24 genes in >=30 of 34", rx(o, r"genes present in >= 30: (\d+)"), rx(n, r"genes present in >= 30: (\d+)"), n)
for c in "123567":
    add("step17i_limma_wf24_n34.R", "WF24 code %s vs code 4 genes padj<0.05" % c, rx(o, r"treatment%s\s+padj<0.05:\s+(\d+)" % c), rx(n, r"treatment%s\s+padj<0.05:\s+(\d+)" % c), R + "limma_WF24_n34_treatment%s_n88.csv" % c)

# ---- step22b-e: adonis R2 and p, betadisper p
def ad(p, start, term, end=None):
    m = re.search(r"^%s\s+\d+\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s+([\d.]+)" % term, sect(p, start, end), re.M)
    return m.groups() if m else ("NOT_FOUND", "NOT_FOUND")
def addad(script, label, op, np_, start, term, end, path):
    a, b = ad(op, start, term, end), ad(np_, start, term, end)
    add(script, label + " R2", a[0], b[0], path); add(script, label + " p", a[1], b[1], path)
o, n = L + "step22b.log", L + "step22b_n88.log"
add("step22b_permanova_fix.R", "genes in >=35 of 44", rx(o, r"genes: (\d+)"), rx(n, r"genes: (\d+)"), n)
addad("step22b_permanova_fix.R", "Bray, whole-animal perm, egg_mass", o, n, "=== A.", "egg_mass", "=== B.", n)
addad("step22b_permanova_fix.R", "Bray, whole-animal perm, treatment", o, n, "=== A.", "treatment", "=== B.", n)
addad("step22b_permanova_fix.R", "Bray, within-animal perm, month", o, n, "=== B.", "month", "=== C.", n)
addad("step22b_permanova_fix.R", "Bray, animal centroids, egg_mass", o, n, "=== C.", "egg_mass", "=== D.", n)
addad("step22b_permanova_fix.R", "Bray, animal centroids, treatment", o, n, "=== C.", "treatment", "=== D.", n)
for lab, s, e in (("treatment", "treatment:", "egg mass:"), ("egg mass", "egg mass:", None)):
    pat = r"^Groups\s+\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)"
    add("step22b_permanova_fix.R", "betadisper p, %s, animal centroids" % lab, rxs(o, s, pat, e), rxs(n, s, pat, e), n)
o, n = L + "step22c.log", L + "step22c_n88.log"
add("step22c_jaccard.R", "WF22 genes variable in presence (>=5 reads)", rx(o, r"genes variable in presence: (\d+)  samples"), rx(n, r"genes variable in presence: (\d+)  samples"), n)
add("step22c_jaccard.R", "WF22 mean genes present per sample", rx(o, r"mean genes present per sample: (\d+)"), rx(n, r"mean genes present per sample: (\d+)"), n)
addad("step22c_jaccard.R", "Jaccard, whole-animal perm, egg_mass", o, n, "=== A.", "egg_mass", "=== B.", n)
addad("step22c_jaccard.R", "Jaccard, whole-animal perm, treatment", o, n, "=== A.", "treatment", "=== B.", n)
addad("step22c_jaccard.R", "Jaccard, animal centroids, egg_mass", o, n, "=== B.", "egg_mass", "=== C.", n)
addad("step22c_jaccard.R", "Jaccard, animal centroids, treatment", o, n, "=== B.", "treatment", "=== C.", n)
add("step22c_jaccard.R", "richness by treatment Kruskal-Wallis p", rx(o, r"p-value = ([\d.]+)"), rx(n, r"p-value = ([\d.]+)"), n)
add("step22c_jaccard.R", "WF24 section D animals (old 36, superseded by 22e; new 34)", "36", "34", n)
addad("step22c_jaccard.R", "WF24 section D Jaccard treatment", o, n, "=== D.", "Model", None, n)
o, n = L + "step22d.log", L + "step22d_n88.log"
add("step22d_cazy_subset.R", "carbohydrate-active representatives", rx(o, r"Pfam family: (\d+)"), rx(n, r"Pfam family: (\d+)"), n)
add("step22d_cazy_subset.R", "carbohydrate genes in >=35 of 44", rx(o, r"of 44: (\d+)"), rx(n, r"of 44: (\d+)"), n)
for cf in ("treatmentSTP1710.7", "treatmentSTP1717.1", "egg_mass2", "egg_mass3", "month2", "month3"):
    pat = r"^\s+%s\s+padj<0.05:\s+(\d+)" % re.escape(cf)
    add("step22d_cazy_subset.R", "limma blocked %s genes padj<0.05" % cf, rx(o, pat), rx(n, pat), R + "limma_cazy_%s_n88.csv" % cf.replace(".", "_"))
addad("step22d_cazy_subset.R", "carbohydrate subset centroids, egg_mass", o, n, "=== B.", "egg_mass", "=== C.", n)
addad("step22d_cazy_subset.R", "carbohydrate subset centroids, treatment", o, n, "=== B.", "treatment", "=== C.", n)
add("step22d_cazy_subset.R", "total carbohydrate abundance Kruskal-Wallis p", rx(o, r"p-value = ([\d.]+)"), rx(n, r"p-value = ([\d.]+)"), n)
o, n = L + "step22e.log", L + "step22e_n88.log"
add("step22e_permanova_wf24_n34.R", "WF24 genes in >=30 of 34", rx(o, r"genes: (\d+)"), rx(n, r"genes: (\d+)"), n)
addad("step22e_permanova_wf24_n34.R", "WF24 Bray treatment", o, n, "=== abundance", "Model", "=== presence", n)
add("step22e_permanova_wf24_n34.R", "WF24 genes variable in presence (>=5 reads)", rx(o, r"genes variable in presence: (\d+)"), rx(n, r"genes variable in presence: (\d+)"), n)
addad("step22e_permanova_wf24_n34.R", "WF24 Jaccard treatment", o, n, "=== presence", "Model", None, n)

# ---- step21b
o, n = L + "step21b.log", L + "step21b_n88.log"
add("step21b_neighbours.py", "unknown families >=3 members", rx(o, r"unknown families >=3 members: (\d+)"), rx(n, r"unknown families >=3 members: (\d+)"), R + "unk_family_neighbours_n88.tsv")
add("step21b_neighbours.py", "families with a domain beside >=50% of members", rx(o, r"in >=50% of members: (\d+) of"), rx(n, r"in >=50% of members: (\d+) of"), R + "unk_family_neighbours_n88.tsv")
def share(p):
    a, b = rx(p, r"in >=50% of members: (\d+) of (\d+)"), rx(p, r"in >=50% of members: (\d+) of (\d+)", 2)
    return "NOT_FOUND" if "NOT_FOUND" in (a, b) else "%.1f" % (100.0 * int(a) / int(b))
add("step21b_neighbours.py", "families with a domain beside >=50% of members %", share(o), share(n), R + "unk_family_neighbours_n88.tsv")
def cover(fam_file, cls_file, sup_col):
    c = collections.Counter(); nU = 0
    for ln in open(fam_file): c[ln.split("\t", 1)[0]] += 1
    with open(cls_file) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if p[1] == "1" and p[5] == "U": nU += 1
    return ("%.2f" % (100.0 * sum(v for v in c.values() if v >= 3) / nU), "%.2f" % (100.0 * sum(1 for v in c.values() if v >= 3) / nU))
co, cn = cover(R + "unk_clusters_50.tsv", R + "rep_classes.tsv", 1), cover(R + "unk_clusters_50_n88.tsv", R + "rep_classes_n88.tsv", 1)
add("step21b_neighbours.py", "supported U genes that are members of families >=3, %", co[0], cn[0], R + "unk_clusters_50_n88.tsv")
add("step21b_neighbours.py", "families >=3 divided by supported U genes, % (README 1.4% uses this ratio)", co[1], cn[1], R + "unk_clusters_50_n88.tsv")

# ---- step12e family-tier union
def fu(p, scope):
    r = row(p, scope)
    return (int(r[1]), int(r[3]), int(r[5])) if r else None
for scope in ("all", "supported"):
    a, b = fu(R + "family_union_dark.tsv", scope), fu(R + "family_union_dark_n88.tsv", scope)
    add("step12e_family_union.py", "50%% families, %s (n)" % scope, a[0] if a else "NOT_FOUND", b[0] if b else "NOT_FOUND", R + "family_union_dark_n88.tsv")
    add("step12e_family_union.py", "50%% families, %s, union-dark (n)" % scope, a[2] if a else "NOT_FOUND", b[2] if b else "NOT_FOUND", R + "family_union_dark_n88.tsv")
    add("step12e_family_union.py", "50%% families, %s, union-dark %%" % scope, pct(a[2], a[0]) if a else "NOT_FOUND", pct(b[2], b[0]) if b else "NOT_FOUND", R + "family_union_dark_n88.tsv")
    add("step12e_family_union.py", "50%% families, %s, eggNOG-dark %%" % scope, pct(a[1], a[0]) if a else "NOT_FOUND", pct(b[1], b[0]) if b else "NOT_FOUND", R + "family_union_dark_n88.tsv")

# ---- step14b AntiFam (--cut_ga)
o, n = L + "step14b.log", L + "step14b_n88.log"
pat = r"spurious sequences flagged: (\d+) of (\d+) \(([\d.]+)%\)"
add("step14b_antifam.sh", "AntiFam-flagged supported U reps (n), --cut_ga", rx(o, pat, 1), rx(n, pat, 1), R + "antifam_hits_n88.tblout")
add("step14b_antifam.sh", "supported U reps screened (n)", rx(o, pat, 2), rx(n, pat, 2), R + "supported_U_reps_n88.faa")
add("step14b_antifam.sh", "AntiFam-flagged %", rx(o, pat, 3), rx(n, pat, 3), R + "antifam_hits_n88.tblout")

# ---- step14c unknown family sizes (n88: step14c clustering filtered, not reclustered)
for lab, key, col in (("unknown families total", "families total", 1), ("unknown families >=3 members", ">=3 members", 1),
                      ("unknown families >=100 members", ">=100 members", 1)):
    a, b = row(R + "unk_family_sizes.tsv", key), row(R + "unk_family_sizes_n88.tsv", key)
    add("step14c_novel_families.sh", lab + ", 50% id c 0.8", a[col] if a else "NOT_FOUND", b[col] if b else "NOT_FOUND", R + "unk_family_sizes_n88.tsv")

# ---- step19a core unknowns by years
o, n = L + "step19a.log", L + "step19a_n88.log"
for y in ("3", "1", "0"):
    pat = r"^%s\s+(\d+)\s+(\d+)" % y
    add("step19a_core_unknown.py", "U genes detected in %s years (n)" % y, rx(o, pat, 1), rx(n, pat, 1), R + "core_unknown_by_years_n88.tsv.gz")
    if y == "1":
        add("step19a_core_unknown.py", "U genes detected in 1 year, mean bp", rx(o, pat, 2), rx(n, pat, 2), n)
    if y != "0":
        add("step19a_core_unknown.py", "%s-year U genes, %% of U reads" % y, rx(o, r"%s years: ([\d.]+)%% of U reads" % y), rx(n, r"%s years: ([\d.]+)%% of U reads" % y), n)

# ---- step20c contig taxonomy
o, n = L + "step20c.log", L + "step20c_n88.log"
add("step20c_gene_tax.py", "stems with contig taxonomy (89 -> 88 metagenomes)", rx(o, r"stems with taxonomy: (\d+ of \d+)"),
    rx(n, r"stems with taxonomy among the 88: (\d+ of \d+)"), n)
def ntax(p):
    v = re.findall(r"^\S+\s+n=(\d+)\s", txt(p), re.M)
    return sum(int(x) for x in v) if v else None
def nsup(p): return rx(p, r"supported genes: (\d+)")
a, b = ntax(o), ntax(n)
add("step20c_gene_tax.py", "genes with contig taxonomy (n)", a or "NOT_FOUND", b or "NOT_FOUND", R + "gene_taxonomy_n88.tsv")
add("step20c_gene_tax.py", "genes with contig taxonomy, % of supported", pct(a, int(nsup(o))) if a else "NOT_FOUND", pct(b, int(nsup(n))) if b else "NOT_FOUND", R + "gene_taxonomy_n88.tsv")
for c in ("K", "KWP", "GU", "EU"):
    for d in ("unclassified", "Bacteria", "Eukaryota"):
        pat = r"^%s\s+n=\d+\s+.*?\b%s:([\d.]+)%%" % (c, d)
        add("step20c_gene_tax.py", "%s genes on classified-file contigs, %% %s" % (c, d), rx(o, pat), rx(n, pat), n)

# ---- step2d / step2e search diagnostics (n88: UHM586.41010 queries dropped from the hit files)
o, n = one(L + "step2d.*.log"), L + "step2d_n88.log"
add("step2d_control_and_tiers.sh", "WF22 positive control recovered %, cov-mode 1", rx(o, r"CONTROL .*= ([\d.]+)%"), rx(n, r"CONTROL .*= ([\d.]+)%"), n)
add("step2d_control_and_tiers.sh", "median aa, new proteins mapped at 95% (cov-mode 1)", rx(o, r"^mapped:\s+n=\d+ median_aa=(\d+)"), rx(n, r"^mapped:\s+n=\d+ median_aa=(\d+)"), n)
add("step2d_control_and_tiers.sh", "median aa, new proteins unmapped", rx(o, r"^unmapped:\s+n=\d+ median_aa=(\d+)"), rx(n, r"^unmapped:\s+n=\d+ median_aa=(\d+)"), n)
o, n = one(L + "step2e.*.log"), L + "step2e_n88.log"
for r_, lab in (("control_cov2", "WF22 positive control recovered %, cov-mode 2"), ("new95_cov2", "new-cohort proteins matching at gene level %, cov-mode 2"),
                ("new50_cov2", "new-cohort proteins matching at family level %, cov-mode 2")):
    pat = r"^%s\s+\d+/\s*\d+ =\s+([\d.]+)%%" % r_
    add("step2e_covmode2.sh", lab, rx(o, pat), rx(n, pat), n)
add("step2e_covmode2.sh", "new-cohort queries (n)", rx(o, r"^new95_cov2\s+\d+/\s*(\d+)"), rx(n, r"^new95_cov2\s+\d+/\s*(\d+)"), n)

# ---- step3 / step4 cluster composition
o, n = one(L + "step3.*.log"), one(L + "step3_n88.*.log")
add("step3_cluster_composition.sh", "95% singletons %", rxs(o, "tier 95", r"singletons: \d+ \(([\d.]+)%\)"), rxs(n, "tier 95", r"singletons: \d+ \(([\d.]+)%\)"), n)
add("step3_cluster_composition.sh", "95% all clusters, three cohorts %", rxs(o, "tier 95", r"WF22\+WF23\+WF24\s+\d+\s+\(\s*([\d.]+)%\)"), rxs(n, "tier 95", r"WF22\+WF23\+WF24\s+\d+\s+\(\s*([\d.]+)%\)"), n)
o, n = one(L + "step4.*.log"), one(L + "step4_n88.*.log")
def tsect(text, start):
    i = text.find(start)
    return text[i:] if i >= 0 else ""
for tier, nxt in (("TIER 95", "TIER 50"), ("TIER 50", None)):
    s = ">=3 members AND >=2 samples"
    for lab, pat in (("supported clusters (n)", r">=3 members AND >=2 samples: (\d+) clusters"),
                     ("supported, three cohorts %", r"WF22\+WF23\+WF24\s+\d+\s+\(\s*([\d.]+)%\)"),
                     ("supported, WF22-exclusive %", r"^\s+WF22\s+\d+\s+\(\s*([\d.]+)%\)")):
        a = re.search(pat, tsect(sect(o, tier, nxt), s), re.M); b = re.search(pat, tsect(sect(n, tier, nxt), s), re.M)
        add("step4_nonsingleton.sh", "%s %s" % (tier.replace("TIER ", "") + "% tier", lab), a.group(1) if a else "NOT_FOUND", b.group(1) if b else "NOT_FOUND", n)
    pat = r"top\s+1\.0% of clusters hold\s+([\d.]+)% of proteins"
    add("step4_nonsingleton.sh", tier.replace("TIER ", "") + "% tier top 1% of clusters hold % of proteins", rxs(o, tier, pat, nxt), rxs(n, tier, pat, nxt), n)

# ---- step5 equal-n comparison (seeded draws; with 35 WF24 stems the same seeds draw different animals)
o, n = one(L + "step5.*.log"), one(L + "step5_n88.*.log")
for k in ("WF24", "WF22\\+WF23\\+WF24", "WF22"):
    pat = r"^\s+%s\s+mean\s+([\d.]+)%%" % k
    add("step5_equaln_rarefaction.sh", "9 per cohort, 10 draws: %s-only clusters %%" % k.replace("\\", "") if "+" not in k else "9 per cohort, 10 draws: three-cohort clusters %", rx(o, pat), rx(n, pat), n)
for c in ("WF22", "WF23", "WF24"):
    pat = r"^ %s: .*\b9:(\d+)" % c
    add("step5_equaln_rarefaction.sh", "%s clusters at 9 own metagenomes (seed 300)" % c, rx(o, pat), rx(n, pat), n)

# ---- step6b batch table
o, n = one(L + "step6b.*.log"), one(L + "step6b_n88.*.log")
def k1(p, coh): return rxs(p, "KEY CONTRAST 1", r"%s: n=\d+ median_clean=[\d,]+ median_proteins=([\d,]+)" % coh, "KEY CONTRAST 2")
def k2(p, rn, f): return rxs(p, "KEY CONTRAST 2", r"%s: n=\d+ median_clean=([\d,]+) median_proteins=([\d,]+)" % rn, None, f)
add("step6b_batch_table.sh", "shared run median proteins, WF23", k1(o, "WF23"), k1(n, "WF23"), R + "batch_table_88.tsv")
add("step6b_batch_table.sh", "shared run median proteins, WF24", k1(o, "WF24"), k1(n, "WF24"), R + "batch_table_88.tsv")
add("step6b_batch_table.sh", "WF24 run 0730 median proteins", k2(o, "UCB_20250730_M006342", 2), k2(n, "UCB_20250730_M006342", 2), R + "batch_table_88.tsv")
def fold(p, f):
    a, b = k2(p, "UCB_20250730_M006342", f), k2(p, "UCB_20250426_M005990", f)
    return "NOT_FOUND" if "NOT_FOUND" in (a, b) else "%.2f" % (num(a) / num(b))
add("step6b_batch_table.sh", "WF24 run 0730 over 0426, median proteins fold", fold(o, 2), fold(n, 2), n)
add("step6b_batch_table.sh", "WF24 run 0730 over 0426, median clean reads fold", fold(o, 1), fold(n, 1), n)

# ---- step13h GU/EU criteria
o, n = L + "step13h_89.log", L + "step13h_n88.log"
for crit, lab in (("E<=1e-5    qcov>=0   scov>=0", "E<=1e-5 any coverage"), ("E<=1e-10   qcov>=50  scov>=50", "E<=1e-10 qcov>=50 scov>=50"),
                  ("E<=1e-20   qcov>=50  scov>=50", "E<=1e-20 qcov>=50 scov>=50")):
    pat = r"^%s\s+(\d+)\s+(\d+)\s+([\d.]+)%%" % re.escape(crit)
    add("step13h_gu_eu_split.sh", "EU (no hit) n, %s" % lab, rx(o, pat, 2), rx(n, pat, 2), n)
    add("step13h_gu_eu_split.sh", "EU %%, %s" % lab, rx(o, pat, 3), rx(n, pat, 3), n)
    if "any" in lab:
        add("step13h_gu_eu_split.sh", "GU (hit) n, %s" % lab, rx(o, pat, 1), rx(n, pat, 1), n)

# ---- step16t / step16u mapping summaries (no remapping)
o, n = L + "step16t.log", L + "step16t_n88.log"
for c in ("WF22", "WF23", "WF24"):
    pat = r"^%s\tn=(\d+)\tmean=([\d.]+)" % c
    add("step16t_qc_table.sh", "%s mean mapped fraction" % c, rx(o, pat, 2), rx(n, pat, 2), R + "mapping_qc_88.tsv")
    add("step16t_qc_table.sh", "%s samples (n)" % c, rx(o, pat, 1), rx(n, pat, 1), R + "mapping_qc_88.tsv")
n = L + "step16u_n88.log"
for tag in ("primary", "mapq10"):
    add("step16u_matrix.py", "genes detected at least once, %s" % tag, rx(n, r"^%s: 89 genes with any count (\d+)" % tag), rx(n, r"^%s: .*n88 rows \d+, genes with any count (\d+)" % tag),
        R + ("count_matrix_primary_n88.tsv.gz" if tag == "primary" else "count_matrix_mapq10_n88.tsv.gz"))

# ---- step16y completeness by class
o, n = L + "step16y.log", L + "step16y_n88.log"
for k in ("K", "U"):
    add("step16y_complete_by_class.py", "%s complete genes mean bp" % k, rx(o, r"^%s complete n=\d+ mean_bp=(\d+)" % k), rx(n, r"^%s complete n=\d+ mean_bp=(\d+)" % k), R + "completeness_by_class_n88.tsv")
Lc = {}
for ln in open(R + "supported_cds_lengths.tsv"):
    g, l = ln.rstrip("\n").split("\t"); Lc[g] = int(l)
def umean(cf):
    s = c = 0
    with open(cf) as fh:
        next(fh)
        for ln in fh:
            q = ln.rstrip("\n").split("\t")
            if q[1] == "1" and q[5] == "U": s += Lc.get(q[0], 0); c += 1
    return "%.0f" % (s / c)
add("step16y_complete_by_class.py", "U all genes mean bp (old computed from rep_classes.tsv)", umean(R + "rep_classes.tsv"), rx(n, r"^U all n=\d+ mean_bp=(\d+)"), n)

# ---- step17b aggregate test (n88: WF24 34 animals)
o, n = L + "step17b_89.log", L + "step17b_n88.log"
for c in ("WF22", "WF23", "WF24"):
    for i, lab in ((0, "raw U%"), (1, "length-normalized U%")):
        ps = lambda p: re.findall(r"--> (?:Kruskal-Wallis H|Mann-Whitney U)=[\d.]+ p=([\d.]+)", sect(p, "=== %s ===" % c, "=== WF2" if c != "WF24" else None))
        a, b = ps(o), ps(n)
        add("step17b_class_shift.py", "%s %s by treatment p" % (c, lab), a[i] if len(a) > i else "NOT_FOUND", b[i] if len(b) > i else "NOT_FOUND", n)

# rows computed above but not cited in README are left out of the table
NOTCITED = re.compile(r"at 8,000,000 proteins|^B_naive (egg_mass|month)|^A_blocked month2|^C_EM3 month|Fungi representatives|^supported (eggnog|pfam|union)|"
                      r"in all 3 years|% of mapped reads|^n samples$|section D|^WF23 genes in|^WF24 genes in >=30|^WF24 genes variable|^(K|KWP) (n|% of supported)$")
rows = [r for r in rows if not NOTCITED.search(r[1]) and not (r[0].startswith("step22b") and r[1].startswith("genes in"))]
out = R + "n88_number_table.tsv"
with open(out, "w") as fh:
    fh.write("script\tquantity\told_value\tnew_value\toutput_file\n")
    for r in rows: fh.write("\t".join(str(x) for x in r) + "\n")
nf = sum(1 for r in rows if "NOT_FOUND" in (r[2], r[3]))
print("rows:", len(rows), " rows with NOT_FOUND:", nf, " wrote", out)
