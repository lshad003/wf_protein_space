# Step 36b A: the 8 structure-confident recurrent EU families whose Step 36 conserved context is "several_variants"
#   (results/step36/structure_vs_context.tsv): one summary row each and a +-5-gene locus plot (Cse2 style of Fig. 4d).
# Members: strict members in the 88 metagenomes (results/step34h/strict_members.tsv, as Step 36). Window: genes -5..+5
#   by Prodigal gene number on the member's contig, both strands (Prodigal headers, as Step 36), oriented so the focal
#   gene points right. Neighbour class and labels exactly as Step 36 (95% rep -> fourway_classes_n88 class4, else
#   rep_classes_n88; K = non-DUF/UPF Pfam names, DUF-only shows its DUF/UPF names, KWP = eggNOG Preferred_name).
# Tracks: one per 95% cluster of the family that has strict members. Arrangement of a member = the ordered tuple of
#   (oriented offset, relative strand, label or class) over -5..+5; each track shows the most common arrangement in that
#   cluster (ties: more genes, then member ID), drawn from the member carrying it (first by ID), with real coordinates.
# Best confident hit: confident_db / confident_target / confident_qtm / confident_target_name of
#   results/step34j_structures_v2.tsv; E = smallest E of that query-target pair in results/step34j/foldseek/<db>.m8
#   (model or its Merizo domains).
# Outputs (refuse to overwrite): results/step36b/several_variant_families.tsv, results/step36b/locus_tracks.tsv,
#   figures/step36b_locus_<rep>.{pdf,png}, figures/step36b_loci_combined.pdf.
import os, sys, re, gzip, glob, collections
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from matplotlib.backends.backend_pdf import PdfPages
W = fs.W; R = W + "/results/"; D = R + "step36b/"
OLD = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly"; NEW = W + "/catalog/by_assembly_new"
OT, OK, OC = D + "several_variant_families.tsv", D + "locus_tracks.tsv", W + "/figures/step36b_loci_combined.pdf"
for p in (OT, OK, OC):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
def contig(g):
    s, rest = g.split("__", 1); c, n = rest.rsplit("_", 1); return s, c, int(n)
SV = [r for r in tsv(R + "step36/structure_vs_context.tsv") if r["context_class"] == "several_variants"]
if len(SV) != 8: raise SystemExit("ABORT: %d several-variant families, expected 8" % len(SV))
F = [r["family_rep"] for r in SV]
for f in F:
    if os.path.exists(W + "/figures/step36b_locus_%s.pdf" % f): raise SystemExit("REFUSING: figure for %s exists" % f)
S88 = {r["stem"] for r in tsv(W + "/metadata/samples_88_treatment.tsv")}
r2f = {r["rep"]: r["family_rep"] for r in tsv(R + "step34c_strict_eu_nmpfams2.tsv") if r["family_rep"] in F}
MEM = collections.defaultdict(list)
for r in tsv(R + "step34h/strict_members.tsv"):
    if r["rep"] in r2f and r["strict"] == "1" and r["member"].split("__", 1)[0] in S88: MEM[r2f[r["rep"]]].append((r["rep"], r["member"]))
st = {r["family_rep"]: r for r in tsv(R + "step34j_structures_v2.tsv") if r["set"] == "EU" and r["family_rep"] in F}
CL = collections.defaultdict(list)
for r in tsv(R + "step36/conserved_labels.tsv"):
    if r["run"] == "EU_family30" and r["unit"] in F: CL[r["unit"]].append(r)
CB = {r["unit"]: r for r in tsv(R + "step36/context_by_family.tsv") if r["run"] == "EU_family30" and r["unit"] in F}
# best-hit E from the Foldseek tables
want = {(st[f]["id"], st[f]["confident_db"], st[f]["confident_target"]) for f in F}
ev = {}
for db in {w[1] for w in want}:
    tg = {(w[0], w[2]) for w in want if w[1] == db}
    for l in open(R + "step34j/foldseek/%s.m8" % db):
        x = l.split("\t", 8); q = re.sub(r"(\.pdb)?(_[A-Z])?$", "", x[0]); q = q.split("__D")[0] if re.search(r"__D\d+$", q) else q
        if (q, x[1]) in tg: k = (q, x[1]); ev[k] = min(ev.get(k, 9e9), float(x[6]))
# genes on member contigs
tc = {contig(m)[:2] for v in MEM.values() for _, m in v}
genes = collections.defaultdict(dict)
for s in sorted({s for s, _ in tc}):
    p = NEW + "/" + s + ".aa.fa.gz"
    if not os.path.exists(p): p = OLD + "/" + s + ".aa.fa.gz"
    with gzip.open(p, "rt") as fh:
        for l in fh:
            if l[0] != ">": continue
            f = l[1:].split(" # "); c, n = f[0].strip().rsplit("_", 1)
            if (s, c) in tc: genes[(s, c)][int(n)] = (int(f[1]), int(f[2]), 1 if f[3] == "1" else -1)
if len(genes) != len(tc): raise SystemExit("ABORT: contigs without headers")
need = {"%s__%s_%d" % (s, c, n) for (s, c), d in genes.items() for n in d}
rep95 = {}
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in need: rep95[m] = r
reps = set(rep95.values()); cls = {}
for r in tsv(R + "fourway_classes_n88.tsv"):
    if r["rep"] in reps: cls[r["rep"]] = r["class4"]
for r in tsv(R + "rep_classes_n88.tsv"):
    if r["rep"] in reps and r["rep"] not in cls: cls[r["rep"]] = r["class"]
pf = collections.defaultdict(set)
for f in glob.glob(R + "pfam_hs/*.tblout"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.split()
        if p[0] in reps: pf[p[0]].add(p[2])
eg = {}
for f in glob.glob(R + "eggnog2/*.emapper.annotations"):
    for l in open(f):
        if l[0] == "#": continue
        p = l.rstrip("\n").split("\t")
        if p[0] in reps: eg[p[0]] = p[8]
def lab(g):
    r = rep95.get(g)
    if r is None: return "no_rep", ""
    c = cls.get(r, "NA")
    if c == "K":
        x = sorted(d for d in pf.get(r, ()) if not (d.startswith("DUF") or d.startswith("UPF")))
        return ("K", "/".join(x[:2]) + ("+" if len(x) > 2 else "")) if x else ("DUF_only", "/".join(sorted(pf.get(r, ()))[:2]))
    if c == "KWP":
        n = eg.get(r, "-"); return ("KWP", n) if n not in ("-", "") else ("KWP_no_name", "")
    return c, ""
COLR = {"focal": "#AA3377", "K": "#E69F00", "KWP": "#56B4E9", "DUF_only": "#EEDD88", "KWP_no_name": "#BBDDEE",
        "GU": "#BBBBBB", "EU": "#E3E3E3", "U": "#E3E3E3", "no_rep": "#FFFFFF", "NA": "#FFFFFF"}
LEG = [("focal", "focal EU gene"), ("K", "K (Pfam)"), ("KWP", "KWP (eggNOG name)"), ("KWP_no_name", "KWP, no name"),
       ("DUF_only", "DUF/UPF only"), ("GU", "GU"), ("EU", "EU or unclassified")]
def window(m):
    s, c, n = contig(m); g = genes[(s, c)]; a0, b0, sd = g[n]; out = []
    for off in range(-5, 6):
        if n + off not in g: continue
        a, b, s2 = g[n + off]; gid = "%s__%s_%d" % (s, c, n + off)
        k, t = ("focal", "") if off == 0 else lab(gid)
        out.append(dict(off=off * sd, rel=s2 * sd, cls=k, lab=t, a=a, b=b, gene=gid))
    out.sort(key=lambda x: x["off"]); return out, sd
def key(w): return tuple((x["off"], x["rel"], x["lab"] or x["cls"]) for x in w)
rows, tracks = [], []
pdf = PdfPages(OC)
for f, sv in zip(F, SV):
    s = st[f]; cb = CB[f]; byc = collections.defaultdict(list)
    for r, m in MEM[f]: byc[r].append(m)
    E = ev.get((s["id"], s["confident_target"]))
    cons = sorted(CL[f], key=lambda r: (-float(r["frac"]), r["label"]))
    rows.append(dict(family_rep=f, query=s["id"], confident_db=s["confident_db"], confident_target=s["confident_target"],
        confident_target_name=s["confident_target_name"], confident_qtm=s["confident_qtm"], confident_evalue="%.2e" % E if E is not None else "NA",
        mean_ca_plddt=s["mean_ca_plddt"], n_member_contigs=cb["n_member_contigs"], n_animals=cb["n_animals"], n_rep95=cb["n_rep95"],
        conserved_labels=";".join("%s %s (%s/%s contigs), modal %+d, %s cluster(s)" % (r["label"], r["frac"], r["n_contigs"], r["n_member_contigs"],
                                  int(r["modal_offset"]), r["n_rep95"]) for r in cons)))
    T = []
    for r in sorted(byc, key=lambda r: (-len(byc[r]), r)):
        ws = {m: window(m) for m in byc[r]}
        kc = collections.Counter(key(w) for w, _ in ws.values())
        best = sorted(kc.items(), key=lambda t: (-t[1], -len(t[0]), min(m for m in byc[r] if key(ws[m][0]) == t[0])))[0][0]
        mm = sorted(m for m in byc[r] if key(ws[m][0]) == best)[0]
        w, sd = ws[mm]; T.append((r, mm, w, sd, kc[best], len(byc[r])))
        tracks.append(dict(family_rep=f, rep95=r, n_members=len(byc[r]), n_with_modal=kc[best], n_arrangements=len(kc), shown_member=mm,
                           arrangement=" ".join("%+d%s:%s" % (x["off"], "" if x["rel"] > 0 else "(rev)", x["lab"] or x["cls"]) for x in w)))
    # plot
    nt = len(T); fig = plt.figure(figsize=(7.09, 0.6 + 0.62 * nt)); ax = fig.add_axes([0.17, 0.08, 0.79, 0.84 - 0.04 / max(nt, 1)])
    xmax = 0
    for i, (r, mm, w, sd, nk, nm) in enumerate(T):
        y = -i; foc = [x for x in w if x["off"] == 0][0]
        def tr(a, b): return ((a - foc["a"], b - foc["a"]) if sd > 0 else (foc["b"] - b, foc["b"] - a))
        xs = [tr(x["a"], x["b"]) for x in w]; lo = min(a for a, _ in xs); hi = max(b for _, b in xs); xmax = max(xmax, hi, -lo)
        ax.plot([lo, hi], [y, y], color="#333333", linewidth=0.5, zorder=1); jl = 0
        for x, (a, b) in zip(w, xs):
            hd = min(120, 0.45 * (b - a)); h = 0.3; d = x["rel"]
            pts = [(a, y - h / 2), (b - hd, y - h / 2), (b, y), (b - hd, y + h / 2), (a, y + h / 2)] if d > 0 else \
                  [(b, y - h / 2), (a + hd, y - h / 2), (a, y), (a + hd, y + h / 2), (b, y + h / 2)]
            ax.add_patch(Polygon(pts, closed=True, facecolor=COLR.get(x["cls"], "#FFFFFF"), edgecolor="#333333", linewidth=0.4, zorder=2))
            t = "EU (focal)" if x["cls"] == "focal" else x["lab"]
            if t:
                lv = jl % 4; jl += 1   # four label rows (above, below, higher, lower) so labels of short neighbours do not overlap
                yl, va = [(y + h / 2 + 0.03, "bottom"), (y - h / 2 - 0.03, "top"), (y + h / 2 + 0.2, "bottom"), (y - h / 2 - 0.2, "top")][lv]
                ax.text((a + b) / 2, yl, t, ha="center", va=va, fontsize=5, fontstyle="italic" if x["cls"] == "KWP" else "normal",
                        color="#AA3377" if x["cls"] == "focal" else "black", fontweight="bold" if x["cls"] == "focal" else "normal")
        ax.text(-0.01, y, "%s\n%d/%d members" % (r, nk, nm), transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5)
    ax.set_xlim(-xmax - 100, xmax + 100); ax.set_ylim(-nt + 0.4, 0.75); ax.axis("off")
    ax.plot([xmax - 1000, xmax], [-nt + 0.48] * 2, color="black", linewidth=0.7); ax.text(xmax - 500, -nt + 0.43, "1 kb", ha="center", va="top", fontsize=5)
    hit = "%s %s, %s (qTM %s, E %s)" % (s["confident_db"], s["confident_target"], s["confident_target_name"][:80], s["confident_qtm"], rows[-1]["confident_evalue"])
    fig.text(0.01, 0.995, "%s  |  %s contigs, %s animals, %s 95%% clusters  |  best confident hit: %s" % (f, cb["n_member_contigs"], cb["n_animals"], cb["n_rep95"], hit),
             ha="left", va="top", fontsize=6)
    fig.legend(handles=[Patch(facecolor=COLR[k], edgecolor="#333333", linewidth=0.4, label=l) for k, l in LEG], loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=7, frameon=False, fontsize=5, handlelength=1.2)
    fs.save(fig, "step36b_locus_%s" % f); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
pdf.close(); print("wrote", OC)
H = list(rows[0].keys())
with open(OT, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in rows: o.write("\t".join(str(r[k]).replace("\t", " ") for k in H) + "\n")
H = list(tracks[0].keys())
with open(OK, "w") as o:
    o.write("\t".join(H) + "\n")
    for r in tracks: o.write("\t".join(str(r[k]) for k in H) + "\n")
for r in rows: print("\n".join("%s: %s" % (k, v) for k, v in r.items()), "\n")
for r in tracks: print("%s | %s | %d/%d (%d arrangements) | %s" % (r["family_rep"], r["rep95"], r["n_with_modal"], r["n_members"], r["n_arrangements"], r["arrangement"]))
print("wrote", OT, OK)
