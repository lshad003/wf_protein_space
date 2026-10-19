# Figure 2: the microbial unknown fraction. a, class shares of the strict catalog (K split into known-function and
# DUF/UPF-only); b, protein length by class (CDS length / 3); c, ORF completeness by class (Prodigal partial code);
# d, read-based prevalence by class (>= 1 read), with the share detected in all three cohorts at >= 1 and >= 5 reads;
# e, member-level audit of EU clusters (Step 34h); f, NMPFams2 comparison (Step 34c, Step 34h).
# Inputs: results/step34a_rep_origin_v2_final.tsv, rep_classes_n88.tsv, fourway_classes_n88.tsv,
# prevalence_primary_n88.tsv.gz, step34b2_catalog_stats_final.tsv, step34h_member_classes.tsv,
# step34c_strict_eu_nmpfams2.tsv, step34c2c_strict_v2_families.tsv. Output: figures/fig2_unknown_fraction.pdf/.png
import sys, csv, gzip
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
fs.setup(); R = fs.W + "/results/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
def tab(path, cols, op=open):
    with op(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = [h.index(c) for c in cols]
        for l in f:
            p = l.rstrip("\n").split("\t"); yield [p[i] for i in ix]
CL = ["K", "KWP", "GU", "EU"]; CC = {"K": "#4477AA", "KWP": "#66CCEE", "GU": "#CCBB44", "EU": "#AA3377"}
cls = {r: c for r, c, s in tab(R + "step34a_rep_origin_v2_final.tsv", ["rep", "class4", "strict_final"]) if s in ("1", "True", "true", "yes")}
if len(cls) != 1359273: raise SystemExit("ABORT: strict_final reps %d" % len(cls))
duf = {r for r, p, nd in tab(R + "rep_classes_n88.tsv", ["rep", "pfam", "pfam_nonduf"]) if p == "1" and nd == "0" and cls.get(r) == "K"}
L = {c: [] for c in CL}; P = {c: {} for c in CL}
for r, ln, co in tab(R + "fourway_classes_n88.tsv", ["rep", "length", "completeness"]):
    c = cls.get(r)
    if c: L[c].append(int(ln) / 3.0); P[c][co] = P[c].get(co, 0) + 1
NS = {c: [] for c in CL}
for r, n in tab(R + "prevalence_primary_n88.tsv.gz", ["gene", "n_samples"], gzip.open):
    c = cls.get(r)
    if c: NS[c].append(int(n))
st = {r["class"]: r for r in rows("step34b2_catalog_stats_final.tsv") if r["set"] == "strict_final"}
NC = {c: sum(1 for v in cls.values() if v == c) for c in CL}
print("classes:", NC, "| DUF/UPF-only K:", len(duf))
print("lengths:", {c: len(L[c]) for c in CL}, "| prevalence rows:", {c: len(NS[c]) for c in CL})
fig = plt.figure(figsize=(7.08, 5.8)); g = fig.add_gridspec(2, 3, width_ratios=[0.75, 1, 1], hspace=0.6, wspace=0.6)
tot = len(cls); ax = fig.add_subplot(g[0, 0])
seg = [("K, known-function domain", NC["K"] - len(duf), CC["K"], None), ("K, DUF/UPF only", len(duf), CC["K"], "////"),
       ("KWP", NC["KWP"], CC["KWP"], None), ("GU", NC["GU"], CC["GU"], None), ("EU", NC["EU"], CC["EU"], None)]
y = 0; cen = []
for lab, n, col, ha in seg:
    h = 100.0 * n / tot; ax.bar(0, h, bottom=y, width=0.6, color=col, hatch=ha, edgecolor="white" if ha else col, linewidth=0)
    cen.append(y + h / 2); y += h
ly = list(cen)
for i in range(len(ly) - 2, -1, -1): ly[i] = min(ly[i], ly[i + 1] - 9)
for (lab, n, col, ha), c0, c1 in zip(seg, cen, ly):
    ax.plot([0.32, 0.45], [c0, c1], color="#777777", linewidth=0.4)
    ax.text(0.48, c1, "%s\n%s (%.2f%%)" % (lab, format(n, ","), 100.0 * n / tot), va="center", fontsize=4.8)
ax.set_xlim(-0.4, 2.8); ax.set_ylim(0, 100); ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
ax.set_ylabel("% of microbial representatives"); fs.letter(ax, "a")
ax = fig.add_subplot(g[0, 1])
for c in CL:
    v = np.sort(np.array(L[c])); ax.plot(v, 100.0 * np.arange(1, len(v) + 1) / len(v), color=CC[c], linewidth=0.9,
                                         label="%s (median %d aa)" % (c, int(np.median(v))))
ax.set_xscale("log"); ax.set_xlabel("protein length (aa)"); ax.set_ylabel("cumulative % of representatives")
ax.set_ylim(0, 100); ax.legend(frameon=False, fontsize=5, loc="lower right"); fs.letter(ax, "b")
ax = fig.add_subplot(g[0, 2]); bot = np.zeros(4)
GRP = [("complete", ["00"], "#333333"), ("one end truncated", ["01", "10"], "#999999"), ("both ends truncated", ["11"], "#DDDDDD")]
for lab, codes, col in GRP:
    v = np.array([100.0 * sum(P[c].get(k, 0) for k in codes) / sum(P[c].values()) for c in CL])
    ax.bar(range(4), v, bottom=bot, color=col, width=0.6, linewidth=0, label=lab); bot += v
for i, c in enumerate(CL): ax.text(i, 101, "%.1f%%" % (100.0 * P[c].get("00", 0) / sum(P[c].values())), ha="center", fontsize=5)
ax.set_xticks(range(4)); ax.set_xticklabels(CL); ax.set_ylim(0, 110); ax.set_ylabel("% of representatives")
ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1); fs.letter(ax, "c")
ax = fig.add_subplot(g[1, 0])
vp = ax.violinplot([NS[c] for c in CL], positions=range(4), widths=0.75, showmedians=True, showextrema=False)
for b, c in zip(vp["bodies"], CL): b.set_facecolor(CC[c]); b.set_edgecolor("none"); b.set_alpha(0.8)
vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(0.8)
for i, c in enumerate(CL):
    ax.text(i, 95, "%.0f%%\n%.0f%%" % (float(st[c]["pct_3cohorts_ge1read"]), float(st[c]["pct_3cohorts_ge5reads"])), ha="center", fontsize=4.8)
ax.text(1.5, 106, "% detected in all 3 cohorts (\u22651 / \u22655 reads)", ha="center", fontsize=4.8, color="#555555")
ax.set_xticks(range(4)); ax.set_xticklabels(CL); ax.set_ylim(0, 110); ax.set_ylabel("metagenomes detected (\u22651 read)"); fs.letter(ax, "d")
mc = rows("step34h_member_classes.tsv"); print("member-audit rows:", len(mc))
any_ = {k: sum(1 for r in mc if r["class_any"] == k) for k in CL}; maj = {k: sum(1 for r in mc if r["class_majority"] == k) for k in CL}
ax = fig.add_subplot(g[1, 1]); x = np.arange(3); w = 0.36
for k, (lab, d, al) in enumerate((("any distinct sequence", any_, 1.0), ("majority of sequences", maj, 0.45))):
    v = [d[c] for c in ("K", "KWP", "GU")]
    ax.bar(x + (k - 0.5) * w, v, w, color=[CC[c] for c in ("K", "KWP", "GU")], alpha=al, linewidth=0, label=lab)
    for xi, vi in zip(x + (k - 0.5) * w, v): ax.text(xi, vi + 8, str(vi), ha="center", fontsize=5)
ax.set_xticks(x); ax.set_xticklabels(["\u2192 K", "\u2192 KWP", "\u2192 GU"]); ax.set_ylabel("EU clusters reassigned")
ax.text(0.02, 0.97, "EU remaining: %s (%.2f%%) any; %s (%.2f%%) majority" % (format(any_["EU"], ","), 100.0 * any_["EU"] / tot,
        format(maj["EU"], ","), 100.0 * maj["EU"] / tot), transform=ax.transAxes, fontsize=4.8, va="top")
ax.set_ylim(0, max(any_["GU"], 1) * 1.3); ax.legend(frameon=False, fontsize=4.8, loc="upper left", bbox_to_anchor=(0.0, 0.88)); fs.letter(ax, "e")
nm = {r["rep"]: r for r in rows("step34c_strict_eu_nmpfams2.tsv")}; cm = {r["rep"]: r for r in mc}
fam = {r["family_rep"] for r in rows("step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"}
srch = [x for x in cm if nm.get(x, {}).get("in_step26_set") == "1"]
SETS = [("all searched", srch), ("still EU, any rule", [x for x in srch if cm[x]["class_any"] == "EU"]),
        ("still EU, majority rule", [x for x in srch if cm[x]["class_majority"] == "EU"]), ("recurrent families", [x for x in srch if nm[x]["family_rep"] in fam])]
ax = fig.add_subplot(g[1, 2])
for i, (lab, s) in enumerate(SETS):
    u = sum(1 for x in s if nm[x]["matched"] == "0"); pu = 100.0 * u / len(s)
    ax.barh(i, pu, color=CC["EU"], height=0.6, linewidth=0); ax.barh(i, 100 - pu, left=pu, color="#BBBBBB", height=0.6, linewidth=0)
    ax.text(2, i, "%.2f%% of %s" % (pu, format(len(s), ",")), va="center", color="white", fontsize=5); print("NMPFams2 %s: %d / %d unmatched" % (lab, u, len(s)))
ax.set_yticks(range(len(SETS))); ax.set_yticklabels([s[0] for s in SETS], fontsize=5); ax.invert_yaxis(); ax.set_xlim(0, 100)
ax.set_xlabel("% without an NMPFams2 match"); fs.letter(ax, "f")
fs.save(fig, "fig2_unknown_fraction")
