# Figure 4. a, best ProstT5/Foldseek E-value per sequence (PDB or CATH50) for recurrent EU families, length-matched known
# proteins and reversed decoys; b, structural match rate by length (Wilson 95% CI); c, ESMFold model of the
# Cse2/CasB-like protein superposed on E. coli CasB (PDB 5H9E chain C), PyMOL render; d, its type I-E cas locus;
# e, Step 34j Foldseek hit class of the ESMFold models (EU family models vs controls) by mean C-alpha pLDDT bin
# (results/step34j/models_by_set_plddt_v2.tsv; ProtNLM-only UniProt names counted as unannotated).
# 2026-10-23: Table S6 (results/supp_fold_candidates.tsv) is no longer written here; scripts/step34j_d2_tables.py
# rebuilds it from the Step 34j representative models.
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
fs.setup(); R = fs.W + "/results/"; F = plt.rcParams["font.sans-serif"][0]
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": F, "mathtext.it": F + ":italic", "mathtext.bf": F + ":bold"})
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
COL = {"queries": "#AA3377", "K": "#4477AA", "decoys": "#999999"}
NAME = {"queries": "Recurrent EU families", "K": "Known, length-matched", "decoys": "Reversed decoys"}
fig = plt.figure(figsize=(7.08, 4.9))
g = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.2], height_ratios=[1, 0.75], hspace=0.6, wspace=0.5)
be = rows("fig_prep_bestE.tsv"); xs = np.linspace(0, 20, 401); ax = fig.add_subplot(g[0, 0])
for s in ("K", "queries", "decoys"):
    v = np.array([-np.log10(float(r["bestE"])) if r["bestE"] != "NA" else -1.0 for r in be if r["set"] == s])
    ax.plot(xs, [100.0 * (v >= x).mean() for x in xs], color=COL[s], linewidth=0.9, label=NAME[s])
ax.axvline(3, color="#555555", linewidth=0.5, linestyle=(0, (2, 2)))
ax.set_xlabel("$-$log$_{10}$ best E-value"); ax.set_ylabel("Sequences at or beyond (%)"); ax.set_xlim(0, 20); ax.set_ylim(0, 70)
ax.legend(frameon=False, fontsize=5, loc="upper right", handlelength=1.2); fs.letter(ax, "a")
bl = rows("step34d2_by_length_fixed.tsv"); SN = sorted({r["set"] for r in bl})
pick = {"queries": next(n for n in SN if "quer" in n), "K": next(n for n in SN if n.startswith("K")), "decoys": next(n for n in SN if "decoy" in n)}
BIN = ["<100", "100-199", "200-399", ">=400"]; ax = fig.add_subplot(g[0, 1])
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)
for s in ("K", "queries", "decoys"):
    rr = {r["len_bin_aa"]: r for r in bl if r["set"] == pick[s]}; k = [int(rr[b]["n_hit"]) for b in BIN]; n = [int(rr[b]["n"]) for b in BIN]
    y = [100.0 * a / b for a, b in zip(k, n)]; ci = [wilson(a, b) for a, b in zip(k, n)]
    ax.errorbar(range(4), y, yerr=[[v - c[0] for v, c in zip(y, ci)], [c[1] - v for v, c in zip(y, ci)]], color=COL[s], marker="o",
                markersize=2.5, linewidth=0.9, capsize=1.2, elinewidth=0.5)
ax.set_xticks(range(4)); ax.set_xticklabels(["<100", "100\u2013199", "200\u2013399", "\u2265400"]); ax.set_xlabel("Protein length (aa)")
ax.set_ylabel("Structural match, E \u2264 10$^{-3}$ (%)"); ax.set_ylim(0, 100); fs.letter(ax, "b")
import re
e29 = None
for r in rows("step34h_cse2_varA.tsv"):
    t = (r["item"] + " " + r["value"]).lower()
    if "pdb" in t and "cath" not in t and "qtm" in t:
        m = re.search(r"qtm\s+([0-9.]+(?:e[-+]?[0-9]+)?)", t)
        if m: e29 = {"qtmscore": m.group(1)}; print("using row:", r["item"], "|", r["value"]); break
if e29 is None: raise SystemExit("ABORT: variant A 5H9E qTM not found in results/step34h_cse2_varA.tsv")
print("Cse2 variant A vs 5H9E qTM:", float(e29["qtmscore"]))
if not 0 < float(e29["qtmscore"]) <= 1: raise SystemExit("ABORT: qTM out of range")
im = plt.imread(R + "fig_render_cse2_varA.png"); yy, xx = np.where(im[:, :, 3] > 0.01); im = im[yy.min():yy.max() + 1, xx.min():xx.max() + 1]
ax = fig.add_subplot(g[0, 2]); ax.imshow(im); ax.axis("off"); fs.letter(ax, "c")
ax.text(0.0, -0.08, "CasB, $\\it{E.\\ coli}$ (PDB 5H9E)", transform=ax.transAxes, fontsize=5.5, color="#777777")
ax.text(0.0, -0.17, "Cse2/CasB-like, this study (TM-score %.2f)" % float(e29["qtmscore"]), transform=ax.transAxes, fontsize=5.5, color="#AA3377")
L = rows("fig_prep_cse2_locus_best.tsv"); print("locus columns:", list(L[0].keys()))
gc = next(c for c in L[0] if c.lower() in ("gene", "gene_id", "id")); cn = next((c for c in L[0] if "cas" in c.lower()), None)
def num(r):
    try: return int(r[gc].rsplit("_", 1)[1])
    except ValueError: return -1
c2 = next(r for r in L if r[gc].endswith("_259")); n2 = num(c2)
lo = [r for r in L if n2 - 4 <= num(r) <= n2 + 3]
neg = c2["strand"].strip() in ("-", "-1"); E0 = max(int(r["end"]) for r in lo); S0 = min(int(r["start"]) for r in lo); xe = E0 - S0
def tr(r):
    a, b = int(r["start"]), int(r["end"]); st = 1 if r["strand"].strip() in ("+", "1") else -1
    return (E0 - b, E0 - a, -st) if neg else (a - S0, b - S0, st)
ax = fig.add_subplot(g[1, 0:2]); ax.plot([0, xe], [0, 0], color="#333333", linewidth=0.5, zorder=1); jl = 0
for r in sorted(lo, key=lambda r: tr(r)[0]):
    a, b, st = tr(r); hd = min(120, 0.45 * (b - a)); h = 0.34
    pts = [(a, -h / 2), (b - hd, -h / 2), (b, 0), (b - hd, h / 2), (a, h / 2)] if st > 0 else [(b, -h / 2), (a + hd, -h / 2), (a, 0), (a + hd, h / 2), (b, h / 2)]
    nm = (r.get(cn) or "").strip() if cn else ""
    nm = "" if nm in ("NA", "-", "none") else nm
    is2 = r is c2; col = "#AA3377" if is2 else ("#E69F00" if nm else "#E3E3E3")
    ax.add_patch(Polygon(pts, closed=True, facecolor=col, edgecolor="#333333", linewidth=0.4, zorder=2))
    if nm or is2:
        yl, va = (h / 2 + 0.08, "bottom") if jl % 2 == 0 else (-h / 2 - 0.08, "top"); jl += 1
        ax.text((a + b) / 2, yl, "cse2-like" if is2 else nm, ha="center", va=va, fontsize=5.5, fontstyle="italic",
                color="#AA3377" if is2 else "black", fontweight="bold" if is2 else "normal")
print("locus genes:", [(num(r), (r.get(cn) or "") if cn else "") for r in lo])
ax.plot([0, 1000], [-0.62, -0.62], color="black", linewidth=0.7); ax.text(500, -0.67, "1 kb", ha="center", va="top", fontsize=5)
ax.set_xlim(-60, xe + 60); ax.set_ylim(-0.8, 0.6); ax.axis("off"); fs.letter(ax, "d")
HC = [("experimental", "PDB or CATH50", "#333333"), ("AFDB_annotated", "AFDB50, annotated", "#777777"),
      ("AFDB_unannotated_only", "AFDB50, unannotated only", "#BBBBBB"), ("ESMAtlas_only", "ESM Atlas only", "#E3E3E3"), ("none", "no hit", "#FFFFFF")]
MB = [r for r in csv.DictReader(open(R + "step34j/models_by_set_plddt_v2.tsv"), delimiter="\t") if r["protnlm_rule"] == "unannotated"]
PB = [(">=70", "\u226570"), ("50-70", "50\u201370"), ("<50", "<50")]
ax = fig.add_subplot(g[1, 2]); xs, xl = [], []
for gi, (st, sl) in enumerate((("EU", "EU"), ("control", "Known"))):
    for bi, (b, bl) in enumerate(PB):
        r = next(r for r in MB if r["set"] == st and r["plddt_bin"] == b); n = int(r["n_models"]); x = gi * 3.6 + bi; bot = 0
        for k, _, c in HC:
            v = 100.0 * int(r[k]) / n; ax.bar(x, v, 0.8, bottom=bot, color=c, edgecolor=COL["queries" if st == "EU" else "K"], linewidth=0.5); bot += v
        ax.text(x, 101.5, format(n, ","), ha="center", va="bottom", fontsize=4.5)
        xs.append(x); xl.append(bl)
        print("Fig. 4e %s pLDDT %s: n %d, %s" % (st, b, n, ", ".join("%s %s" % (k, r[k]) for k, _, _ in HC)))
ax.set_xticks(xs); ax.set_xticklabels(xl, fontsize=5); ax.set_ylim(0, 108); ax.set_ylabel("Models (%)")
for gi, (st, sl) in enumerate((("EU", "EU families"), ("control", "Known"))):
    ax.text(gi * 3.6 + 1, -24, sl, ha="center", va="top", fontsize=5.5, color=COL["queries" if st == "EU" else "K"])
ax.text(3.1, -13, "pLDDT", ha="center", va="top", fontsize=5)
ax.legend(handles=[Patch(facecolor=c, edgecolor="#333333", linewidth=0.4, label=l) for _, l, c in HC], frameon=False, fontsize=5,
          loc="upper left", bbox_to_anchor=(1.0, 1.0), handlelength=1.0); fs.letter(ax, "e")
fs.save(fig, "fig4_structure")
