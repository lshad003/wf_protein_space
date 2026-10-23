# Figure S3. Step 34j ESMFold models of the recurrent EU families (4,250 inputs) and length-matched known controls (4,509).
# a, mean C-alpha pLDDT per model; b, models with a Foldseek hit in any of PDB, CATH50, AFDB50, ESMAtlas30 by pLDDT bin,
# at E <= 1e-3 and at E <= 1e-3 with qTM >= 0.5; c, the same per database, all models (AFDB50 annotated: one of the 25
# best AFDB50 hits has an annotated UniProt entry, ProtNLM-only names counted as unannotated).
# Inputs: results/step34j_structures_v2.tsv, results/step34j/hit_share_by_db_v2.tsv. Output: figures/figS3_structures.*
import sys, csv
import numpy as np
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts")
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
fs.setup(); R = fs.W + "/results/"
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
COL = {"EU": "#AA3377", "control": "#4477AA"}; NAME = {"EU": "Recurrent EU families", "control": "Known, length-matched"}
S = [r for r in rows("step34j_structures_v2.tsv") if r["mean_ca_plddt"] != "NA"]
HS = rows("step34j/hit_share_by_db_v2.tsv")
def get(c, db, s, b): return next(r for r in HS if r["criterion"] == c and r["database"] == db and r["set"] == s and r["plddt_bin"] == b)
fig, axs = plt.subplots(1, 3, figsize=(7.08, 2.2), gridspec_kw=dict(width_ratios=[1, 1, 1.6], wspace=0.45))
ax = axs[0]; bins = np.arange(20, 101, 2.5)
for s in ("control", "EU"):
    v = [float(r["mean_ca_plddt"]) for r in S if r["set"] == s]
    ax.hist(v, bins=bins, histtype="step", color=COL[s], linewidth=0.9, label="%s (n = %s)" % (NAME[s], format(len(v), ",")))
for x in (50, 70): ax.axvline(x, color="#555555", linewidth=0.5, linestyle=(0, (2, 2)))
ax.set_ylim(0, 680); ax.set_xlabel("Mean Cα pLDDT"); ax.set_ylabel("Models"); ax.legend(frameon=False, fontsize=5, loc="upper left"); fs.letter(ax, "a")
ax = axs[1]; PB = [("<50", "<50"), ("50-70", "50–70"), (">=70", "≥70")]; x = np.arange(3)
for s in ("EU", "control"):
    for c, ls in (("E<=1e-3", "-"), ("E<=1e-3,qTM>=0.5", (0, (2, 1.5)))):
        y = [float(get(c, "any", s, b)["pct"]) for b, _ in PB]
        ax.plot(x, y, color=COL[s], linestyle=ls, marker="o", markersize=2.5, linewidth=0.9)
    ax.text(2.12, float(get("E<=1e-3", "any", s, ">=70")["pct"]) + (4 if s == "control" else 0), "Known" if s == "control" else "EU", color=COL[s], fontsize=5, va="center")
ax.set_xticks(x); ax.set_xticklabels([l for _, l in PB]); ax.set_xlabel("Mean Cα pLDDT"); ax.set_ylabel("Models with a hit, any database (%)")
ax.set_ylim(0, 100); ax.set_xlim(-0.3, 2.6)
ax.legend(handles=[plt.Line2D([], [], color="#333333", linewidth=0.9, label="E ≤ 10$^{-3}$"),
                   plt.Line2D([], [], color="#333333", linewidth=0.9, linestyle=(0, (2, 1.5)), label="and qTM ≥ 0.5")], frameon=False, fontsize=5, loc="upper left")
fs.letter(ax, "b")
ax = axs[2]; DBL = [("PDB", "PDB"), ("CATH50", "CATH50"), ("AFDB50", "AFDB50"), ("AFDB50_annotated", "AFDB50\nannotated"), ("ESMAtlas30", "ESM\nAtlas30"), ("any", "any")]
x = np.arange(len(DBL)); w = 0.2
for k, (s, c) in enumerate((("EU", "E<=1e-3"), ("EU", "E<=1e-3,qTM>=0.5"), ("control", "E<=1e-3"), ("control", "E<=1e-3,qTM>=0.5"))):
    y = [float(get(c, db, s, "all")["pct"]) for db, _ in DBL]
    ax.bar(x + (k - 1.5) * w, y, w * 0.95, color=COL[s] if c == "E<=1e-3" else "white", edgecolor=COL[s], linewidth=0.6)
ax.set_xticks(x); ax.set_xticklabels([l for _, l in DBL], fontsize=5); ax.set_ylabel("Models with a hit (%)"); ax.set_ylim(0, 100)
ax.legend(handles=[Patch(facecolor=COL["EU"], label="EU, E ≤ 10$^{-3}$"), Patch(facecolor="white", edgecolor=COL["EU"], label="EU, and qTM ≥ 0.5"),
                   Patch(facecolor=COL["control"], label="Known, E ≤ 10$^{-3}$"), Patch(facecolor="white", edgecolor=COL["control"], label="Known, and qTM ≥ 0.5")],
          frameon=False, fontsize=5, loc="upper left", ncol=2)
fs.letter(ax, "c")
fs.save(fig, "figS3_structures")
