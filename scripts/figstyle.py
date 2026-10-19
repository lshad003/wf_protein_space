# Shared figure style: 7.0 in wide, 6/7/5.5 pt type, Arial or the closest available sans; PDF text kept editable.
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
def setup():
    avail = {f.name for f in fm.fontManager.ttflist}
    name = next(n for n in ("Arial", "Liberation Sans", "Helvetica", "Nimbus Sans", "DejaVu Sans") if n in avail)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [name], "font.size": 6,
        "axes.titlesize": 7, "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False})
    print("font:", name)
def letter(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=7, fontweight="bold", va="bottom", ha="left")
def save(fig, stem):
    os.makedirs(W + "/figures", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig("%s/figures/%s.%s" % (W, stem, ext), dpi=300, bbox_inches="tight")
    print("wrote figures/%s.pdf and .png" % stem)
