# Shared figure style. Okabe-Ito, fixed class colours, used by every figure.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
CLS = ["K","KWP","GU","EU"]
COL = {"K":"#0072B2","KWP":"#56B4E9","GU":"#E69F00","EU":"#D55E00"}
COH = {"WF22":"#009E73","WF23":"#CC79A7","WF24":"#F0E442"}
LAB = {"K":"K\nPfam domain","KWP":"KWP\neggNOG only",
       "GU":"GU\nin nr, no function","EU":"EU\nno database hit"}
plt.rcParams.update({
    "font.size":8,"axes.titlesize":9,"axes.labelsize":8,
    "xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,
    "axes.spines.top":False,"axes.spines.right":False,
    "figure.dpi":300,"savefig.bbox":"tight"})
MM = 1/25.4
def size(w_mm,h_mm): return (w_mm*MM, h_mm*MM)
