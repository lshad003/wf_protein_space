# Render for Figure 4c: ESMFold model of the Cse2/CasB-like protein (UHM45.35758__k141_100597_236) superposed on
# E. coli CasB (PDB 5H9E chain C) by CE alignment; cartoon, ray-traced, transparent background.
# Run on a compute node: pymol -cq fig_render_cse2.py. Output: results/fig_render_cse2_varA.png
from pymol import cmd
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
cmd.load(W + "/results/esmfold_pdb/UHM45.35758__k141_100597_236.pdb", "cse2")
cmd.load(W + "/results/structures/5h9e.cif", "casb")
cmd.remove("casb and not chain C"); cmd.remove("solvent or inorganic")
print("cealign (cse2 onto CasB):", cmd.cealign("casb", "cse2"))
cmd.hide("everything"); cmd.show("cartoon"); cmd.set("cartoon_fancy_helices", 1); cmd.set("cartoon_discrete_colors", 1)
cmd.set_color("wfmag", [0.667, 0.200, 0.467]); cmd.color("grey80", "casb"); cmd.color("wfmag", "cse2")
cmd.set("cartoon_transparency", 0.3, "casb")
cmd.bg_color("white"); cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0); cmd.set("depth_cue", 0)
cmd.set("ray_trace_mode", 1); cmd.set("ray_trace_gain", 0.3); cmd.set("antialias", 2); cmd.set("specular", 0.15)
cmd.orient("cse2 and ss H+S"); cmd.zoom("cse2 or casb", 2)
cmd.png(W + "/results/fig_render_cse2_varA.png", width=2400, height=2400, dpi=600, ray=1)
print("wrote results/fig_render_cse2_varA.png")
