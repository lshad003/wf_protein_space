# Step 31a: input for the structure screen of the single-variant defence-context families (Step 30b).
# Selection: rows of results/step30b_single_variant_context.tsv with set "unknown" whose most common +-3 neighbour set
#   (modal_pm3_set) carries an annotation whose only category in results/step30b_category_map.tsv is
#   "CRISPR/Cas and other defence" (the rule step30b used for categories_in_context; cross-checked against that
#   column). Aborts unless exactly 84.
# Sequences: family rep from results/prostt5_main_eu15580.faa. Terminal '*' removed as in step29a (ESMFold alphabet has
#   no stop); aborts on an internal stop. Length = residues without '*'.
# Reps folded in Step 29 (results/esmfold_pdb/<rep>.pdb) are copied to results/def84_step29_pdb/ and not refolded;
#   the rest go to results/def84_fold.faa for step31b.
# Outputs (refuses to overwrite): results/def84_input.faa (all 84), results/def84_labels.tsv, results/def84_fold.faa,
#   results/def84_step29_pdb/
import os, re, shutil
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
POS = "UHM636.41030__k141_332084_5"; DEF = "CRISPR/Cas and other defence"
OUT = [R + "def84_input.faa", R + "def84_labels.tsv", R + "def84_fold.faa", R + "def84_step29_pdb"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
CR = []
with open(R + "step30b_category_map.tsv") as fh:
    next(fh)
    for l in fh:
        c, rx = l.rstrip("\n").split("\t")[:2]
        if c != "other": CR.append((c, re.compile(rx)))
def is_def(a):
    hits = {c for c, rx in CR if rx.search(a.split(":", 1)[1])}
    return hits == {DEF}
sel, disagree = [], 0
with open(R + "step30b_single_variant_context.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t")))
        if f["set"] != "unknown": continue
        ms = f["modal_pm3_set(member_frac_pm3)"]
        items = [] if ms == "-" else [x.rpartition("(") for x in ms.split(";")]
        d = [(a, fr.rstrip(")")) for a, _, fr in items if is_def(a)]
        if bool(d) != (DEF in f["categories_in_context"].split(";")): disagree += 1
        if d: sel.append((f["family_rep"], int(f["n_animals"]), f["euk_flag"], ";".join("%s(%s)" % x for x in d)))
print("selected %d | disagreements with categories_in_context %d" % (len(sel), disagree))
if len(sel) != 84 or disagree: raise SystemExit("ABORT: expected exactly 84 with no disagreement")
seq, cur = {}, None
for l in open(R + "prostt5_main_eu15580.faa"):
    if l[0] == ">": cur = l[1:].split()[0]; seq[cur] = []
    else: seq[cur].append(l.strip())
seq = {k: "".join(v).rstrip("*") for k, v in seq.items()}
miss = [r for r, *_ in sel if r not in seq]
if miss: raise SystemExit("ABORT: %d reps not in prostt5_main_eu15580.faa: %s" % (len(miss), miss[:5]))
stop = [r for r, *_ in sel if "*" in seq[r]]
if stop: raise SystemExit("ABORT: internal stop in %s" % stop[:5])
os.mkdir(OUT[3]); prev = []
with open(OUT[0], "w") as fa, open(OUT[1], "w") as lb, open(OUT[2], "w") as ff:
    lb.write("rep\tlength\tn_animals\teuk_flag\tdefence_neighbours(member_frac_pm3)\tstep29_pdb\n")
    for r, na, eu, d in sel:
        s = seq[r]; fa.write(">%s\n%s\n" % (r, s))
        p = R + "esmfold_pdb/" + r + ".pdb"
        if os.path.exists(p): shutil.copy2(p, OUT[3]); prev.append(r)
        else: ff.write(">%s\n%s\n" % (r, s))
        lb.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (r, len(s), na, eu, d, "yes" if r in prev else "no"))
L = sorted(len(seq[r]) for r, *_ in sel)
print("lengths: min %d, median %d, max %d | euk-flagged %d" % (L[0], L[len(L) // 2], L[-1], sum(x[2] == "yes" for x in sel)))
print("positive control %s included: %s" % (POS, "yes" if POS in {x[0] for x in sel} else "no"))
print("already folded in Step 29 (copied to %s): %d" % (OUT[3], len(prev)))
for r in prev: print("  " + r)
print("to fold (results/def84_fold.faa): %d" % (84 - len(prev)))
