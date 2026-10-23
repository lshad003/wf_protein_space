# Step 34j B4b: accessions b4 requests that the gap list (scripts/step34j_b4_uniprot_gap.py) does not hold.
# The gap script ranks the AFDB50 hits of a query unit by (E, -bits, target name); b4 ranks by (E, -bits) in file order,
#   so at a tie on the 25th hit the two keep different targets. This rebuilds b4's requested set exactly
#   (scripts/step34j_b4_classify.py) and lists accessions that are neither in a saved live batch nor in
#   results/step34j/uniprot_unisave/unisave_raw.jsonl.
# Output (refuses to overwrite): results/step34j/uniprot_gap/b4_tie_missing.tsv (input to step34j_b4a_unisave_fetch.py $1).
import os, re, glob, json, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; D = W + "/results/step34j/"; OUT = D + "uniprot_gap/b4_tie_missing.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
units = {os.path.basename(x)[:-4] for x in glob.glob(D + "esmfold_pdb/*.pdb")} | {os.path.basename(x)[:-4] for x in glob.glob(D + "domains/*.pdb")}
def unit(q):
    for c in (q, re.sub(r"\.pdb$", "", q), re.sub(r"(\.pdb)?_[A-Z]$", "", q)):
        if c in units: return c
    raise SystemExit("ABORT: cannot map query %s" % q)
def acc(t):
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", t); return m.group(1) if m else None
byu = collections.defaultdict(list)
for l in open(D + "foldseek/AFDB50.m8"):
    f = l.rstrip("\n").split("\t"); e = float(f[6])
    if e <= 1e-3: byu[unit(f[0])].append((e, float(f[7]), f[1]))
want = set()
for v in byu.values():
    for e, b, t in sorted(v, key=lambda h: (h[0], -h[1]))[:25]:
        if acc(t): want.add(acc(t))
live = set()
for f in glob.glob(D + "uniprot/batch_*.tsv"): live |= {l.split("\t")[0] for l in list(open(f))[1:]}
gap = {l.strip() for l in list(open(D + "uniprot_gap/missing_accessions.tsv"))[1:]}
uns = set()
for l in open(D + "uniprot_unisave/unisave_raw.jsonl"):
    try: uns.add(json.loads(l)["acc"])
    except ValueError: pass
need = sorted(want - live - uns)
print("b4 requested %d; live %d; in gap list %d; gap-list accessions b4 does not request %d; to fetch %d" % (len(want), len(want & live), len(want & gap), len(gap - want), len(need)))
with open(OUT, "w") as o:
    o.write("accession\n")
    for a in need: o.write(a + "\n")
print("wrote", OUT)
