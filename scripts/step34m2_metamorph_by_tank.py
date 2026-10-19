# Step 34m2: metamorphosis date per rearing tank for the two within-clutch tank contrasts of Step 34m
#   (WF24 egg mass 4: tanks 4.1, 4.2; WF23 egg mass 3: tanks 5.0, 6.0).
# Input: metadata/wf23_wf24_sequenced_metadata.tsv (column 7, metamorphosis_date, ISO dates); excluded animals from
#   metadata/wf24_excluded.tsv (UHM586.41010, UHM590.41012). Median of an even n is the midpoint of the two middle dates.
# Output: results/step34m2_metamorph_by_tank.tsv (n, min, median, max, distinct dates); refuses to overwrite.
import os, datetime
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
OUT = W + "/results/step34m2_metamorph_by_tank.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
ex = {l.split("\t")[0] for l in list(open(W + "/metadata/wf24_excluded.tsv"))[1:]}
if ex != {"UHM586.41010", "UHM590.41012"}: raise SystemExit("ABORT: unexpected exclusions %s" % ex)
WANT = {("WF24", "4"): ("4.1", "4.2"), ("WF23", "3"): ("5.0", "6.0")}
g = {}
with open(W + "/metadata/wf23_wf24_sequenced_metadata.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    if h[6] != "metamorphosis_date": raise SystemExit("ABORT: column 7 is %s" % h[6])
    for l in fh:
        f = dict(zip(h, l.rstrip("\n").split("\t")))
        if f["stem"] in ex or (f["wf_year"], f["egg_mass"]) not in WANT: continue
        g.setdefault((f["wf_year"], f["egg_mass"], f["tank"]), []).append(datetime.date.fromisoformat(f["metamorphosis_date"]))
rows = []
for (y, e), tanks in WANT.items():
    for t in tanks:
        v = sorted(g.get((y, e, t), []))
        if not v: raise SystemExit("ABORT: no animals for %s EM%s tank %s" % (y, e, t))
        n = len(v); med = v[n // 2] if n % 2 else v[n // 2 - 1] + (v[n // 2] - v[n // 2 - 1]) / 2
        rows.append((y, e, t, n, v[0], med, v[-1], len(set(v))))
with open(OUT, "w") as o:
    o.write("wf_year\tegg_mass\ttank\tn_animals\tmin_date\tmedian_date\tmax_date\tdistinct_dates\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
for r in rows: print("\t".join(map(str, r)))
print("wrote", OUT)
