# Checks the Section 3.6 sentence on month-1 sequencing runs: WF22 month x egg mass x run.
import csv, re, sys, collections
B = "/bigdata/stajichlab/lshad003/wf_protein_space/results/batch_table_89.tsv"
D = "/bigdata/stajichlab/lshad003/wf_protein_space/metadata/wf22_design.tsv"

def load(p):
    with open(p) as f:
        r = csv.DictReader(f, delimiter="\t")
        return r.fieldnames, list(r)

def stemcol(cols, rows):
    for c in cols:
        if re.match(r"^UHM\d+\.\d+$", (rows[0].get(c) or "").strip()):
            return c
    return cols[0]

bc, br = load(B)
dc, dr = load(D)
print("batch columns:", bc)
print("design columns:", dc)
bs, ds = stemcol(bc, br), stemcol(dc, dr)
runs = [c for c in bc if "run" in c.lower()]
mon = [c for c in dc if "month" in c.lower()]
egg = [c for c in dc if "egg" in c.lower()]
print("stem:", bs, "/", ds, "| run:", runs, "| month:", mon, "| egg:", egg)
if not (runs and mon and egg):
    sys.exit("STOP: run, month or egg column not found; paste the lines above")
run, m, e = runs[0], mon[0], egg[0]
bidx = dict((row[bs].strip(), row) for row in br)
miss = [row[ds] for row in dr if row[ds].strip() not in bidx]
print("design rows:", len(dr), "| not in batch table:", len(miss), miss)
t1, t2 = collections.Counter(), collections.Counter()
for row in dr:
    b = bidx.get(row[ds].strip())
    r = b[run] if b else "NA"
    t1[(row[m], r)] += 1
    t2[(row[m], row[e], r)] += 1
print("\nmonth\trun\tn")
for k in sorted(t1):
    print(k[0], k[1], t1[k], sep="\t")
print("\nmonth\tegg_mass\trun\tn")
for k in sorted(t2):
    print(k[0], k[1], k[2], t2[k], sep="\t")
