# Step 34c5: NMPFams2 match status of the recurrent EU families at the family level (2026-10-23).
# The 3,868 of 4,127 (93.72%) of Step 34c2c and Fig. 2f counts searched 95% representatives, not families.
# Search set: complete ORFs >= 35 aa (in_step26_set == 1). Two family-level rules:
#   any_searched: family unmatched if none of its searched representatives matched (families with >= 1 searched rep)
#   family_rep:   family unmatched if its family representative was searched and did not match
# Inputs: results/step34c_strict_eu_nmpfams2.tsv, step34c2c_strict_v2_families.tsv (recurrent, not excluded),
# step34b2_nmpfams2_filter.tsv. Output: results/step34c5_family_nmpfams2.tsv (one row per recurrent family).
import csv, os
R = "/bigdata/stajichlab/lshad003/wf_protein_space/results/"
OUT = R + "step34c5_family_nmpfams2.tsv"
if os.path.exists(OUT): raise SystemExit("ABORT: %s exists" % OUT)
def rows(f): return list(csv.DictReader(open(R + f), delimiter="\t"))
fam = [r["family_rep"] for r in rows("step34c2c_strict_v2_families.tsv") if r["recurrent"] == "1" and r["excluded"] == "0"]
fset = set(fam); print("recurrent families:", len(fam))
flt = {r["rep"]: r for r in rows("step34b2_nmpfams2_filter.tsv")}
nm = [r for r in rows("step34c_strict_eu_nmpfams2.tsv") if r["family_rep"] in fset]
by = {}
for r in nm: by.setdefault(r["family_rep"], []).append(r)
srch = [r for r in nm if r["in_step26_set"] == "1"]
print("representatives in recurrent families %d | searched %d | unmatched %d (%.2f%%)"
      % (len(nm), len(srch), sum(r["matched"] == "0" for r in srch), 100.0 * sum(r["matched"] == "0" for r in srch) / len(srch)))
n = {"any": 0, "any_u": 0, "own": 0, "own_u": 0}
with open(OUT, "w") as o:
    o.write("family_rep\tn_reps\tn_reps_searched\tn_reps_matched\tfamily_rep_F1_complete\tfamily_rep_F2_ge35aa\t"
            "family_rep_searched\tfamily_rep_matched\tunmatched_any_searched\tunmatched_family_rep\n")
    for f in fam:
        rs = by.get(f, []); s = [r for r in rs if r["in_step26_set"] == "1"]; m = sum(r["matched"] == "1" for r in s)
        own = [r for r in rs if r["rep"] == f]
        if len(own) != 1: raise SystemExit("ABORT: family rep %s found %d times" % (f, len(own)))
        own = own[0]; osr = own["in_step26_set"] == "1"
        if osr != (flt[f]["F1_complete"] == "1" and flt[f]["F2_ge35aa"] == "1"): raise SystemExit("ABORT: filter mismatch " + f)
        ua = ("1" if m == 0 else "0") if s else "NA"; uo = ("1" if own["matched"] == "0" else "0") if osr else "NA"
        if s: n["any"] += 1; n["any_u"] += ua == "1"
        if osr: n["own"] += 1; n["own_u"] += uo == "1"
        o.write("\t".join(map(str, [f, len(rs), len(s), m, flt[f]["F1_complete"], flt[f]["F2_ge35aa"],
                                    int(osr), own["matched"] if osr else "NA", ua, uo])) + "\n")
print("family level, any searched representative: unmatched %d of %d (%.2f%%)" % (n["any_u"], n["any"], 100.0 * n["any_u"] / n["any"]))
print("family level, family representative only:  unmatched %d of %d (%.2f%%)" % (n["own_u"], n["own"], 100.0 * n["own_u"] / n["own"]))
ok = (n["any_u"], n["any"], n["own_u"], n["own"]) == (3501, 3668, 3484, 3646)
print("reproduces 3,501/3,668 and 3,484/3,646:", "yes" if ok else "NO")
print("wrote", OUT)
if not ok: raise SystemExit(1)
