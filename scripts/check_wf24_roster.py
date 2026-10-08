import csv, collections, glob
BASE = "/bigdata/stajichlab/lshad003/wf_protein_space/"
roster = {
 "1": "UHM583 UHM584 UHM608 UHM610 UHM648", "2": "UHM581 UHM585 UHM615 UHM624 UHM632",
 "3": "UHM577 UHM578 UHM591 UHM596 UHM636", "4": "UHM586 UHM598 UHM617 UHM634 UHM649",
 "5": "UHM582 UHM587 UHM616 UHM626 UHM633", "6": "UHM575 UHM594 UHM607 UHM644 UHM666",
 "7": "UHM574 UHM580 UHM611 UHM628 UHM656"}
jason = {a: c for c, s in roster.items() for a in s.split()}
def load(f):
    return {r["stem"]: r for r in csv.DictReader(open(BASE + "metadata/" + f), delimiter="\t")}
main = load("samples_89_treatment.tsv")
forj = load("samples_89_for_jason.tsv")
print("for_jason columns:", list(next(iter(forj.values())).keys()))
wf24 = [s for s, r in main.items() if r["cohort"] == "WF24"]
print("WF24 rows:", len(wf24), "| missing from for_jason:", [s for s in wf24 if s not in forj])
def info(s):
    m, j = main[s], forj.get(s, {})
    return "treatment=%s resolved=%s egg_mass_cohort=%s jason_file_treatment=%s" % (
        m["treatment"], j.get("treatment_resolved", "NA"), m["egg_mass_cohort"], j.get("treatment", "NA"))
print("\n== in catalog, NOT on roster ==")
for s in wf24:
    if s.split(".")[0] not in jason: print(s, info(s))
print("\n== on roster, NOT in catalog ==")
have = {s.split(".")[0] for s in wf24}
for a, c in sorted(jason.items()):
    if a not in have: print(a, "code", c)
print("\n== code mismatches (only where a recorded value is a code 1-7) ==")
for s in wf24:
    a = s.split(".")[0]
    codes = [v for v in (main[s]["treatment"], forj.get(s, {}).get("treatment", "")) if v.strip() in roster]
    for v in codes:
        if a in jason and v.strip() != jason[a]: print(s, "recorded=" + v, "jason=" + jason[a])
print("\n== code to resolved strain, per Jason's roster ==")
pairs = collections.Counter((jason[s.split(".")[0]], forj.get(s, {}).get("treatment_resolved", "NA"))
                            for s in wf24 if s.split(".")[0] in jason)
for k, n in sorted(pairs.items()): print(k, n)
print("\n== files containing 7724 ==")
for f in sorted(glob.glob(BASE + "metadata/*")):
    try:
        for i, line in enumerate(open(f), 1):
            if "7724" in line: print(f, i, line.rstrip()[:120])
    except Exception as e:
        print("skip", f, e)
