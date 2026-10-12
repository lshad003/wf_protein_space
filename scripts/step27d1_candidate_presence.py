# Step 27d1: sample presence of the 61 full-length Basidiobolus candidate targets (step27c).
# Targets: results/basid_candidates_context.tsv (col 1, 95% reps; qcov >= 0.8 and tcov >= 0.8).
# Members: results/clusters_95.tsv (rep, member; all 89 samples). Sample = member stem (before "__").
# Cohort/treatment: metadata/samples_88_treatment.tsv, WF24 codes via metadata/wf24_treatment_key.tsv (4 = Control);
#   stems only in metadata/wf24_excluded.tsv (UHM586.41010, recorded code 4 control, contaminated, excluded from
#   n88) get treatment "excluded_control" and are counted separately from the 88-sample controls.
# Year from cohort name: WF22 2022, WF23 2023, WF24 2024.
# Presence here means a gene in the 95% cluster was predicted from that sample's assembly; not read-level presence.
# Read-only on every _n88 file.
import os, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results"
M = W + "/metadata"
OUT = R + "/basid_candidates_presence.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")

ctx = {}
for i, l in enumerate(open(R + "/basid_candidates_context.tsv")):
    if i == 0: H = {c: j for j, c in enumerate(l.rstrip("\n").split("\t"))}; continue
    f = l.rstrip("\n").split("\t"); ctx[f[0]] = f
print("targets:", len(ctx))

key = {}
for i, l in enumerate(open(M + "/wf24_treatment_key.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t"); key[f[0]] = f[1]
meta = {}
for i, l in enumerate(open(M + "/samples_88_treatment.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    t = key[f[2]] if f[1] == "WF24" else f[2]
    meta[f[0]] = (f[1], t, "control" if t.lower() == "control" else "inoculated")
for i, l in enumerate(open(M + "/wf24_excluded.tsv")):
    if i == 0: continue
    f = l.rstrip("\n").split("\t")
    meta.setdefault(f[0], ("WF24", "excluded_control", "excluded"))
YR = {"WF22": "2022", "WF23": "2023", "WF24": "2024"}

mem = collections.defaultdict(list)
for l in open(R + "/clusters_95.tsv"):
    r, m = l.rstrip("\n").split("\t")[:2]
    if r in ctx: mem[r].append(m)
miss = [t for t in ctx if t not in mem]
if miss: raise SystemExit("ABORT: targets not found as reps in clusters_95.tsv: %d" % len(miss))

HDR = ["target", "context", "rep_stem", "rep_treatment", "n_queries", "n_members", "n_samples", "samples",
       "cohorts", "years", "n_control_samples", "control_samples", "n_excluded_samples", "n_inoculated_samples",
       "in_any_control", "treatments", "members"]
rows = []
for t, ms in mem.items():
    st = sorted({m.split("__", 1)[0] for m in ms})
    for s in st:
        if s not in meta: raise SystemExit("ABORT: stem without metadata: " + s)
    ctl = [s for s in st if meta[s][2] == "control"]
    exc = [s for s in st if meta[s][2] == "excluded"]
    ino = [s for s in st if meta[s][2] == "inoculated"]
    c = ctx[t]
    rows.append([t, c[H["context"]], c[H["stem"]], c[H["treatment"]], c[H["n_queries"]], str(len(ms)), str(len(st)),
                 ",".join(st), ",".join(sorted({meta[s][0] for s in st})),
                 ",".join(sorted({YR[meta[s][0]] for s in st})), str(len(ctl)), ",".join(ctl) or "NA",
                 str(len(exc)), str(len(ino)), "yes" if ctl else "no",
                 ",".join("%s:%s" % (s, meta[s][1]) for s in st), ",".join(sorted(ms))])
rows.sort(key=lambda r: (-int(r[6]), r[0]))
with open(OUT, "w") as fh:
    fh.write("\t".join(HDR) + "\n")
    for r in rows: fh.write("\t".join(r) + "\n")

C = {c: i for i, c in enumerate(HDR)}
print("\ntarget\tcontext\tn_samples\tyears\tn_control\tn_excluded\tn_inoculated\tcontrol_samples")
for r in rows:
    print("\t".join(r[C[c]] for c in ("target", "context", "n_samples", "years", "n_control_samples",
                                      "n_excluded_samples", "n_inoculated_samples", "control_samples")))
nc = sum(r[C["in_any_control"]] == "yes" for r in rows)
print("\ntargets in any 88-sample control:", nc, "of", len(rows))
print("targets in the excluded control UHM586.41010:", sum(r[C["n_excluded_samples"]] != "0" for r in rows))
print("targets in any control incl. excluded:",
      sum(r[C["in_any_control"]] == "yes" or r[C["n_excluded_samples"]] != "0" for r in rows))
print("samples per target:", dict(sorted(collections.Counter(int(r[C["n_samples"]]) for r in rows).items())))
print("targets by number of years:", dict(collections.Counter(len(r[C["years"]].split(",")) for r in rows)))
print("\nwrote", OUT)
