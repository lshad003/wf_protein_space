# Step 36b B: Pfam near misses among the 4,509 recurrent EU families.
# Search: HMMER 3.4 hmmsearch --max -E 10 --domE 10 vs Pfam-A 38.2 (scripts/step36b_hmm.sh) of every distinct member
#   sequence (results/step36b/members.tsv, uniq_map.tsv): family reps, all members of the families' 95% clusters in the
#   88 metagenomes, strict flagged. Scores from the tblout: full-sequence bit score and E, best 1-domain bit score.
# GA: the model's GA line in Pfam-A.hmm (sequence GA, domain GA). A sequence passes GA as under --cut_ga: full score
#   >= sequence GA and best domain >= domain GA. Gap of a failing sequence = max(sequence GA - full score,
#   domain GA - best domain); gap_full and gap_domain are reported separately.
# Near miss (criterion): the family representative fails GA on a Pfam model whose name does not start DUF or UPF, with
#   gap <= 3 bits. One row per family and model. Members passing GA: among all members in the 88 and among strict ones.
# Structure (results/step34j_structures_v2.tsv; member-validated families only, else not_modelled): hit_class,
#   best_target / best_target_name, structure_confident_annotated, confident_target / confident_target_name.
# Context: Step 36 conserved labels of the family (results/step36/context_by_family.tsv, EU_family30).
# Strong candidate (proposed, for manual review): (a) >= 1 member passes GA on the model, or (b) the structure hit text
#   or a conserved neighbour names the same family or system as the model: same Step 34j keyword theme
#   (scripts/step34j_themes.py) on the model name + description vs the hit name or the neighbour's Pfam name +
#   description, or a shared specific word (lowercased words of >= 4 letters outside a generic list).
# Positive control: Cse2 (UHM45.35758__k141_100597_236, CRISPR_Cse2) must be a near miss with members passing GA;
#   otherwise no output and exit 2.
# Outputs (refuse to overwrite): results/step36b/near_miss_pfam.tsv, results/step36b/near_miss_strong.tsv,
#   results/step36b/rep_best_pfam.tsv (best-scoring non-DUF model per family rep, any gap).
import os, sys, re, glob, collections
sys.path.insert(0, "/bigdata/stajichlab/lshad003/wf_protein_space/scripts"); import step34j_themes as TH
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; D = R + "step36b/"
PF = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/"
ON, OS, OB = D + "near_miss_pfam.tsv", D + "near_miss_strong.tsv", D + "rep_best_pfam.tsv"
for p in (ON, OS, OB):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
GAP = 3.0; CSE2 = "UHM45.35758__k141_100597_236"
def tsv(p):
    with open(p) as fh:
        h = next(fh).rstrip("\n").split("\t")
        for l in fh: yield dict(zip(h, l.rstrip("\n").split("\t")))
# Pfam models: GA and descriptions
ga, acc = {}, {}; nm = None
for l in open(PF + "Pfam-A.hmm"):
    if l.startswith("NAME"): nm = l.split()[1]
    elif l.startswith("ACC"): acc[nm] = l.split()[1]
    elif l.startswith("GA"): x = l.split(); ga[nm] = (float(x[1]), float(x[2].rstrip(";")))
de = {}; cur = None
for l in open(PF + "Pfam-A.hmm.dat", errors="replace"):
    if l.startswith("#=GF ID"): cur = l[10:].strip()
    elif l.startswith("#=GF DE") and cur: de[cur] = l[10:].strip()
print("Pfam models %d, with GA %d, with DE %d" % (len(acc), len(ga), len(de)), flush=True)
# members
M = list(tsv(D + "members.tsv")); u = {}
for r in tsv(D + "uniq_map.tsv"): u[r["member"]] = r["uid"]
fam_mem = collections.defaultdict(list)
for r in M: fam_mem[r["family_rep"]].append(r)
fams = sorted(fam_mem); print("families %d, member rows %d, uids %d" % (len(fams), len(M), len(set(u.values()))), flush=True)
# hits: uid -> model -> (full, E, bestdom)
need = set(u.values()); H = collections.defaultdict(dict); nrow = 0
ch = sorted(glob.glob(D + "hmm/C.*.tblout"))
if len(ch) != len(glob.glob(D + "chunks/C.*.faa")): raise SystemExit("ABORT: %d tblout for %d chunks" % (len(ch), len(glob.glob(D + "chunks/C.*.faa"))))
for f in ch:
    for l in open(f):
        if l[0] == "#": continue
        p = l.split(); nrow += 1
        H[p[0]][p[2]] = (float(p[5]), float(p[4]), float(p[8]))
print("hit rows %d (E <= 10), sequences with a hit %d of %d" % (nrow, len(H), len(need)), flush=True)
def passes(h, m): return h[0] >= ga[m][0] and h[2] >= ga[m][1]
def gaps(h, m): return ga[m][0] - h[0], ga[m][1] - h[2]
isduf = lambda m: m.startswith("DUF") or m.startswith("UPF")
# structure, context, member-validated
st = {r["family_rep"]: r for r in tsv(R + "step34j_structures_v2.tsv") if r["set"] == "EU"}
mc = collections.defaultdict(list)
for r in tsv(R + "step34h_member_classes.tsv"):
    if r["family_rep"] in fam_mem: mc[r["family_rep"]].append(r["class_any"])
mval = {f: int(bool(mc[f]) and all(c == "EU" for c in mc[f])) for f in fams}
cx = {r["unit"]: r for r in tsv(R + "step36/context_by_family.tsv") if r["run"] == "EU_family30"}
GEN = set("""domain protein family containing putative like terminal region subunit uncharacterized uncharacterised structure
crystal from with resolution bacterial type associated predicted probable related conserved hypothetical function unknown
binding repeat fold superfamily system complex factor component chain small large alpha beta helix helical motif
cryo putative possible domains proteins homolog homologue sequence specific dependent""".split())
def words(t): return {w for w in re.findall(r"[a-z][a-z0-9]{3,}", (t or "").lower()) if w not in GEN}
def match(model, texts):
    mt = model + " " + de.get(model, ""); mth, mw = TH.themes(mt), words(mt.replace("_", " "))
    ev = []
    for src, t in texts:
        if not t or t in ("-", "NA"): continue
        th = mth & TH.themes(t); w = mw & words(t.replace("_", " "))
        if th: ev.append("%s theme %s (%s)" % (src, "/".join(sorted(th)), t[:60]))
        elif w: ev.append("%s word %s (%s)" % (src, "/".join(sorted(w)), t[:60]))
    return ev
rows, best = [], []
for f in fams:
    rid = u.get(f)
    if rid is None: raise SystemExit("ABORT: family rep not searched " + f)
    rh = H.get(rid, {}); allm = [r for r in fam_mem[f] if r["strict"] != "rep_not_in_88"]; strict = [r for r in allm if r["strict"] == "1"]
    cand = []
    for m, h in rh.items():
        if isduf(m) or m not in ga: continue
        if passes(h, m): cand.append((m, h, 0.0, 0.0, 0.0, True)); continue
        gf, gd = gaps(h, m); cand.append((m, h, max(gf, gd), gf, gd, False))
    if cand:
        b = min(cand, key=lambda c: (c[2], -c[1][0])); best.append((f, b[0], b[1][0], b[1][2], ga[b[0]][0], ga[b[0]][1], b[2], int(b[5])))
    s = st.get(f, {}); c = cx.get(f, {})
    for m, h, g, gf, gd, ok in cand:
        if ok or g > GAP: continue
        pa = sum(1 for r in allm if u[r["member"]] in H and m in H[u[r["member"]]] and passes(H[u[r["member"]]][m], m))
        ps = sum(1 for r in strict if u[r["member"]] in H and m in H[u[r["member"]]] and passes(H[u[r["member"]]][m], m))
        hit_txt = s.get("confident_target_name") if s.get("structure_confident_annotated") == "1" else s.get("best_target_name")
        labs = [x.split(":")[1] for x in c.get("conserved_labels", "-").split(";") if x.startswith("Pfam:") or x.startswith("eggNOG:")]
        texts = [("structure", hit_txt)] + [("neighbour", "%s %s" % (x, de.get(x, ""))) for x in labs]
        ev = match(m, texts)
        rows.append(dict(family_rep=f, model=m, pfam_acc=acc.get(m, "NA"), model_description=de.get(m, ""), rep_full_bits=h[0], rep_best_domain_bits=h[2],
            rep_full_evalue=h[1], ga_seq=ga[m][0], ga_dom=ga[m][1], gap=round(g, 2), gap_full=round(gf, 2), gap_domain=round(gd, 2),
            n_members_pass_ga=pa, n_members=len(allm), n_strict_pass_ga=ps, n_strict=len(strict), member_validated=mval[f],
            hit_class=s.get("hit_class", "not_modelled"), structure_confident=s.get("structure_confident_annotated", "NA"),
            best_target=s.get("best_target", "NA"), best_target_name=s.get("best_target_name", "NA"),
            confident_target=s.get("confident_target", "NA"), confident_target_name=s.get("confident_target_name", "NA"),
            step36_context=c.get("context_class", "NA"), step36_conserved_labels=c.get("conserved_labels", "NA"),
            strong_a_members_pass=int(pa > 0), strong_b_named=int(bool(ev)), match_evidence="; ".join(ev) or "-"))
# positive control
cs = [r for r in rows if r["family_rep"] == CSE2 and r["model"] == "CRISPR_Cse2"]
if not cs or not cs[0]["strong_a_members_pass"]:
    print("POSITIVE CONTROL FAILED:", cs); raise SystemExit(2)
c = cs[0]
print("Cse2 recovered: rep full %.1f (domain %.1f), GA %.1f/%.1f, gap %.2f (full %.2f, domain %.2f); members passing GA %d of %d, strict %d of %d; evidence: %s"
      % (c["rep_full_bits"], c["rep_best_domain_bits"], c["ga_seq"], c["ga_dom"], c["gap"], c["gap_full"], c["gap_domain"],
         c["n_members_pass_ga"], c["n_members"], c["n_strict_pass_ga"], c["n_strict"], c["match_evidence"]))
rows.sort(key=lambda r: (-(r["strong_a_members_pass"] + r["strong_b_named"]), -r["n_members_pass_ga"], r["gap"], r["family_rep"], r["model"]))
K = list(rows[0].keys())
def w(p, R_):
    with open(p, "w") as o:
        o.write("\t".join(K) + "\n")
        for r in R_: o.write("\t".join(str(r[k]).replace("\t", " ") for k in K) + "\n")
strong = [r for r in rows if r["strong_a_members_pass"] or r["strong_b_named"]]
w(ON, rows); w(OS, strong)
with open(OB, "w") as o:
    o.write("family_rep\tmodel\trep_full_bits\trep_best_domain_bits\tga_seq\tga_dom\tgap\trep_passes_ga\n")
    for x in best: o.write("\t".join(str(v) for v in x) + "\n")
nf = lambda R_: len({r["family_rep"] for r in R_})
print("\nfamily reps passing GA on a non-DUF model (should be 0 for EU reps): %d" % sum(1 for x in best if x[7]))
print("near misses (gap <= %.0f bits, non-DUF/UPF): %d family-model pairs in %d families (member-validated %d)"
      % (GAP, len(rows), nf(rows), nf([r for r in rows if r["member_validated"] == 1])))
for lab, R_ in (("strong, any", strong), ("(a) members pass GA", [r for r in rows if r["strong_a_members_pass"]]),
                ("(b) structure or neighbour names the model", [r for r in rows if r["strong_b_named"]]),
                ("(a) and (b)", [r for r in rows if r["strong_a_members_pass"] and r["strong_b_named"]])):
    print("  %-44s %4d pairs, %4d families" % (lab, len(R_), nf(R_)))
print("\ntop strong candidates:")
for r in strong[:15]:
    print("  %s %s gap %.2f (full %.1f / GA %.1f) members pass %d/%d (strict %d/%d) mv %d | %s %s | %s | %s"
          % (r["family_rep"], r["model"], r["gap"], r["rep_full_bits"], r["ga_seq"], r["n_members_pass_ga"], r["n_members"], r["n_strict_pass_ga"],
             r["n_strict"], r["member_validated"], r["hit_class"], (r["confident_target_name"] if r["structure_confident"] == "1" else r["best_target_name"])[:50],
             r["step36_conserved_labels"][:80], r["match_evidence"][:120]))
print("wrote", ON, OS, OB)
