# Step 34c2: genus labels of the strict EU 30% families (descriptive only).
# Units: recurrent strict EU families (>= 10 animals AND >= 2 cohorts; results/step34c_strict_eu_families.tsv), the
#   non-recurrent strict EU families as baseline (a), and a random 5,000 strict K reps (seed 34) as baseline (b).
# Members: every member protein of the unit's 95% reps (results/clusters_95_n88.tsv); family -> reps from
#   results/step34c_strict_eu_clusters30.tsv. Strict = class4 and origin Tiara-prokaryotic in results/step34a_rep_origin.tsv.
# Contig taxonomy: results/gene_taxonomy_n88.tsv holds one row per supported rep only, so member contigs are read from
#   the same per-stem source step20c used (results_scaffold_classify_mmseqs/<stem>/<stem>_uniref50_lca.tsv[.gz]:
#   contig, taxid, rank, name). Cross-check: every rep's label must equal its gene_taxonomy_n88.tsv row.
# Genus as in step30b: name when rank is genus; first word of the name for ranks below genus ("Candidatus X" kept as
#   two words; "uncultured" gives none); otherwise none.
# Majority genus = most common genus among the unit's genus-labelled members (ties broken by name); its share is given.
# A unit is in a flagged genus when its majority genus is flagged; "any flagged" = at least one member in a flagged genus.
# Output: results/step34c2_recurrent_genus.tsv (one row per unit, all three sets); refuses to overwrite.
import os, gzip, random, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"
R = W + "/results/"
FEC = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/results_scaffold_classify_mmseqs"
OUT = R + "step34c2_recurrent_genus.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
FLAG = ["Acetobacter", "Gluconobacter", "Komagataeibacter", "Commensalibacter", "Lactobacillus", "Lactiplantibacillus",
        "Levilactobacillus", "Fructilactobacillus", "Leuconostoc", "Wolbachia", "Providencia"]
FS = set(FLAG)
BELOW = {"species", "strain", "species group", "species subgroup", "subspecies", "subgenus", "forma", "varietas",
         "isolate", "serotype", "serogroup", "biotype", "genotype", "forma specialis", "morph", "pathogroup"}
def genus(rank, name):
    if rank == "genus": return name
    if rank in BELOW:
        w = name.split()
        if not w or w[0] == "uncultured": return None
        return " ".join(w[:2]) if w[0] == "Candidatus" and len(w) > 1 else w[0]
    return None
print("genus rule (step30b): rank genus -> name; rank below genus -> first word of name (Candidatus X kept, uncultured none); else none")

sK = []; sEU = set()
with open(R + "step34a_rep_origin.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[2] != "Tiara-prokaryotic": continue
        if f[1] == "EU": sEU.add(f[0])
        elif f[1] == "K": sK.append(f[0])
print("strict EU reps %d | strict K reps %d" % (len(sEU), len(sK)), flush=True)
if len(sEU) != 41392 or len(sK) != 1120953: raise SystemExit("ABORT: strict counts differ from Step 34a")
rec = {}
with open(R + "step34c_strict_eu_families.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); rec[f[0]] = f[8] == "1"
reps_of = collections.defaultdict(list)
for l in open(R + "step34c_strict_eu_clusters30.tsv"):
    a, b = l.rstrip("\n").split("\t"); reps_of[a].append(b)
print("families %d | recurrent %d" % (len(rec), sum(rec.values())), flush=True)
if sum(rec.values()) != 4561 or set(reps_of) != set(rec): raise SystemExit("ABORT: family files disagree with Step 34c")
random.seed(34)
kpick = sorted(random.sample(sorted(sK), 5000))
units = []  # (set, unit id, reps)
for f in sorted(rec): units.append(("recurrent_EU_family" if rec[f] else "nonrecurrent_EU_family", f, reps_of[f]))
for r in kpick: units.append(("random_K_rep", r, [r]))
need = {r for u in units for r in u[2]}
mem = collections.defaultdict(list)
with open(R + "clusters_95_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if r in need: mem[r].append(m)
if set(mem) != need: raise SystemExit("ABORT: reps without members")
def contig(g):
    s, rest = g.split("__", 1); return s, rest.rsplit("_", 1)[0]
want = collections.defaultdict(set)
for ms in mem.values():
    for m in ms:
        s, c = contig(m); want[s].add(c)
for r in need:
    s, c = contig(r); want[s].add(c)
print("member proteins %d | stems %d | contigs to label %d" % (sum(len(v) for v in mem.values()), len(want),
      sum(len(v) for v in want.values())), flush=True)
lab = {}; nofile = []
for s in sorted(want):
    f = "%s/%s/%s_uniref50_lca.tsv" % (FEC, s, s)
    if not os.path.exists(f): f += ".gz"
    if not os.path.exists(f): nofile.append(s); continue
    w = want[s]
    for l in (gzip.open(f, "rt") if f.endswith(".gz") else open(f)):
        p = l.rstrip("\n").split("\t")
        if len(p) < 10 or p[0] not in w: continue
        lab[s + "__" + p[0]] = (p[2], p[3])
print("stems without an LCA file:", nofile, "| contigs with an LCA row %d" % len(lab), flush=True)
gt = {}
with open(R + "gene_taxonomy_n88.tsv") as fh:
    next(fh)
    for l in fh:
        p = l.rstrip("\n").split("\t")
        if p[0] in need: gt[p[0]] = (p[2], p[3])
bad = 0; both = 0
for r in need:
    s, c = contig(r); t = lab.get(s + "__" + c)
    if (t is None) != (r not in gt): bad += 1
    elif t is not None:
        both += 1; bad += t != gt[r]
print("cross-check vs gene_taxonomy_n88.tsv on %d reps: %d with a row, %d disagreements" % (len(need), both, bad), flush=True)
if bad: raise SystemExit("ABORT: member-contig labels disagree with gene_taxonomy_n88.tsv")
def gen_of(m):
    s, c = contig(m); t = lab.get(s + "__" + c)
    return genus(*t) if t else None
rows = []
for st, u, rs in units:
    ms = [m for r in rs for m in mem[r]]
    withrow = sum(1 for m in ms if lab.get("__".join(contig(m))) is not None)
    g = collections.Counter(x for x in (gen_of(m) for m in ms) if x)
    ng = sum(g.values())
    if ng:
        top = sorted(g.items(), key=lambda x: (-x[1], x[0]))[0]
        maj, share = top[0], top[1] / ng
    else: maj, share = "NA", 0.0
    anyf = any(x in FS for x in g)
    rows.append((st, u, len(rs), len(ms), withrow, ng, len(g), maj, share, int(maj in FS), int(anyf)))
with open(OUT, "w") as fh:
    fh.write("set\tunit\tn_reps\tn_members\tmembers_with_lca_row\tmembers_with_genus\tn_genera\tmajority_genus\tmajority_share\tmajority_flagged\tany_member_flagged\n")
    for r in rows: fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f\t%d\t%d\n" % r)
SETS = ("recurrent_EU_family", "nonrecurrent_EU_family", "random_K_rep")
print("\nset                       units  >=1 member LCA row   >=1 member genus   majority flagged (of genus-labelled)   any member flagged   majority share >= 0.5")
for st in SETS:
    X = [r for r in rows if r[0] == st]; L = [r for r in X if r[5] > 0]
    mf = sum(r[9] for r in L); af = sum(r[10] for r in L); half = sum(1 for r in L if r[8] >= 0.5)
    print("%-24s %6d  %6d (%5.1f%%)   %6d (%5.1f%%)   %6d (%5.2f%%)                     %6d (%5.2f%%)    %6d (%5.1f%%)" % (
        st, len(X), sum(1 for r in X if r[4] > 0), 100.0 * sum(1 for r in X if r[4] > 0) / len(X), len(L), 100.0 * len(L) / len(X),
        mf, 100.0 * mf / len(L), af, 100.0 * af / len(L), half, 100.0 * half / len(L)))
print("\ntop 25 majority genera, recurrent families (n with a genus label); share in baselines (of genus-labelled units)")
cnt = {st: collections.Counter(r[7] for r in rows if r[0] == st and r[5] > 0) for st in SETS}
den = {st: sum(cnt[st].values()) for st in SETS}
print("%-4s %-28s %7s %7s   %14s   %14s" % ("rank", "genus", "rec", "rec%", "nonrec%", "K5000%"))
for i, (g, n) in enumerate(cnt[SETS[0]].most_common(25), 1):
    print("%-4d %-28s %7d %6.2f%%   %13.2f%%   %13.2f%%%s" % (i, g, n, 100.0 * n / den[SETS[0]], 100.0 * cnt[SETS[1]][g] / den[SETS[1]],
          100.0 * cnt[SETS[2]][g] / den[SETS[2]], "   FLAGGED" if g in FS else ""))
print("\nflagged genera as majority genus (units): recurrent / non-recurrent / K5000")
for g in FLAG:
    print("  %-22s %6d %6d %6d" % (g, cnt[SETS[0]][g], cnt[SETS[1]][g], cnt[SETS[2]][g]))
print("wrote", OUT, "rows:", len(rows))
