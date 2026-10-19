# Step 34j A2: non-redundant genomic contexts and context conservation, EU recurrent families and known controls.
# Input: results/step34j/context_members.tsv, context_neighbours.tsv, rep_annot.tsv (step34j_a1_gather.py).
# Context of a member = its +-5 neighbour genes (oriented to the member). Non-redundant (NR) contexts: two contexts share
#   a background if their neighbour 95%-rep lists are identical or share >= 80% (multiset overlap / larger list size);
#   greedy, contexts sorted by number of neighbours (desc) then member ID, each joins the first seed it matches.
#   Members on single-gene contigs (no neighbour) have no context and are counted separately. Each NR context is
#   represented by its seed; all scores are over NR contexts only.
# Neighbour units: 50% cluster of the neighbour gene (results/clusters_50_n88.tsv; neighbours in the member's own 50%
#   cluster are excluded as tandem copies), Pfam accession, KEGG KO and KEGG module of the neighbour's 95% rep.
#   Fraction of NR contexts carrying a unit anywhere within +-5.
# Score: conservation = fraction of NR contexts carrying the top 50%-cluster neighbour (top neighbour). Strand
#   conservation = fraction of those contexts with the modal relative strand; distance conservation = fraction at the
#   modal oriented offset (closest copy per context); median gap in bp.
# Percentile: mid-rank percentile of a family's conservation among controls in the same NR-context bin
#   (1, 2, 3-4, 5-9, >= 10). Conserved context (criterion): >= 3 NR contexts AND conservation >= 0.5.
# Outputs (refuse to overwrite): results/step34j/context_scores.tsv (EU families and controls),
#   results/step34j/context_nr_seeds.tsv (set, family, seed member, contexts collapsed into it).
import os, collections
D = "/bigdata/stajichlab/lshad003/wf_protein_space/results/step34j/"; O = D + "context_scores.tsv"; OSD = D + "context_nr_seeds.tsv"
for p in (O, OSD):
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
ann = {}
with open(D + "rep_annot.tsv") as fh:
    next(fh)
    for l in fh:
        f = l.rstrip("\n").split("\t"); ann[f[0]] = f[1:]
def units(r):
    a = ann.get(r)
    if not a: return set(), set(), set()
    pf = {x for x in a[0].split(",") if x != "-"}
    ko = {x for x in a[1].split(",") if x not in ("-", "")}
    md = {x for x in a[2].split(",") if x not in ("-", "")}
    return pf, ko, md
mem = {}; fam_members = collections.defaultdict(list)
with open(D + "context_members.tsv") as fh:
    h = next(fh).rstrip("\n").split("\t")
    for l in fh:
        d = dict(zip(h, l.rstrip("\n").split("\t"))); k = (d["set"], d["family"])
        mem[(k, d["member"])] = d; fam_members[k].append(d["member"])
nb = collections.defaultdict(list)
with open(D + "context_neighbours.tsv") as fh:
    next(fh)
    for l in fh:
        st, fam, m, off, ng, rs, dist, r95, r50 = l.rstrip("\n").split("\t")
        nb[((st, fam), m)].append((int(off), r95, r50, rs, int(dist)))
own50 = {}
need = {m for (_, m) in mem}
with open("/bigdata/stajichlab/lshad003/wf_protein_space/results/clusters_50_n88.tsv") as fh:
    for l in fh:
        r, m = l.rstrip("\n").split("\t")
        if m in need: own50[m] = r
print("members %d, own 50%% cluster found %d" % (len(need), len(own50)), flush=True)
def sim(a, b):
    ca, cb = collections.Counter(a), collections.Counter(b)
    return sum((ca & cb).values()) / max(len(a), len(b))
BINS = ((1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10 ** 9, ">=10"))
def nbin(n):
    for a, b, s in BINS:
        if a <= n <= b: return s
    return "0"
rows = []
osd = open(OSD, "w"); osd.write("set\tfamily\tseed_member\tn_contexts\n")
for k, ms in fam_members.items():
    ctx = []
    for m in ms:
        v = sorted(nb.get((k, m), []))
        if v: ctx.append((m, v))
    ctx.sort(key=lambda x: (-len(x[1]), x[0]))
    seeds = []
    for m, v in ctx:
        rl = [x[1] for x in v]
        for s in seeds:
            if sim(rl, s[2]) >= 0.8: s[3] += 1; break
        else: seeds.append([m, v, rl, 1])
    for s in seeds: osd.write("%s\t%s\t%s\t%d\n" % (k[0], k[1], s[0], s[3]))
    n = len(seeds); c50, cpf, cko, cmd = collections.Counter(), collections.Counter(), collections.Counter(), collections.Counter()
    annotated = 0
    for m, v, _, _ in seeds:
        o50 = own50.get(m); s50, spf, sko, smd = set(), set(), set(), set()
        for off, r95, r50, rs, dist in v:
            if r50 != "NA" and r50 != o50: s50.add(r50)
            pf, ko, md = units(r95); spf |= pf; sko |= ko; smd |= md
        c50.update(s50); cpf.update(spf); cko.update(sko); cmd.update(smd)
        annotated += bool(spf or sko)
    top = c50.most_common(1)[0] if c50 else ("-", 0)
    strand_c = dist_c = med_gap = top_ann = "NA"
    if top[1]:
        rel, offs, gaps = collections.Counter(), collections.Counter(), []
        for m, v, _, _ in seeds:
            hits = [x for x in v if x[2] == top[0]]
            if hits:
                x = min(hits, key=lambda x: (abs(x[0]), x[0])); rel[x[3]] += 1; offs[x[0]] += 1; gaps.append(x[4])
        strand_c = "%s:%.3f" % (rel.most_common(1)[0][0], rel.most_common(1)[0][1] / top[1])
        dist_c = "%+d:%.3f" % (offs.most_common(1)[0][0], offs.most_common(1)[0][1] / top[1])
        gaps.sort(); med_gap = str(gaps[len(gaps) // 2])
        r95s = collections.Counter(x[1] for m, v, _, _ in seeds for x in v if x[2] == top[0])
        a = ann.get(r95s.most_common(1)[0][0])
        top_ann = ";".join(y for y in (a[0], a[3], a[4]) if y != "-") if a else "-"
        top_ann = top_ann or "-"
    fmt = lambda c: ",".join("%s(%.2f)" % (x, y / n) for x, y in c.most_common(3)) if n else "-"
    rows.append([k[0], k[1], len(ms), len(ctx), len(ms) - len(ctx), n, nbin(n), "%.3f" % (annotated / n) if n else "NA",
                 top[0], "%.3f" % (top[1] / n) if n else "NA", top_ann, strand_c, dist_c, med_gap,
                 fmt(cpf) or "-", fmt(cko) or "-", fmt(cmd) or "-"])
osd.close()
ctl = collections.defaultdict(list)
for r in rows:
    if r[0] == "control" and r[9] != "NA": ctl[r[6]].append(float(r[9]))
def pct(x, b):
    v = ctl.get(b)
    if not v: return "NA"
    return "%.1f" % (100.0 * (sum(y < x for y in v) + 0.5 * sum(y == x for y in v)) / len(v))
with open(O, "w") as o:
    o.write("set\tfamily\tn_strict_members\tn_members_with_context\tn_members_no_neighbour\tn_nr_contexts\tnr_bin\t"
            "frac_nr_annotated\ttop_neighbour50\tconservation\ttop_neighbour_annotation\tstrand_conservation\t"
            "distance_conservation\tmedian_gap_bp\ttop_pfam\ttop_kegg_ko\ttop_kegg_module\tpercentile_vs_controls\tconserved_context\n")
    for r in rows:
        p = pct(float(r[9]), r[6]) if r[9] != "NA" else "NA"
        cons = int(r[5] >= 3 and r[9] != "NA" and float(r[9]) >= 0.5)
        o.write("\t".join(map(str, r + [p, cons])) + "\n")
for st in ("EU", "control"):
    R = [r for r in rows if r[0] == st]
    print("%s: units %d | with >= 1 NR context %d | >= 3 NR %d | conserved (>= 3 NR, >= 0.5) %d" % (st, len(R), sum(r[5] >= 1 for r in R),
          sum(r[5] >= 3 for r in R), sum(r[5] >= 3 and float(r[9]) >= 0.5 for r in R if r[9] != "NA")))
for b in [x[2] for x in BINS]:
    v = sorted(ctl.get(b, []))
    if v: print("control conservation, NR bin %s: n %d, median %.3f, 90th %.3f" % (b, len(v), v[len(v) // 2], v[int(0.9 * (len(v) - 1))]))
print("wrote", O)
