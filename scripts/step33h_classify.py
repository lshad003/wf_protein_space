# Step 33h classify: where the EU reps on fly contigs fall in the D. melanogaster genome (dm6, GCF_000001215.4).
# Sets: (1) every EU supported rep (results/fourway_classes_n88.tsv class4 EU, as defined incl. UHM586.41010) whose own
#   contig is fly95 in results/step33g_contig_fly.tsv; (2) the family reps of the fly-majority families among the 15,580
#   (results/step33g_family_fly.tsv majority_fly95 = 1); a family rep whose own contig is not fly95 is reported as such.
# Set-2 reps are direct-30% family reps (step26e1), not necessarily 95% reps: their sequence comes from
#   results/prostt5_main_eu15580.faa when absent from catalog/db/LsPS_AA_95_rep.fasta.
# Prodigal coordinates (start, end, strand, partial) from the per-sample Prodigal headers ("k141_X_N # start # end #
#   strand # ID=..;partial=.."): catalog/by_assembly_new/<stem>.aa.fa.gz if present, else
#   /bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly/<stem>.aa.fa.gz. Each source
#   protein must equal the catalog sequence (catalog/db/LsPS_AA_95_rep.fasta; terminal '*' ignored); aborts otherwise.
# Projection: the rep's contig alignment in results/step33g_paf/<stem>.fly.paf (minimap2 -c, cg:Z) with the largest
#   overlap with the ORF; >= 50% of the ORF must lie in its aligned query range, else "not projected". The ORF ends
#   (clipped to the alignment) are mapped through the CIGAR (positions in query insertions map to the current target
#   base); projected strand = ORF strand x alignment strand.
# Classes (overlap threshold MIN = 30 bp, or the projected length if shorter; first match wins):
#   CDS same strand, in frame | CDS same strand, other frame | CDS opposite strand | UTR or non-coding exon |
#   intron (inside a gene, pseudogene or transcript span, no exon overlap >= MIN) | intergenic.
#   Frame: at the first base p shared with a same-strand CDS, the ORF codon position (from the ORF 5' end, clipped
#   offset included; indels between the anchor and p ignored) is compared with the CDS codon position from its GFF
#   phase; in frame if equal for any overlapping same-strand CDS.
# Repeat flags (separate): fraction of the projected span in soft-masked (lowercase) dm6 sequence; overlap >= MIN with
#   a GFF mobile_genetic_element.
# Outputs (refuse to overwrite): results/step33h_rep_class.tsv, results/step33h_queries.faa (projected-set proteins
#   without '*', for DIAMOND)
import os, re, gzip, bisect, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"; FLY = W + "/catalog/fly/"
OUT = [R + "step33h_rep_class.tsv", R + "step33h_queries.faa"]
for p in OUT:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
PRO_NEW = W + "/catalog/by_assembly_new/%s.aa.fa.gz"
PRO_OLD = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly/%s.aa.fa.gz"
MIN = 30
def contig(g): s, rest = g.split("__", 1); return s, rest[:rest.rfind("_")]
fly = set()
for i, l in enumerate(open(R + "step33g_contig_fly.tsv")):
    if i:
        f = l.rstrip("\n").split("\t")
        if f[6] == "1": fly.add((f[0], f[1]))
eu = []
with open(R + "fourway_classes_n88.tsv") as fh:
    next(fh)
    for l in fh:
        r, c = l.split("\t", 2)[:2]
        if c == "EU" and contig(r) in fly: eu.append(r)
famrep = [l.split("\t")[0] for i, l in enumerate(open(R + "step33g_family_fly.tsv")) if i and l.rstrip("\n").split("\t")[4] == "1"]
fam_onfly = [r for r in famrep if contig(r) in fly]
print("set 1, EU supported reps on fly95 contigs: %d" % len(eu))
print("set 2, reps of fly-majority families: %d; own contig fly95: %d; not: %d" % (len(famrep), len(fam_onfly), len(famrep) - len(fam_onfly)), flush=True)
S1, S2 = set(eu), set(fam_onfly); want = S1 | S2; bystem = collections.defaultdict(set)
for r in want: bystem[r.split("__", 1)[0]].add(r)
cat = {}; cur = None
for l in open(W + "/catalog/db/LsPS_AA_95_rep.fasta"):
    if l[0] == ">":
        h = l[1:].split()[0]; cur = h if h in want else None
        if cur: cat[h] = []
    elif cur: cat[cur].append(l.strip())
cat = {k: "".join(v).rstrip("*") for k, v in cat.items()}
# family reps (step26e1 direct 30% clustering of member proteins) need not be 95% reps: take the rest from
# results/prostt5_main_eu15580.faa (the 15,580 family rep sequences, terminal * removed here)
n95 = len(cat); cur = None
for l in open(R + "prostt5_main_eu15580.faa"):
    if l[0] == ">":
        h = l[1:].split()[0]; cur = h if (h in want and h not in cat) else None
        if cur: cat[h] = ""
    elif cur: cat[cur] += l.strip()
cat = {k: v.rstrip("*") for k, v in cat.items()}
print("catalog sequences: %d from LsPS_AA_95_rep.fasta, %d from prostt5_main_eu15580.faa" % (n95, len(cat) - n95), flush=True)
if set(cat) != want: raise SystemExit("ABORT: %d reps without a catalog sequence" % len(want - set(cat)))
orf, srcs = {}, collections.Counter()
for s, ids in bystem.items():
    p = PRO_NEW % s if os.path.exists(PRO_NEW % s) else PRO_OLD % s
    if not os.path.exists(p): raise SystemExit("ABORT: no Prodigal file for " + s)
    srcs["by_assembly_new" if p == PRO_NEW % s else "Fecal/Proteins/by_assembly"] += 1
    cur, seq = None, []
    for l in list(gzip.open(p, "rt")) + [">"]:
        if l[0] == ">":
            if cur:
                q = "".join(seq).rstrip("*")
                if q != cat[cur]: raise SystemExit("ABORT: %s Prodigal sequence differs from the catalog" % cur)
            cur, seq = None, []
            if len(l) > 1:
                f = l[1:].split(" # "); r = s + "__" + f[0]
                if r in ids:
                    cur = r; pa = re.search(r"partial=(\d\d)", f[4]).group(1)
                    orf[r] = (int(f[1]) - 1, int(f[2]) - 1, 1 if f[3] == "1" else -1, pa)
        elif cur: seq.append(l.strip())
if set(orf) != want: raise SystemExit("ABORT: %d reps without Prodigal coordinates" % len(want - set(orf)))
print("Prodigal coordinates for %d reps (sequences equal the catalog) | sources: %s" % (len(orf), dict(srcs)), flush=True)
need = collections.defaultdict(set)
for r in want: need[r.split("__", 1)[0]].add(contig(r)[1])
aln = collections.defaultdict(list)
for s, cs in need.items():
    for l in open(R + "step33g_paf/%s.fly.paf" % s):
        f = l.rstrip("\n").split("\t")
        if f[0] in cs:
            cg = re.search(r"\tcg:Z:(\S+)", l)
            if not cg: raise SystemExit("ABORT: PAF row without cg:Z")
            aln[(s, f[0])].append((int(f[2]), int(f[3]) - 1, 1 if f[4] == "+" else -1, f[5], int(f[7]), int(f[8]) - 1, re.findall(r"(\d+)([MIDNSHP=X])", cg.group(1))))
def mapq(a, x):
    qs, qe, st, tn, ts, te, ops = a
    q = qs if st == 1 else qe; t = ts
    for n, op in ops:
        n = int(n)
        if op in "M=X":
            if st == 1 and q <= x < q + n: return t + (x - q)
            if st == -1 and q - n < x <= q: return t + (q - x)
            q += n * st; t += n
        elif op == "I":
            if (st == 1 and q <= x < q + n) or (st == -1 and q - n < x <= q): return t
            q += n * st
        elif op in "DN": t += n
    raise ValueError("position %d outside alignment" % x)
seqs, cur = {}, None
for l in gzip.open(FLY + "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz", "rt"):
    if l[0] == ">": cur = l[1:].split()[0]; seqs[cur] = []
    else: seqs[cur].append(l.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
feat = collections.defaultdict(lambda: collections.defaultdict(list))
SPAN = {"gene", "pseudogene", "mRNA", "lnc_RNA", "antisense_RNA", "miRNA", "tRNA", "snoRNA", "primary_transcript", "rRNA", "ncRNA", "snRNA", "pre_miRNA", "transcript", "RNase_MRP_RNA", "RNase_P_RNA", "SRP_RNA"}
for l in gzip.open(FLY + "GCF_000001215.4_genomic.gff.gz", "rt"):
    if l[0] == "#": continue
    f = l.split("\t")
    if len(f) < 9: continue
    ty = f[2]; k = "CDS" if ty == "CDS" else "exon" if ty == "exon" else "mge" if ty == "mobile_genetic_element" else "span" if ty in SPAN else None
    if k: feat[k][f[0]].append((int(f[3]) - 1, int(f[4]) - 1, 1 if f[6] == "+" else -1, int(f[7]) if f[7] in "012" else 0))
IDX = {}
for k, d in feat.items():
    for c, v in d.items():
        v.sort(); IDX[(k, c)] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v) + 1)
def hits(k, c, lo, hi):
    if (k, c) not in IDX: return []
    st, v, ml = IDX[(k, c)]; i = bisect.bisect_left(st, lo - ml); out = []
    while i < len(v) and v[i][0] <= hi:
        a = v[i]; ov = min(a[1], hi) - max(a[0], lo) + 1
        if ov > 0: out.append((a, ov))
        i += 1
    return out
CL = ["CDS same strand, in frame", "CDS same strand, other frame", "CDS opposite strand", "UTR or non-coding exon", "intron", "intergenic", "not projected"]
res = {}
for r in sorted(want):
    s, c = contig(r); g0, g1, gs, pa = orf[r]; L = g1 - g0 + 1
    best, bo = None, 0
    for a in aln.get((s, c), []):
        ov = min(a[1], g1) - max(a[0], g0) + 1
        if ov > bo: best, bo = a, ov
    if best is None or bo < 0.5 * L: res[r] = ("not projected", "", -1, -1, 0, 0.0, 0, pa, L); continue
    c0, c1 = max(g0, best[0]), min(g1, best[1])
    ta, tb = mapq(best, c0), mapq(best, c1); lo, hi = min(ta, tb), max(ta, tb); ps = gs * best[2]; chrom = best[3]
    q5 = g0 if gs == 1 else g1; anc = c0 if gs == 1 else c1; d = abs(anc - q5); tan = mapq(best, anc)
    m = min(MIN, hi - lo + 1)
    cds = hits("CDS", chrom, lo, hi); same = [(a, ov) for a, ov in cds if a[2] == ps and ov >= m]; opp = [(a, ov) for a, ov in cds if a[2] != ps and ov >= m]
    if same:
        inf = False
        for a, ov in same:
            p = max(lo, a[0]); rp = (d + p - tan) % 3 if ps == 1 else (d + tan - p) % 3
            cp = (p - a[0] - a[3]) % 3 if a[2] == 1 else (a[1] - p - a[3]) % 3
            inf |= rp == cp
        k = CL[0] if inf else CL[1]
    elif opp: k = CL[2]
    elif any(ov >= m for _, ov in hits("exon", chrom, lo, hi)): k = CL[3]
    elif any(ov >= m for _, ov in hits("span", chrom, lo, hi)): k = CL[4]
    else: k = CL[5]
    sq = seqs[chrom][lo:hi + 1]; mf = sum(ch.islower() for ch in sq) / len(sq)
    te = int(any(ov >= m for _, ov in hits("mge", chrom, lo, hi)))
    res[r] = (k, chrom, lo, hi, ps, mf, te, pa, L)
with open(OUT[0], "w") as o:
    o.write("rep\tin_set1_fly_EU\tin_set2_flymajority_famrep\tclass\tchrom\tstart0\tend0\tproj_strand\tsoftmasked_frac\tmge_overlap\tprodigal_partial\torf_nt\n")
    for r in sorted(want):
        k, ch, lo, hi, ps, mf, te, pa, L = res[r]
        o.write("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%.3f\t%d\t%s\t%d\n" % (r, r in S1, r in S2, k, ch or "NA", lo, hi, {1: "+", -1: "-", 0: "NA"}[ps], mf, te, pa, L))
with open(OUT[1], "w") as o:
    for r in sorted(want): o.write(">%s\n%s\n" % (r, cat[r]))
print("classes (all projected-set reps):", dict(collections.Counter(v[0] for v in res.values())))
for p in OUT: print("wrote", p)
