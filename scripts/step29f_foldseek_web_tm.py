# Step 29f: TM-scores from the Foldseek web-server JSON downloads (3Di/AA, 2026-10-14). The JSON holds Foldseek's
# alignment and C-alpha coordinates for query and target but no TM-score. TM-score is computed on Foldseek's fixed
# alignment (seeded Kabsch superposition with iterative refinement, standard d0), normalised by query (qTM) and
# target (tTM) length. Top 25 hits per database by E-value.
import glob, json, os, numpy as np
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; D = W + "/results/foldseek_web4"
id2seq, cur = {}, None
for l in open(W + "/results/esmfold_input.faa"):
    l = l.strip()
    if not l: continue
    if l[0] == ">": cur = l[1:].split()[0]; id2seq[cur] = []
    else: id2seq[cur].append(l)
seq2id = {"".join(v).rstrip("*"): k for k, v in id2seq.items()}
def coords(s): return np.array([float(x) for x in str(s).split(",")]).reshape(-1, 3)
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T; return R, qc - pc @ R.T
def d0(L): return max(0.5, 1.24 * (L - 15) ** (1 / 3) - 1.8) if L > 21 else 0.5
def tm(P, Q, L):
    n = len(P); D0 = d0(L); best = 0.0
    if n < 3: return 0.0
    for fl in sorted({max(4, n // k) for k in (1, 2, 4, 8)}, reverse=True):
        for s in range(0, n - fl + 1, max(1, fl // 2)):
            idx = np.arange(s, s + fl)
            for _ in range(20):
                R, t = kabsch(P[idx], Q[idx]); d = np.linalg.norm(P @ R.T + t - Q, axis=1)
                best = max(best, float(np.sum(1 / (1 + (d / D0) ** 2)) / L)); cut = D0; new = np.where(d < cut)[0]
                while len(new) < 3 and cut < 50: cut += 0.5; new = np.where(d < cut)[0]
                if np.array_equal(new, idx): break
                idx = new
    return best
def pairs(h, qC):
    tC = coords(h["tCa"]); qi, ti = h["qStartPos"] - 1, h["dbStartPos"] - 1; P, Q = [], []
    for a, b in zip(h["qAln"], h["dbAln"]):
        if a != "-" and b != "-": P.append(qC[qi]); Q.append(tC[ti])
        qi += a != "-"; ti += b != "-"
    return np.array(P), np.array(Q)
out = open(W + "/results/step29f_foldseek_web_tm.tsv", "w")
out.write("protein\tjson\tdb\ttarget\tdescription\ttaxon\tprob\tevalue\tseqid\tq_start\tq_end\tq_len\tt_start\tt_end\tt_len\tn_aligned\tqTM\ttTM\n")
for f in sorted(glob.glob(D + "/Foldseek_*.json")):
    d = json.load(open(f))[0]; seq = d["query"]["sequence"].rstrip("*"); pid = seq2id.get(seq)
    if pid is None: raise SystemExit("ABORT: no esmfold_input sequence matches " + os.path.basename(f))
    qC = coords(d["query"]["qCa"])
    if len(qC) != len(seq): raise SystemExit("ABORT: CA count %d != length %d in %s" % (len(qC), len(seq), f))
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    if tm(qC, qC @ rot.T + 5, len(seq)) < 0.999: raise SystemExit("ABORT: self-test failed")
    print("\n== %s <- %s | %d aa | mode %s" % (pid, os.path.basename(f), len(seq), d.get("mode")))
    for r in d["results"]:
        al = r["alignments"]; hits = []
        for v in (al.values() if isinstance(al, dict) else al): hits += v if isinstance(v, list) else [v]
        top = sorted(hits, key=lambda h: float(h["eval"]))[:25]; shown = 0
        for h in top:
            P, Q = pairs(h, qC); qtm, ttm = tm(P, Q, len(seq)), tm(P, Q, int(h["dbLen"]))
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n" % (pid, os.path.basename(f), r["db"], h["target"],
                str(h.get("description", "")).replace("\t", " "), h.get("taxName", ""), h["prob"], h["eval"], h["seqId"], h["qStartPos"], h["qEndPos"],
                h["qLen"], h["dbStartPos"], h["dbEndPos"], h["dbLen"], len(P), qtm, ttm))
            if shown < 2: print("  %-14s %-28s prob %.2f E %.1e qTM %.2f tTM %.2f | %s" % (r["db"], str(h["target"])[:28], float(h["prob"]), float(h["eval"]), qtm, ttm, str(h.get("description", ""))[:45])); shown += 1
        if not top: print("  %-14s no hits" % r["db"])
out.close()
