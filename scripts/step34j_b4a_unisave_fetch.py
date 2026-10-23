# Step 34j B4a: close the UniProt gap of b4 (job 30171185: 31,660 entries returned for 88,365 requested AFDB accessions;
#   missing list results/step34j/uniprot_gap/missing_accessions.tsv, scripts/step34j_b4_uniprot_gap.py).
# A. Secondary accessions: every entry in the saved b4 responses (results/step34j/uniprot/batch_*.tsv) is re-queried
#    (UniProt REST accessions endpoint, 500 per call, format=json without a fields filter: the TSV has no secondary
#    accession column, and 'sec_acc' is refused with HTTP 400), so a requested secondary accession that UniProt answered
#    with its primary can be mapped back. Raw responses gzipped: uniprot_unisave/sec_acc/batch_*.json.gz.
# B. UniSave: for each missing accession, GET /unisave/<acc>?format=json (version list), then the flat file of the highest
#    entryVersion (GET /unisave/<acc>?format=txt&versions=<n>); if the version list is empty or 404, GET
#    /unisave/<acc>/status. Raw bodies appended to uniprot_unisave/unisave_raw.jsonl (one line per accession; resumable:
#    accessions already recorded with a terminal status are skipped). Optional $1: a further accession list (header line
#    first), fetched in addition to the gap list.
# Rate: 3 worker threads sharing one limiter, at most 4 requests/s overall. Retries: 6, backoff 5-160 s, honouring
#   Retry-After, on network errors, 429 and 5xx; 400/404 are terminal and recorded.
# Exit 1 if any accession ends in status "failed" or any secondary batch fails.
import os, sys, json, glob, gzip, time, threading, urllib.request, urllib.error, concurrent.futures as CF
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; D = W + "/results/step34j/"; O = D + "uniprot_unisave/"
os.makedirs(O + "sec_acc", exist_ok=True)
UA = "wf_protein_space step34j b4a (academic; urllib)"
lock = threading.Lock(); last = [0.0]; MIN_GAP = 0.25
def get(url):
    for attempt in range(7):
        with lock:
            w = last[0] + MIN_GAP - time.time()
            if w > 0: time.sleep(w)
            last[0] = time.time()
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60) as r:
                return r.status, r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", "replace")
            if e.code in (400, 404): return e.code, body
            ra = e.headers.get("Retry-After")
            wait = int(ra) if ra and ra.isdigit() else 5 * 2 ** attempt
        except Exception as e:
            body = repr(e); wait = 5 * 2 ** attempt
        if attempt == 6: return -1, body
        time.sleep(wait)
# A
ret = []
for f in sorted(glob.glob(D + "uniprot/batch_*.tsv")):
    ret += [l.split("\t")[0] for l in list(open(f))[1:]]
ret = sorted(set(ret)); print("A: saved entries %d" % len(ret), flush=True)
badA = 0
for bi in range(0, len(ret), 500):
    out = O + "sec_acc/batch_%05d.json.gz" % (bi // 500)
    if os.path.exists(out) and os.path.getsize(out) > 0: continue
    st, body = get("https://rest.uniprot.org/uniprotkb/accessions?accessions=%s&format=json&size=500" % ",".join(ret[bi:bi + 500]))
    if st != 200: badA += 1; print("A: batch %d status %d" % (bi // 500, st), flush=True); continue
    with gzip.open(out + ".tmp", "wt") as g: g.write(body)
    os.rename(out + ".tmp", out)
print("A: done, failed batches %d" % badA, flush=True)
# B
miss = [l.strip() for l in open(D + "uniprot_gap/missing_accessions.tsv")][1:]
if len(sys.argv) > 1:   # 2026-10-23: optional extra list (results/step34j/uniprot_gap/b4_tie_missing.tsv, step34j_b4b_tie_accessions.py)
    miss += [x for x in (l.strip() for l in list(open(sys.argv[1]))[1:]) if x and x not in set(miss)]
RAW = O + "unisave_raw.jsonl"; done = set()
if os.path.exists(RAW):
    for l in open(RAW):
        try: r = json.loads(l)
        except ValueError: continue
        if r["status"] != "failed": done.add(r["acc"])
todo = [a for a in miss if a not in done]
print("B: missing %d, already done %d, to fetch %d" % (len(miss), len(done), len(todo)), flush=True)
wl = threading.Lock(); fh = open(RAW, "a"); cnt = {"n": 0}
def one(a):
    r = {"acc": a, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    st, body = get("https://rest.uniprot.org/unisave/%s?format=json" % a)
    r["json_status"], r["json_raw"] = st, body
    vers = []
    if st == 200:
        try: vers = [x["entryVersion"] for x in json.loads(body).get("results", [])]
        except (ValueError, KeyError, TypeError): pass
    if st == -1: r["status"] = "failed"
    elif vers:
        v = max(vers); r["version"] = v
        st2, b2 = get("https://rest.uniprot.org/unisave/%s?format=txt&versions=%d" % (a, v))
        r["txt_status"], r["txt_raw"] = st2, b2
        r["status"] = "ok" if st2 == 200 else "failed" if st2 == -1 else "txt_error"
    else:
        st3, b3 = get("https://rest.uniprot.org/unisave/%s/status" % a)
        r["status_status"], r["status_raw"] = st3, b3
        r["status"] = "failed" if st3 == -1 else "no_versions"
    with wl:
        fh.write(json.dumps(r) + "\n"); fh.flush(); cnt["n"] += 1
        if cnt["n"] % 1000 == 0: print("B: %d / %d %s" % (cnt["n"], len(todo), time.strftime("%H:%M:%S")), flush=True)
with CF.ThreadPoolExecutor(3) as ex: list(ex.map(one, todo))
fh.close()
final = {}
for l in open(RAW):
    try: r = json.loads(l)
    except ValueError: continue   # a line cut by a cancelled run
    final[r["acc"]] = r["status"]
import collections
c = collections.Counter(final.get(a, "absent") for a in miss)
print("B: final status over %d missing: %s" % (len(miss), dict(c)), flush=True)
sys.exit(1 if badA or c["failed"] or c["absent"] else 0)
