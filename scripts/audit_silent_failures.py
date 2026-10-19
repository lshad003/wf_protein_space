# Audit: past jobs that Slurm recorded as COMPLETED although their log shows a failure (2026-10-19; companion to
#   scripts/audit_exit_codes.py). Masked exit codes (a trailing echo or date) let a crashed job end COMPLETED.
# Input: argv[1] = sacct export (JobID|JobName|State|ExitCode|WorkDir|SubmitLine, -X -P, all jobs since 2026-07-01 plus
#   queued jobs; array tasks keyed as <job>_<task>). Logs: logs/* whose name carries a job id (.<id>.log, .<id>_<task>.log;
#   also .out/.err); logs without a job id are counted and skipped.
# Hits, per log line: printed non-zero exit ("exit 1", "exit code 2", "search exit 1 ..."), a Python traceback header,
#   or a word ending in "Error" (Error, KeyError, MemoryError ...). Only jobs whose State is COMPLETED are reported.
# Expected output (heuristic): file names written by the job's script or the python/R scripts it calls (strings that
#   end in a data extension, taken from lines that assign OUT*/O*, redirect with >, open(..., "w"), cp, write.table,
#   fwrite, saveRDS, to_csv); each resolved under results/ or the project root; exists yes/no; named in CHATINDEX.md
#   yes/no (basename match).
# Reviewed verdicts (2026-10-19, by reading each flagged log, the script, the rerun and CHATINDEX.md) are in REVIEW below,
#   keyed by script; reviewed_output is checked live for existence (and non-zero size) and for its name in CHATINDEX.md.
# Output: results/audit_silent_failures.tsv; refuses to overwrite.
import os, re, sys, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; L = W + "/logs/"; SC = W + "/scripts/"
OUT = W + "/results/audit_silent_failures.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
jobs = {}; rec = None
start = re.compile(r"^(\d+(?:_\d+)?)(?:_\[[^|]*\])?\|")
for l in open(sys.argv[1]):
    m = start.match(l)
    if m:
        f = l.rstrip("\n").split("|", 5)
        rec = {"state": f[2].split()[0] if f[2] else "", "exit": f[3], "sub": f[5] if len(f) > 5 else ""}
        jobs.setdefault(m.group(1), rec)
    elif rec is not None: rec["sub"] += "\n" + l.rstrip("\n")
def script_of(sub):
    t = sub.split()
    if t and t[0] == "sbatch":
        for x in t[1:]:
            if x.endswith(".sh"): return os.path.basename(x)
        if "--wrap" in sub: return "(--wrap)"
    return "NA"
chat = open(W + "/CHATINDEX.md").read()
DAL = "/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Ls_MAG_C/dallas/"
REVIEW = {
 "dallas74_quant167.sh": (DAL + "MAG167_quant/UHM102.10840/bin_abundance_table.tab", "MAG167_quant",
   "expected: metaWRAP python2 splitter fails on numpy and the script falls back to scripts/dallas75_split_bins_py3.py; "
   "first-run tables (29206748) were empty, rerun 29210315 wrote all 50 non-empty; documented in CHATINDEX"),
 "dallas57_drep1048.sh": (DAL + "drep1048_full50/dereplicated_genomes", "drep1048_full50",
   "benign: conda deactivate message at start; dRep finished (167 representatives)"),
 "step2b_subsample_diagnostic.sh": (W + "/results/step2b_hits.m8", "step2b_hits.m8",
   "real failure (MMseqs2 prefilter died); output never written and not cited; diagnostic superseded by step2e"),
 "step30b_eu15580_neighbours.sh": (W + "/results/step30b_family_summary.tsv", "step30b_family_summary.tsv",
   "real failure in the report after the files were written; rerun 29959130 reproduced every number; documented in CHATINDEX"),
 "step31c_def84_summary.sh": (W + "/results/step31_def84_summary.tsv", "step31_def84_summary.tsv",
   "real failure (ABORT, CA count != length) masked by trailing date; rerun 29977327 (wrapper now exits with rc) wrote the "
   "cited file; documented in CHATINDEX"),
 "step34j_b3a_merizo.sh": (W + "/results/step34j/domains", "step34j_b3a2_domains.py",
   "real failure (KeyError 'filename') masked by trailing date; domains written 2026-10-19 by scripts/step34j_b3a2_domains.py; "
   "documented in the Step 34j fix block"),
 "step34j_envs2.sh": ("NA", "step34j_envs2",
   "environment build: DefenseFinder version check failed; completed by step34j_envs3.sh (job 30050175); no data output"),
 "step8e_test.sh": ("NA", "step8e",
   "test of emapper --annotate_hits_table, printed FAIL (0 rows); origin of the CLAUDE.md emapper rule; no data output"),
}
def reviewed(sc):
    if sc not in REVIEW: return ("not reviewed", "NA", "NA", "NA")
    path, key, verdict = REVIEW[sc]
    if path == "NA": ex = "NA (no data output)"
    elif os.path.isdir(path): ex = "yes (%d entries)" % len(os.listdir(path))
    elif os.path.exists(path): ex = "yes (%d bytes)" % os.path.getsize(path)
    else: ex = "no"
    return (verdict, path, ex, "yes" if key in chat else "no")
EXT = r"(?:tsv|txt|faa|fasta|fna|fa|m8|gz|csv|pdf|png|svg|rds|RDS|paf|bam|tblout|domtblout|json|done|log)"
WR = re.compile(r"(\bOUT\w*\s*=|\bO\w*\s*=|>\s*\S|open\([^)]*['\"]w['\"]|\bcp\s|write\.table|fwrite|saveRDS|to_csv|writeLines)")
NM = re.compile(r"""['"/$+ ]((?:results/|catalog/)?[\w.\-]+(?:/[\w.\-]+)*\.%s)\b""" % EXT)
def outputs(script):
    files = []; seen = set()
    if script in ("NA", "(--wrap)") or not os.path.exists(SC + script): return files
    todo = [SC + script]
    for l in open(SC + script):
        for c in re.findall(r"scripts/([\w.\-]+\.(?:py|R))", l): todo.append(SC + c)
    for p in todo:
        if p in seen or not os.path.exists(p): continue
        seen.add(p)
        for l in open(p, errors="replace"):
            if l.lstrip().startswith("#") or not WR.search(l): continue
            for n in NM.findall(l):
                if n.endswith(".log") or "$" in n: continue
                if n not in files: files.append(n)
    return files
def resolve(n):
    for c in (W + "/" + n, W + "/results/" + n):
        if os.path.exists(c): return "yes"
    return "no"
PAT = [("exit_nonzero", re.compile(r"(?:^|[^A-Za-z./_\-$])exit(?: code| status)?[ :=]+([1-9]\d*)\b")),
       ("traceback", re.compile(r"Traceback \(most recent call last\)")),
       ("error_word", re.compile(r"\b[A-Za-z]*Error\b"))]
idre = re.compile(r"[._\-](\d{6,}(?:_\d+)?)\.(?:log|out|err)$")
rows = []; noid = 0; nolog = 0; states = collections.Counter()
for f in sorted(os.listdir(L)):
    m = idre.search(f)
    if not m: noid += 1; continue
    j = m.group(1); r = jobs.get(j)
    if r is None and "_" in j: r = jobs.get(j.split("_")[0])
    if r is None: nolog += 1; continue
    states[r["state"]] += 1
    if r["state"] != "COMPLETED": continue
    hits = collections.OrderedDict()
    with open(L + f, errors="replace") as fh:
        for no, line in enumerate(fh, 1):
            for nm, p in PAT:
                if p.search(line): hits.setdefault(nm, []).append((no, line.strip()[:200]))
    if not hits: continue
    sc = script_of(r["sub"]); outs = outputs(sc)
    ex = ["%s:%s" % (n, resolve(n)) for n in outs[:12]]
    ci = [os.path.basename(n) for n in outs if os.path.basename(n) in chat]
    for nm, h in hits.items():
        rows.append((j, f, sc, nm, len(h), "line %d: %s" % h[0], ";".join(ex) if ex else "NA",
                     "yes: " + ";".join(sorted(set(ci))[:8]) if ci else "no") + reviewed(sc))
with open(OUT, "w") as o:
    o.write("job\tlog\tscript\thit_type\tn_hit_lines\tfirst_hit\tauto_outputs_exist\tauto_output_named_in_CHATINDEX\t"
            "reviewed_verdict\treviewed_output\treviewed_output_exists\treviewed_named_in_CHATINDEX\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("logs without a job id %d; with a job id but no sacct record %d; states of matched logs %s" % (noid, nolog, dict(states)))
print("COMPLETED jobs with hits: %d; rows %d" % (len(set(r[0] for r in rows)), len(rows)))
print("by hit type:", dict(collections.Counter(r[3] for r in rows)))
print("wrote", OUT)
