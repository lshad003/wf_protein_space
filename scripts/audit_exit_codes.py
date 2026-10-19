# Audit: do sbatch scripts exit with their main command's status? (2026-10-19, after Step 34j Merizo job 30050461 crashed
#   but ended with date, so Slurm recorded COMPLETED and afterok dependencies released.)
# Scripts audited: scripts/*.sh with #SBATCH headers, plus any scripts/*.sh named in an sbatch submit line.
# Submit lines from sacct (SubmitLine; all jobs of this user since 2026-07-01), exported to the file given as argv[1]
#   with: sacct -X -n -P -u $USER -S 2026-07-01 -E now -o JobID,JobName,State,ExitCode,WorkDir,SubmitLine, followed by
#   the same fields for every job in squeue (sacct -j <ids>), because jobs pending on a dependency have no eligible time
#   and are missing from the -S/-E window.
# Parsing (heuristic, bash not executed): full-line comments and heredoc bodies skipped (the line opening a heredoc is
#   kept); lines split on ";" outside quotes. A command followed by "|| { ...; exit N; }" (N != 0) is kept as a guarded
#   command (its failure propagates); a one-line test of $? or $rc that exits non-zero, or "[ $? = 0 ] && ...", is
#   kept as a status check; other guard lines (one-line if...fi, or || / && with exit) are ignored;
#   structural tokens (then, else, do, done, fi, {, }) ignored. Trivial commands: echo printf date hostname rm cp mv ls wc
#   du cat mkdir sleep module cd export touch head tail grep cut sort uniq test [ nproc free df pwd sync, assignments
#   without a command substitution, and literal exit inside a block.
# main command = last non-trivial statement; final = last statement executed at the end of the file.
# masked = no when final is the main command, the main command is guarded, a status check directly follows the main
#   command (or its rc=$?; a later check of a variable set by V=$? right after the main command also counts), or final is "exit $V" with "V=$?" (or "exit $?") directly after the main command;
#   otherwise yes. Scripts whose main command sits in a loop report the last iteration only (noted).
# dependency chain = yes when a job of the script was submitted with --dependency/-d (downstream) or a job of the script
#   was named in another job's dependency (upstream).
# Output: results/audit_exit_codes.tsv; refuses to overwrite.
import os, re, sys, collections
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; SC = W + "/scripts/"
OUT = W + "/results/audit_exit_codes.tsv"
if os.path.exists(OUT): raise SystemExit("REFUSING: " + OUT + " exists")
# --- submit lines
jobs = {}; rec = None
start = re.compile(r"^(\d+)(_[\[\]\d,\-%]+)?\|")
for l in open(sys.argv[1]):
    m = start.match(l)
    if m:
        f = l.rstrip("\n").split("|", 5)
        rec = {"id": m.group(1), "state": f[2], "workdir": f[4], "sub": f[5] if len(f) > 5 else ""}
        jobs.setdefault(m.group(1), rec)
    elif rec is not None:
        rec["sub"] += "\n" + l.rstrip("\n")
def script_of(sub):
    toks = sub.split()
    if not toks or toks[0] != "sbatch": return None
    for t in toks[1:]:
        if t.endswith(".sh"): return os.path.basename(t)
    return None
job_script = {}; down = collections.defaultdict(int); upids = set()
dep_re = re.compile(r"(?:--dependency[= ]|(?<!\S)-d\s*)(\S+)")
for jid, r in jobs.items():
    s = script_of(r["sub"])
    if s: job_script[jid] = s
    for d in dep_re.findall(r["sub"]):
        ids = re.findall(r"(\d{5,})", d)
        if ids and s: down[s] += 1
        upids.update(ids)
up = collections.Counter(job_script[j] for j in upids if j in job_script)
nsub = collections.Counter(job_script.values())
# --- script parsing
TRIV = set("echo printf date hostname rm cp mv ls wc du cat mkdir sleep module cd export touch head tail grep cut sort "
           "uniq test [ [[ nproc free df pwd sync true ml".split())
STRUCT = {"then", "else", "do", "done", "fi", "{", "}", "esac", ";;"}
def split_semis(line):
    out, cur, q, i = [], "", None, 0
    while i < len(line):
        c = line[i]
        if q:
            cur += c
            if c == q: q = None
        elif c in "'\"": q = c; cur += c
        elif c == ";" and line[i:i + 2] != ";;": out.append(cur); cur = ""
        else: cur += c
        i += 1
    out.append(cur)
    return [x.strip() for x in out if x.strip()]
def parse(path):
    raw = open(path).read().split("\n"); lines = []; i = 0; buf = ""; start_no = None
    while i < len(raw):
        l = raw[i]; i += 1
        if buf == "": start_no = i
        if l.rstrip().endswith("\\"): buf += l.rstrip()[:-1] + " "; continue
        l = buf + l; buf = ""
        s = l.strip()
        if not s or s.startswith("#"): continue
        lines.append((start_no, s))
        m = re.search(r"<<-?\s*['\"]?(\w+)['\"]?", s)
        if m:
            tag = m.group(1)
            while i < len(raw) and raw[i].strip() != tag: i += 1
            i += 1
    st = []; depth = 0
    for no, s in lines:
        if re.match(r"^(\[\[?|if \[\[?|test)\s.*\$(\?|\{?\w*rc\w*\}?)(?!\w)", s) and (re.search(r"\bexit\s+[1-9]", s) or re.search(r"\$\?\s*(=|-eq)\s*0\s*\]+\s*&&", s)):
            st.append((no, s, "status", depth)); continue
        g = re.match(r"^(.*?)\s*\|\|\s*\{?.*\bexit\s+[1-9]", s)
        if g and g.group(1).split() and g.group(1).split()[0] not in TRIV and not g.group(1).startswith("if "):
            st.append((no, g.group(1).strip(), "guarded", depth)); continue
        if (s.startswith("if ") and re.search(r"\bfi\b", s)) or (re.search(r"\|\||&&", s) and re.search(r"\bexit\b", s)):
            continue
        for x in split_semis(s):
            head = x.split()[0]
            if head in ("if", "for", "while", "until", "case"):
                depth += 1
                if head == "if": continue
                if head in ("for", "while", "until"): st.append((no, x, "loop", depth)); continue
            if head in ("fi", "done", "esac"): depth -= 1; continue
            if head in STRUCT or x in STRUCT: continue
            if x.startswith("then ") or x.startswith("do ") or x.startswith("else "):
                x = x.split(None, 1)[1]; head = x.split()[0]
            st.append((no, x, "", depth))
    return st
def trivial(x, depth):
    head = x.split()[0]
    if re.match(r"^\w+=\$\?$", x): return True
    if re.match(r"^exit\b", x): return True
    x = re.sub(r"^(\w+=[^\s'\"]*\s+)+(?=\S)", "", x); head = x.split()[0]
    if re.match(r"^\w+=", x) and "$(" not in x and "`" not in x: return True
    if re.match(r"^\w+=\$\(\s*(ls|wc|grep|cat|date|echo|basename|dirname|awk|sed|head|tail|cut|sort|nproc)\b", x): return True
    return head in TRIV
rows = []
cands = sorted(f for f in os.listdir(SC) if f.endswith(".sh"))
for f in cands:
    txt = open(SC + f).read()
    if "#SBATCH" not in txt and f not in nsub: continue
    st = parse(SC + f)
    work = [s for s in st if s[2] != "loop"]
    main_i = max([i for i, s in enumerate(work) if s[2] != "status" and not trivial(s[1], s[3])], default=None)
    final = work[-1] if work else None
    note = []
    if main_i is None:
        masked = "unknown"; main = "NA"; note.append("no non-trivial command found")
    else:
        mn = work[main_i]; main = "line %d: %s" % (mn[0], mn[1][:200])
        if mn[3] > 0 and any(s[2] == "loop" for s in st): note.append("main command inside a block or loop")
        after = work[main_i + 1:main_i + 3]
        chk = (after and after[0][2] == "status") or (len(after) == 2 and re.match(r"^\w+=\$\?$", after[0][1]) and after[1][2] == "status")
        if after and re.match(r"^\w+=\$\?$", after[0][1]):
            v = after[0][1].split("=")[0]
            chk = chk or any(x[2] == "status" and re.search(r"\$\{?%s\}?(?!\w)" % v, x[1]) for x in work[main_i + 2:])
        if final is mn or mn[2] == "guarded" or chk: masked = "no"
        else:
            m = re.match(r"^exit\s+\$\{?(\w+|\?)\}?$", final[1])
            nxt = work[main_i + 1] if main_i + 1 < len(work) else None
            if m and nxt is not None and (nxt is final and m.group(1) == "?" or re.match(r"^%s=\$\?$" % re.escape(m.group(1)), nxt[1])):
                masked = "no"
            else: masked = "yes"
    dep = "yes" if (down.get(f) or up.get(f)) else "no"
    role = ",".join(x for x, c in (("downstream", down.get(f)), ("upstream", up.get(f))) if c) or "-"
    rows.append((f, main, masked, dep, role, nsub.get(f, 0), "line %d: %s" % (final[0], final[1][:120]) if final else "NA", "; ".join(note)))
rows.sort(key=lambda r: (r[3] != "yes", r[2] != "yes", r[0]))
with open(OUT, "w") as o:
    o.write("script\tmain_command_line\tmasked\tused_in_dependency_chain\tchain_role\tjobs_submitted_since_20260701\tfinal_statement\tnote\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
c = collections.Counter((r[3], r[2]) for r in rows)
print("scripts audited %d; (chain, masked) counts: %s" % (len(rows), dict(c)))
print("wrote", OUT)
