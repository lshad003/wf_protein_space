# Step 28h inputs: the Step 28 ProstT5/Foldseek query sets with the stop character removed.
# Steps 28d/28f/28g wrote catalog sequences unchanged, so every query ended in '*' (decoys began with '*'),
# and ProstT5 embedded it. Here the terminal '*' is removed as in step29a; aborts on an internal stop.
#  results/prostt5_main_eu15580_nostar.faa          from results/prostt5_main_eu15580.faa (15,580)
#  results/prostt5_pilot_eu740_nostar.faa           from results/prostt5_pilot_eu740.faa (740)
#  results/prostt5_calib_known_lenmatched_nostar.faa from results/prostt5_calib_known_lenmatched.faa (740)
#  results/prostt5_decoy_eu740rev_nostar.faa        the stripped pilot sequences reversed, headers "rev_" (as step28f_sets)
# Same IDs and order as the originals. Refuses to overwrite.
import os
W = "/bigdata/stajichlab/lshad003/wf_protein_space"; R = W + "/results/"
SETS = (("prostt5_main_eu15580", 15580), ("prostt5_pilot_eu740", 740), ("prostt5_calib_known_lenmatched", 740))
DEC = R + "prostt5_decoy_eu740rev_nostar.faa"
for p in [R + s + "_nostar.faa" for s, _ in SETS] + [DEC]:
    if os.path.exists(p): raise SystemExit("REFUSING: " + p + " exists")
def read(faa):
    seq, cur = {}, None
    for l in open(faa):
        if l[0] == ">": cur = l[1:].split()[0]; seq[cur] = []
        else: seq[cur].append(l.strip())
    return {k: "".join(v) for k, v in seq.items()}
def write(p, items):
    with open(p, "w") as fh:
        for k, s in items: fh.write(">%s\n%s\n" % (k, s))
    print("wrote %s %d" % (p, len(items)))
strip = {}
for s, n in SETS:
    q = read(R + s + ".faa")
    if len(q) != n: raise SystemExit("ABORT: %s has %d sequences, expected %d" % (s, len(q), n))
    st = {k: v.rstrip("*") for k, v in q.items()}
    bad = [k for k, v in st.items() if "*" in v]
    if bad: raise SystemExit("ABORT: internal stop in %s: %s" % (s, bad[:5]))
    print("%s: %d sequences, %d had a terminal '*'" % (s, n, sum(len(q[k]) != len(st[k]) for k in q)))
    write(R + s + "_nostar.faa", list(st.items())); strip[s] = st
dec = [("rev_" + k, v[::-1]) for k, v in strip["prostt5_pilot_eu740"].items()]
old = read(R + "prostt5_decoy_eu740rev.faa")
if [k for k, _ in dec] != list(old): raise SystemExit("ABORT: decoy IDs or order differ from prostt5_decoy_eu740rev.faa")
print("decoys identical to the original decoys after removing a leading '*': %d of %d"
      % (sum(old[k].lstrip("*") == v for k, v in dec), len(dec)))
write(DEC, dec)
