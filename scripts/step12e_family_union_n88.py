# Step 12e n88: family-tier (50%) dark fraction under the eggNOG+Pfam union on
# the 88-metagenome catalog. Families and member metagenomes from
# results/clusters_50_n88.tsv; known sets restricted to the n88 representatives
# (results/known_union_ids_n88.txt, eggnog_hit_ids_n88.txt).
# Family known = any member protein is in the given known set (95% rep IDs).
# Support = >=3 members from >=2 metagenomes, matching step11.
# Cross-check: families, supported families and eggNOG-dark supported families
# must equal step11_n88 (results/family_dark_fraction_n88.tsv).
import sys, csv
W="/bigdata/stajichlab/lshad003/wf_protein_space"
union=set(open(W+"/results/known_union_ids_n88.txt").read().split())
egg=set(open(W+"/results/eggnog_hit_ids_n88.txt").read().split())
print("union ids:",len(union),"egg ids:",len(egg),flush=True)
MULTI="\x00M"
cnt={};first={};uk=set();ek=set();n=0
for line in open(W+"/results/clusters_50_n88.tsv"):
    rep,mem=line.rstrip("\n").split("\t")[:2]
    stem=mem.split("__",1)[0]
    cnt[rep]=cnt.get(rep,0)+1
    f=first.get(rep)
    if f is None: first[rep]=sys.intern(stem)
    elif f!=MULTI and f!=stem: first[rep]=MULTI
    if mem in union: uk.add(rep)
    if mem in egg: ek.add(rep)
    n+=1
    if n%10000000==0: print("lines:",n,flush=True)
fams=list(cnt)
sup=[r for r in fams if cnt[r]>=3 and first[r]==MULTI]
T=len(fams); S=len(sup)
eds=sum(1 for r in sup if r not in ek)
uds=sum(1 for r in sup if r not in uk)
eda=sum(1 for r in fams if r not in ek)
uda=sum(1 for r in fams if r not in uk)
out=open(W+"/results/family_union_dark_n88.tsv","w")
out.write("# n88. known: any member in eggNOG set (egg cols) or eggNOG+Pfam union (union cols)\n")
out.write("# support: >=3 members from >=2 metagenomes, 50%% tier\n")
out.write("scope\tfamilies\tegg_known\tegg_dark\tunion_known\tunion_dark\n")
for row in [("all",T,T-eda,eda,T-uda,uda),("supported",S,S-eds,eds,S-uds,uds)]:
    out.write("%s\t%d\t%d\t%d\t%d\t%d\n"%row); print("%s\t%d\t%d\t%d\t%d\t%d"%row)
out.close()
print("pct all: egg_dark %.2f%%, union_dark %.2f%%"%(100.0*eda/T,100.0*uda/T))
print("pct supported: egg_dark %.2f%%, union_dark %.2f%%"%(100.0*eds/S,100.0*uds/S))
ref={r[0]:r for r in csv.reader(open(W+"/results/family_dark_fraction_n88.tsv"),delimiter="\t")}
a=ref["all families"]; s=ref[">=3 members and >=2 metagenomes"]
ok=(T==int(a[1]) and S==int(s[1]) and eds==int(s[1])-int(s[2]))
print("CROSS-CHECK vs step11_n88 (%s, %s, %d):"%(a[1],s[1],int(s[1])-int(s[2])),"PASS" if ok else "MISMATCH, table written for inspection only")
sys.exit(0 if ok else 1)
