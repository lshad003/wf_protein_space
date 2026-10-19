# Step 34j: keyword themes for the proposed concordance between a structural hit and genomic context (imported by
# step34j_a4_context_merge.py and step34j_c1_matrix.py). A theme is assigned when a regex matches the hit text (PDB
# entry COMPOUND, CATH superfamily name, AFDB UniProt protein name) or a neighbour annotation (Pfam ID and
# description, eggNOG preferred name and description, KEGG KO). Generic enzyme words are deliberately not themes.
# Proposed call, for manual review:
#   supports  = a hit theme is carried by >= 30% of non-redundant contexts (>= 2 NR contexts), or the matching mobile/
#               defence evidence (geNomad virus or provirus for phage; geNomad plasmid for mobile element; DefenseFinder
#               system within +-5 genes for defence) covers >= 30% of member contigs;
#   conflicts = the hit has a theme, the context is conserved (>= 3 NR, top neighbour >= 0.5), the top neighbour has a
#               theme, no hit theme is shared with it, and no hit theme reaches 10% of NR contexts;
#   neutral   = everything else (hit without a theme, unannotated or unconserved context).
import re
THEMES = [
    ("CRISPR-Cas", r"crispr|\bcas\d|\bcas[a-z]?\d|\bcse\d|\bcsy\d|\bcsm\d|\bcmr\d|cascade"),
    ("phage", r"phage|capsid|portal protein|\btail\b|tail fib|terminase|baseplate|virion|holin|endolysin|major head|prohead|tape measure"),
    ("mobile element", r"transposase|integrase|recombinase|resolvase|insertion (element|sequence)|\bIS\d|conjugat|relaxase|\bmob[A-Z]\b|plasmid|excisionase"),
    ("toxin-antitoxin", r"toxin|antitoxin|addiction|\brelE\b|\bparE\b|\bhic[AB]\b|\bmazF\b|\bvap[BC]\b|PIN domain|\bhipA\b|\byoeB\b"),
    ("defence/R-M", r"restriction|DNA methyl|methylase|modification methyl|abortive|abi[A-Z]|defen[cs]e|retron|gabija|thoeris|BREX|DISARM|zorya|septu|wadjet|lamassu|argonaute"),
    ("secretion/pilus/adhesion", r"secretion|secretin|\bpil[A-Z]\b|pilus|pilin|fimbri|type IV|type VI|T6SS|sortase|adhesin|\bHcp\b|VgrG|autotransporter"),
    ("flagellum/chemotaxis", r"flagell|motility|chemotaxis|\bmot[AB]\b|\bche[A-Z]\b"),
    ("transport/membrane", r"transporter|permease|\bABC\b|efflux|channel|porin|symporter|antiporter|\bTonB\b|outer membrane|MFS"),
    ("transcription regulation", r"transcription|regulator|\bHTH\b|helix-turn-helix|sigma factor|sigma-\d|repressor|activator|response regulator|DNA-binding"),
    ("DNA replication/repair", r"helicase|DNA polymerase|primase|recombination|repair|nuclease|topoisomerase|gyrase|DNA ligase|replication"),
    ("cell envelope/division", r"peptidoglycan|murein|cell wall|lipopolysaccharide|\bLPS\b|penicillin-binding|transpeptidase|cell division|\bfts[A-Z]\b|septum|S-layer|teichoic"),
    ("translation/ribosome", r"ribosom|\btRNA\b|tRNA |translation|elongation factor"),
    ("carbohydrate active", r"glycoside hydrolase|glycosyl hydrolase|carbohydrate|polysaccharide|amylase|xylan|cellul|chitin|glycosyltransferase|pectin|mannan|galactan"),
    ("sporulation", r"sporulation|spore|\bspo[0-9A-Z]"),
]
CRX = [(n, re.compile(p, re.I)) for n, p in THEMES]
def themes(text):
    if not text or text in ("-", "NA"): return set()
    return {n for n, r in CRX if r.search(text)}
MGE_THEME = {"virus": "phage", "provirus": "phage", "plasmid": "mobile element", "defence": "defence/R-M"}
def concordance(hit_text, nr_theme_frac, n_nr, conserved, top_nb_text, mge_frac):
    """hit_text: annotation of the structural hit; nr_theme_frac: theme -> fraction of NR contexts;
    mge_frac: {'virus','provirus','plasmid','defence'} -> fraction of member contigs; returns (call, evidence)."""
    ht = themes(hit_text)
    if not ht: return "neutral", "hit has no theme"
    sup = [(t, f) for t, f in nr_theme_frac.items() if t in ht and f >= 0.3 and n_nr >= 2]
    msup = [(k, f) for k, f in mge_frac.items() if MGE_THEME[k] in ht and f >= 0.3]
    if sup or msup:
        ev = ["%s in %.0f%% of %d NR contexts" % (t, 100 * f, n_nr) for t, f in sorted(sup, key=lambda x: -x[1])]
        ev += ["%s on %.0f%% of member contigs" % (k, 100 * f) for k, f in msup]
        return "supports", "hit theme %s; %s" % ("/".join(sorted(ht)), "; ".join(ev))
    tt = themes(top_nb_text)
    if conserved and tt and not (tt & ht) and all(nr_theme_frac.get(t, 0) < 0.1 for t in ht):
        return "conflicts", "hit theme %s; conserved top neighbour theme %s (%s)" % ("/".join(sorted(ht)), "/".join(sorted(tt)), top_nb_text[:80])
    best = max(((t, nr_theme_frac.get(t, 0)) for t in ht), key=lambda x: x[1])
    return "neutral", "hit theme %s; best context support %s %.0f%% of %d NR contexts" % ("/".join(sorted(ht)), best[0], 100 * best[1], n_nr)
