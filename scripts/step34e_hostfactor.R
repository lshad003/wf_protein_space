# Step 34e: host-factor tests on one gene set (usage: Rscript step34e_hostfactor.R <label> <counts.tsv.gz>).
# Designs reused unchanged (prevalence filters, normalisation, models, permutation schemes):
#   step22b_permanova_fix_n88.R  WF22 Bray-Curtis on Hellinger TMM-CPM, genes in >= 35 of 44; egg_mass + treatment by
#                                whole-animal permutation on the balanced subset (42 samples, 14 animals); month
#                                permuted within animal; animal centroids (15); betadisper on centroids
#   step22c_jaccard_n88.R        WF22 Jaccard presence (>= 5 reads, variable genes), same terms; richness Kruskal
#                                (month within animal added here for presence; not in step22c)
#   step22e_permanova_wf24_n34_n88.R  WF24 34 animals (metadata/wf24_excluded.tsv), abundance genes in >= 30, presence
#   step17g_limma_wf22_n88.R     WF22 limma-voom TMM, ~ treatment + egg_mass + month, genes in >= 35, naive and
#                                blocked on animal (duplicateCorrelation, two rounds), EM3 STP1710.7 vs control
#   step17h_limma_wf2324_n88.R   WF23 limma-voom ~ treatment (ref Control), genes in >= 8 of 9
#   step17i_limma_wf24_n34_n88.R WF24 limma-voom ~ treatment (ref 4), genes in >= 30 of 34
# Added here: WF23 community tests (UHM520.7734 vs Control; abundance genes in >= 8 of 9, presence >= 5 reads; free
#   permutation, 126 distinct labelings), and WF24 with all 35 animals as sensitivity (UHM590.41012 included).
# Normalisation (TMM) is computed within the gene set analysed. Treatment and egg mass are constant within WF22
#   animals: under within-animal permutation they are untestable (recorded as such), not null.
# Outputs: results/step34e_<label>_permanova.tsv, results/step34e_<label>_limma.tsv,
#   results/step34e_limma/<label>_<model>_<coef>.csv.gz
suppressMessages({library(vegan); library(edgeR); library(limma)})
a <- commandArgs(TRUE); LAB <- a[1]; CF <- a[2]
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"; R <- paste0(W, "/results/")
O1 <- paste0(R, "step34e_", LAB, "_permanova.tsv"); O2 <- paste0(R, "step34e_", LAB, "_limma.tsv")
if (file.exists(O1) || file.exists(O2)) stop("REFUSING: outputs exist for ", LAB)
dir.create(paste0(R, "step34e_limma"), showWarnings=FALSE)
cat("gene set:", LAB, "| counts:", CF, "\n")
cnt0 <- read.delim(gzfile(CF), row.names=1, check.names=FALSE)
cat("genes in set:", nrow(cnt0), " columns:", ncol(cnt0), "\n")
P <- list(); L <- list()
addp <- function(cohort, metric, test, res, note="") {
  r <- as.data.frame(res)
  for (t in rownames(r)) if (!(t %in% c("Residual", "Total")))
    P[[length(P) + 1]] <<- data.frame(set=LAB, cohort=cohort, metric=metric, test=test, term=t, df=r[t, "Df"],
      R2=round(r[t, "R2"], 5), F=round(r[t, "F"], 4), p=r[t, "Pr(>F)"], note=note)
}
addu <- function(cohort, metric, test, term, note)
  P[[length(P) + 1]] <<- data.frame(set=LAB, cohort=cohort, metric=metric, test=test, term=term, df=NA, R2=NA, F=NA, p=NA, note=note)
addd <- function(cohort, metric, term, an)
  P[[length(P) + 1]] <<- data.frame(set=LAB, cohort=cohort, metric=metric, test="betadisper_centroids", term=term,
    df=an$Df[1], R2=NA, F=round(an$`F value`[1], 4), p=an$`Pr(>F)`[1], note="dispersion")
rep <- function(fit, cohort, model, ngenes) {
  for (cf in colnames(coef(fit))) {
    if (cf == "(Intercept)") next
    t <- topTable(fit, coef=cf, number=Inf, sort.by="none")
    n5 <- sum(t$adj.P.Val < 0.05, na.rm=TRUE); n5f <- sum(t$adj.P.Val < 0.05 & abs(t$logFC) > 1, na.rm=TRUE)
    cat(sprintf("  %-10s %-10s %-22s genes %7d  padj<0.05: %6d   +|LFC|>1: %6d\n", cohort, model, cf, ngenes, n5, n5f))
    L[[length(L) + 1]] <<- data.frame(set=LAB, cohort=cohort, model=model, coef=cf, genes_tested=ngenes, padj_lt_0.05=n5, padj_lt_0.05_lfc_gt1=n5f)
    write.csv(t[order(t$adj.P.Val), ], gzfile(paste0(R, "step34e_limma/", LAB, "_", cohort, "_", model, "_", gsub("[^A-Za-z0-9]", "_", cf), ".csv.gz")))
  }
}
hell <- function(cn) decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(cn)), method="TMM"))), "hellinger")
pres <- function(cn) { pa <- t((cn >= 5) * 1); pa[, colSums(pa) > 0 & colSums(pa) < nrow(pa), drop=FALSE] }
# ---------------- WF22
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE); rownames(meta) <- meta$stem
for (v in c("treatment", "egg_mass", "month", "animal")) meta[[v]] <- factor(meta[[v]])
meta$treatment <- relevel(meta$treatment, ref="control")
c22 <- cnt0[, meta$stem, drop=FALSE]
a22 <- c22[rowSums(c22 > 0) >= 35, , drop=FALSE]
cat("\nWF22 abundance genes in >= 35 of 44:", nrow(a22), "\n")
tb <- table(meta$animal); kb <- names(tb)[tb == max(tb)]
mb <- meta[meta$animal %in% kb, ]; mb$animal <- droplevels(mb$animal)
hB <- how(nperm=999, plots=Plots(strata=mb$animal, type="free"), within=Within(type="none"))
hW <- how(nperm=999, plots=Plots(strata=meta$animal, type="none"), within=Within(type="free"))
am <- unique(meta[, c("animal", "treatment", "egg_mass")]); rownames(am) <- am$animal
for (metric in c("abundance_bray_hellinger", "presence_jaccard_ge5")) {
  if (metric == "abundance_bray_hellinger") { x <- hell(a22); d <- vegdist(x, method="bray") }
  else { x <- pres(c22); d <- vegdist(x, method="jaccard", binary=TRUE); cat("WF22 presence variable genes:", ncol(x), "\n") }
  db <- as.dist(as.matrix(d)[mb$stem, mb$stem])
  set.seed(1); r <- adonis2(db ~ egg_mass + treatment, data=mb, permutations=hB, by="terms"); print(r)
  addp("WF22", metric, "whole_animal_perm_balanced_42", r)
  set.seed(1); r <- adonis2(d ~ month, data=meta, permutations=hW, by="terms"); print(r)
  addp("WF22", metric, "within_animal_perm_44", r)
  for (t in c("egg_mass", "treatment")) addu("WF22", metric, "within_animal_perm_44", t, "untestable: constant within animal")
  cen <- t(sapply(levels(meta$animal), function(a) colMeans(x[meta$animal == a, , drop=FALSE])))
  dc <- vegdist(cen, method="bray"); amc <- am[rownames(cen), ]
  set.seed(1); r <- adonis2(dc ~ egg_mass + treatment, data=amc, permutations=999, by="terms"); print(r)
  addp("WF22", metric, "animal_centroids_15", r)
  addd("WF22", metric, "treatment", anova(betadisper(dc, amc$treatment)))
  addd("WF22", metric, "egg_mass", anova(betadisper(dc, amc$egg_mass)))
  if (metric == "presence_jaccard_ge5") {
    kw <- kruskal.test(rowSums(x) ~ meta$treatment)
    P[[length(P) + 1]] <- data.frame(set=LAB, cohort="WF22", metric=metric, test="richness_kruskal_44", term="treatment",
      df=kw$parameter, R2=NA, F=round(kw$statistic, 4), p=kw$p.value, note="genes detected per sample; F column = chi-squared")
  }
}
cat("\nWF22 limma-voom\n")
d <- calcNormFactors(DGEList(counts=as.matrix(a22)), method="TMM")
des <- model.matrix(~ treatment + egg_mass + month, data=meta)
v <- voom(d, des); rep(eBayes(lmFit(v, des)), "WF22", "B_naive", nrow(a22))
v <- voom(d, des); cf <- duplicateCorrelation(v, des, block=meta$animal)
v <- voom(d, des, block=meta$animal, correlation=cf$consensus); cf <- duplicateCorrelation(v, des, block=meta$animal)
cat("refined within-animal correlation:", round(cf$consensus, 4), "\n")
L[[length(L) + 1]] <- data.frame(set=LAB, cohort="WF22", model="A_blocked", coef="within_animal_correlation", genes_tested=nrow(a22),
  padj_lt_0.05=round(cf$consensus, 4), padj_lt_0.05_lfc_gt1=NA)
rep(eBayes(lmFit(v, des, block=meta$animal, correlation=cf$consensus)), "WF22", "A_blocked", nrow(a22))
m3 <- meta[meta$egg_mass == "3" & meta$treatment %in% c("control", "STP1710.7"), ]
m3$treatment <- droplevels(m3$treatment); m3$month <- droplevels(m3$month)
d3 <- calcNormFactors(DGEList(counts=as.matrix(a22[, m3$stem, drop=FALSE])), method="TMM")
des3 <- model.matrix(~ treatment + month, data=m3)
v3 <- voom(d3, des3); c3 <- duplicateCorrelation(v3, des3, block=m3$animal)
v3 <- voom(d3, des3, block=m3$animal, correlation=c3$consensus); c3 <- duplicateCorrelation(v3, des3, block=m3$animal)
rep(eBayes(lmFit(v3, des3, block=m3$animal, correlation=c3$consensus)), "WF22", "C_EM3", nrow(a22))
# ---------------- WF23
mt <- read.delim(paste0(W, "/metadata/samples_88_treatment.tsv"), stringsAsFactors=FALSE); rownames(mt) <- mt$stem
m23 <- mt[mt$cohort == "WF23", ]; m23$treatment <- relevel(factor(m23$treatment), ref="Control")
c23 <- cnt0[, m23$stem, drop=FALSE]; a23 <- c23[rowSums(c23 > 0) >= 8, , drop=FALSE]
cat("\nWF23 genes in >= 8 of 9:", nrow(a23), "\n")
set.seed(1); r <- adonis2(vegdist(hell(a23), method="bray") ~ treatment, data=m23, permutations=999); print(r)
addp("WF23", "abundance_bray_hellinger", "free_perm_9", r, "126 distinct labelings")
x <- pres(c23); set.seed(1); r <- adonis2(vegdist(x, method="jaccard", binary=TRUE) ~ treatment, data=m23, permutations=999); print(r)
addp("WF23", "presence_jaccard_ge5", "free_perm_9", r, "126 distinct labelings")
d <- calcNormFactors(DGEList(counts=as.matrix(a23)), method="TMM"); des <- model.matrix(~ treatment, data=m23)
rep(eBayes(lmFit(voom(d, des), des)), "WF23", "limma", nrow(a23))
# ---------------- WF24 (34 primary, 35 sensitivity)
ex <- read.delim(paste0(W, "/metadata/wf24_excluded.tsv"), stringsAsFactors=FALSE)
for (k in c("34", "35")) {
  m24 <- mt[mt$cohort == "WF24", ]
  if (k == "34") m24 <- m24[!(m24$stem %in% ex$stem), ]
  stopifnot(nrow(m24) == as.integer(k))
  m24$treatment <- relevel(factor(m24$treatment), ref="4"); co <- paste0("WF24_n", k)
  c24 <- cnt0[, m24$stem, drop=FALSE]; mp <- if (k == "34") 30 else 31
  a24 <- c24[rowSums(c24 > 0) >= mp, , drop=FALSE]
  cat("\n", co, " genes in >= ", mp, ": ", nrow(a24), "\n", sep="")
  set.seed(1); r <- adonis2(vegdist(hell(a24), method="bray") ~ treatment, data=m24, permutations=999); print(r)
  addp(co, "abundance_bray_hellinger", "free_perm", r, paste0("genes in >= ", mp))
  x <- pres(c24); set.seed(1); r <- adonis2(vegdist(x, method="jaccard", binary=TRUE) ~ treatment, data=m24, permutations=999); print(r)
  addp(co, "presence_jaccard_ge5", "free_perm", r)
  d <- calcNormFactors(DGEList(counts=as.matrix(a24)), method="TMM"); des <- model.matrix(~ treatment, data=m24)
  rep(eBayes(lmFit(voom(d, des), des)), co, "limma", nrow(a24))
}
write.table(do.call(rbind, P), O1, sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, L), O2, sep="\t", quote=FALSE, row.names=FALSE)
cat("wrote", O1, "and", O2, "\nDONE step34e", LAB, "\n")
