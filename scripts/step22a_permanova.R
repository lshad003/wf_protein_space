# Step 22a: community-level tests of treatment on the full gene matrix.
# 1. PERMANOVA with egg mass and month fitted first, permutations restricted
#    within animal so the longitudinal design is respected.
# 2. db-RDA constrained on treatment, conditioning out egg mass and month.
#    This finds the axis that best separates treatments, rather than picking
#    an unconstrained axis after the fact.
# 3. Dispersion check, since PERMANOVA can be driven by spread not location.
suppressMessages({library(vegan); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
set.seed(1)

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
meta$treatment <- factor(meta$treatment); meta$egg_mass <- factor(meta$egg_mass)
meta$month <- factor(meta$month); meta$animal <- factor(meta$animal)

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)
cnt <- cnt[, meta$stem, drop=FALSE]
cnt <- cnt[rowSums(cnt > 0) >= 35, , drop=FALSE]
cat("genes:", nrow(cnt), " samples:", ncol(cnt), "\n")

# CPM then Hellinger, standard for compositional gene abundance
x <- t(cpm(calcNormFactors(DGEList(counts=as.matrix(cnt)), method="TMM")))
x <- decostand(x, "hellinger")
d <- vegdist(x, method="bray")

cat("\n=== 1. PERMANOVA, egg mass and month fitted before treatment ===\n")
h <- how(nperm=999, blocks=meta$animal)
print(adonis2(d ~ egg_mass + month + treatment, data=meta,
              permutations=h, by="terms"))

cat("\n=== 2. dispersion by treatment (should be non-significant) ===\n")
print(anova(betadisper(d, meta$treatment)))

cat("\n=== 3. db-RDA constrained on treatment, egg mass and month removed ===\n")
m <- dbrda(d ~ treatment + Condition(egg_mass + month), data=meta)
print(m)
cat("\nvariance explained by treatment axis:",
    round(100*m$CCA$tot.chi/m$tot.chi, 2), "%\n")
cat("\npermutation test of the constrained axes:\n")
print(anova(m, permutations=h))
cat("\nper-axis test:\n")
print(anova(m, by="axis", permutations=h))

cat("\n=== 4. WF24, no repeated sampling ===\n")
m24 <- read.delim(paste0(W,"/metadata/samples_89_treatment.tsv"), stringsAsFactors=FALSE)
m24 <- m24[m24$cohort=="WF24", ]; rownames(m24) <- m24$stem
m24$treatment <- factor(m24$treatment)
c24 <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)[, m24$stem, drop=FALSE]
c24 <- c24[rowSums(c24 > 0) >= 30, , drop=FALSE]
cat("genes:", nrow(c24), " samples:", ncol(c24), "\n")
x24 <- decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(c24)), method="TMM"))), "hellinger")
d24 <- vegdist(x24, method="bray")
print(adonis2(d24 ~ treatment, data=m24, permutations=999))
cat("\nDONE step22a\n")
