# Step 22b: WF22 community tests with the correct permutation scheme.
# Treatment and egg mass are constant within animal, so they must be tested by
# permuting whole animals. Month varies within animal and is tested by
# permuting timepoints inside each animal. Step 22a permuted freely within
# animal blocks, which never moves a between-animal label; its egg mass and
# treatment p-values are void.
suppressMessages({library(vegan); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
set.seed(1)

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
meta$treatment <- factor(meta$treatment); meta$egg_mass <- factor(meta$egg_mass)
meta$month <- factor(meta$month); meta$animal <- factor(meta$animal)

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)[, meta$stem, drop=FALSE]
cnt <- cnt[rowSums(cnt > 0) >= 35, , drop=FALSE]
x <- decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(cnt)), method="TMM"))), "hellinger")
d <- vegdist(x, method="bray")
cat("genes:", nrow(cnt), " samples:", ncol(cnt), " animals:", nlevels(meta$animal), "\n")

cat("\n=== A. between-animal terms: permute whole animals, balanced subset ===\n")
# vegan requires a balanced design to permute strata, and one animal has 2
# samples rather than 3, so that animal is dropped for this test only.
tb <- table(meta$animal)
kb <- names(tb)[tb == max(tb)]
mb <- meta[meta$animal %in% kb, ]
mb$animal <- droplevels(mb$animal)
db <- as.dist(as.matrix(d)[mb$stem, mb$stem])
cat("balanced subset:", nrow(mb), "samples from", nlevels(mb$animal), "animals\n")
hB <- how(nperm=999, plots=Plots(strata=mb$animal, type="free"),
          within=Within(type="none"))
print(adonis2(db ~ egg_mass + treatment, data=mb, permutations=hB, by="terms"))

cat("\n=== B. within-animal term: permute timepoints inside animals ===\n")
hW <- how(nperm=999, plots=Plots(strata=meta$animal, type="none"),
          within=Within(type="free"))
print(adonis2(d ~ month, data=meta, permutations=hW, by="terms"))

cat("\n=== C. animal-level test, one centroid per animal ===\n")
# collapses the repeated measures entirely: 15 independent units
am <- unique(meta[, c("animal","treatment","egg_mass")])
rownames(am) <- am$animal
cen <- t(sapply(levels(meta$animal), function(a) colMeans(x[meta$animal==a, , drop=FALSE])))
dc <- vegdist(cen, method="bray")
am <- am[rownames(cen), ]
cat("units:", nrow(am), "\n")
print(adonis2(dc ~ egg_mass + treatment, data=am, permutations=999, by="terms"))

cat("\n=== D. dispersion by treatment and egg mass (animal centroids) ===\n")
cat("treatment:\n"); print(anova(betadisper(dc, am$treatment)))
cat("egg mass:\n");  print(anova(betadisper(dc, am$egg_mass)))
cat("\nDONE step22b\n")
