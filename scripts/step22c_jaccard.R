# Step 22c: presence/absence version of the community test. Asks whether
# treatment changes which genes are present, rather than how abundant they are.
# Same permutation logic as step22b: whole animals for between-animal terms.
suppressMessages({library(vegan); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
set.seed(1)

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
for (v in c("treatment","egg_mass","month","animal")) meta[[v]] <- factor(meta[[v]])

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)[, meta$stem, drop=FALSE]
# detection threshold stated: a gene counts as present at >=5 reads
pa <- t((cnt >= 5) * 1)
pa <- pa[, colSums(pa) > 0 & colSums(pa) < nrow(pa), drop=FALSE]
cat("genes variable in presence:", ncol(pa), " samples:", nrow(pa), "\n")
cat("mean genes present per sample:", round(mean(rowSums(pa))), "\n")
d <- vegdist(pa, method="jaccard", binary=TRUE)

cat("\n=== A. between-animal terms, balanced subset, whole-animal permutation ===\n")
tb <- table(meta$animal); kb <- names(tb)[tb == max(tb)]
mb <- meta[meta$animal %in% kb, ]; mb$animal <- droplevels(mb$animal)
db <- as.dist(as.matrix(d)[mb$stem, mb$stem])
hB <- how(nperm=999, plots=Plots(strata=mb$animal, type="free"),
          within=Within(type="none"))
print(adonis2(db ~ egg_mass + treatment, data=mb, permutations=hB, by="terms"))

cat("\n=== B. animal centroids, 15 independent units ===\n")
cen <- t(sapply(levels(meta$animal), function(a) colMeans(pa[meta$animal==a, , drop=FALSE])))
am <- unique(meta[, c("animal","treatment","egg_mass")]); rownames(am) <- am$animal
am <- am[rownames(cen), ]
print(adonis2(vegdist(cen, method="bray") ~ egg_mass + treatment, data=am,
              permutations=999, by="terms"))

cat("\n=== C. richness: number of genes detected per sample ===\n")
r <- rowSums(pa)
print(kruskal.test(r ~ meta$treatment))
print(tapply(r, meta$treatment, function(z) round(c(median=median(z), min=min(z), max=max(z)))))

cat("\n=== D. WF24 presence/absence, 36 independent animals ===\n")
m24 <- read.delim(paste0(W,"/metadata/samples_89_treatment.tsv"), stringsAsFactors=FALSE)
m24 <- m24[m24$cohort=="WF24", ]; rownames(m24) <- m24$stem
m24$treatment <- factor(m24$treatment)
c24 <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)[, m24$stem, drop=FALSE]
p24 <- t((c24 >= 5) * 1)
p24 <- p24[, colSums(p24) > 0 & colSums(p24) < nrow(p24), drop=FALSE]
cat("genes variable in presence:", ncol(p24), "\n")
print(adonis2(vegdist(p24, method="jaccard", binary=TRUE) ~ treatment,
              data=m24, permutations=999))
cat("\nDONE step22c\n")
