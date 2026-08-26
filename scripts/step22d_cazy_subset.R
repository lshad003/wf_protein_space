# Step 22d: treatment tested on the carbohydrate-active subset only.
# No CAZy annotation exists for this catalogue, so the subset is defined by
# Pfam family name, which is a looser proxy and is stated as such.
# Rationale: a treatment effect confined to a small functional class could be
# diluted across 368,236 genes and missed by a whole-catalogue test.
suppressMessages({library(vegan); library(edgeR); library(limma)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
set.seed(1)

# Pfam name patterns for carbohydrate-active families
pat <- "^Glyco_hydro|^Glycos_transf|^CBM|^Chitin_bind|^Polysacc|^Cellulase|^Chitinase|^Glyco_tran|^Glyco_trans"
dom <- read.table(pipe(paste0("cat ", W, "/results/pfam_hs/*.tblout | grep -v '^#' | awk '$3 ~ /",
   "Glyco_hydro|Glycos_transf|CBM_|Chitin_bind|Polysacc|Cellulase|Chitinase/ {print $1}' | sort -u")),
   col.names="rep", stringsAsFactors=FALSE)
cz <- unique(dom$rep)
cat("representatives with a carbohydrate-active Pfam family:", length(cz), "\n")

meta <- read.delim(paste0(W,"/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
rownames(meta) <- meta$stem
for (v in c("treatment","egg_mass","month","animal")) meta[[v]] <- factor(meta[[v]])
meta$treatment <- relevel(meta$treatment, ref="control")

cnt <- read.delim(gzfile(paste0(W,"/results/count_matrix_primary.tsv.gz")),
                  row.names=1, check.names=FALSE)[, meta$stem, drop=FALSE]
cnt <- cnt[rownames(cnt) %in% cz, , drop=FALSE]
cnt <- cnt[rowSums(cnt > 0) >= 35, , drop=FALSE]
cat("carbohydrate genes present in >=35 of 44:", nrow(cnt), "\n")

cat("\n=== A. per-gene, limma blocked on animal ===\n")
d <- calcNormFactors(DGEList(counts=as.matrix(cnt)), method="TMM")
des <- model.matrix(~ treatment + egg_mass + month, data=meta)
v <- voom(d, des)
cf <- duplicateCorrelation(v, des, block=meta$animal)
v <- voom(d, des, block=meta$animal, correlation=cf$consensus)
cf <- duplicateCorrelation(v, des, block=meta$animal)
cat("within-animal correlation:", round(cf$consensus,4), "\n")
fit <- eBayes(lmFit(v, des, block=meta$animal, correlation=cf$consensus))
for (c1 in colnames(coef(fit))) {
  if (c1 == "(Intercept)") next
  t <- topTable(fit, coef=c1, number=Inf, sort.by="none")
  cat(sprintf("  %-24s padj<0.05: %5d\n", c1, sum(t$adj.P.Val < 0.05, na.rm=TRUE)))
  write.csv(t[order(t$adj.P.Val), ],
    paste0(W,"/results/limma_cazy_", gsub("[^A-Za-z0-9]","_",c1), ".csv"))
}

cat("\n=== B. community test on the carbohydrate subset, animal centroids ===\n")
x <- decostand(t(cpm(d)), "hellinger")
cen <- t(sapply(levels(meta$animal), function(a) colMeans(x[meta$animal==a, , drop=FALSE])))
am <- unique(meta[, c("animal","treatment","egg_mass")]); rownames(am) <- am$animal
am <- am[rownames(cen), ]
print(adonis2(vegdist(cen, method="bray") ~ egg_mass + treatment, data=am,
              permutations=999, by="terms"))

cat("\n=== C. total carbohydrate gene abundance per sample ===\n")
tot <- colSums(cpm(d)[, , drop=FALSE])
print(kruskal.test(tot ~ meta$treatment))
print(round(tapply(tot, meta$treatment, median)))
cat("\nDONE step22d\n")
