# PCoA coordinates for Figure 5: WF22 Bray-Curtis on Hellinger TMM-CPM, genes in >= 35 of 44, exactly as
# step34e_hostfactor.R, for strict_final and strict_final_EU. % variance is over the positive eigenvalues.
# Output: results/fig_pcoa_wf22.tsv
# Optional arguments (Step 34s): <counts dir> <output file>; reads step34e_counts_<set>.tsv.gz from <counts dir> and
#   refuses to overwrite <output file>. Default (no arguments) unchanged.
suppressMessages({library(vegan); library(edgeR)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
CD <- paste0(W, "/results"); OUT <- paste0(W, "/results/fig_pcoa_wf22.tsv"); aa <- commandArgs(TRUE)
if (length(aa) > 1) { CD <- sub("/+$", "", aa[1]); OUT <- aa[2]; if (file.exists(OUT)) stop("REFUSING: output exists") }
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE); rownames(meta) <- meta$stem
hell <- function(cn) decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(cn)), method="TMM"))), "hellinger")
out <- NULL
for (lab in c("strict_final", "strict_final_EU")) {
  cnt <- read.delim(gzfile(sprintf("%s/step34e_counts_%s.tsv.gz", CD, lab)), row.names=1, check.names=FALSE)
  c22 <- cnt[, meta$stem, drop=FALSE]
  if (ncol(c22) != 44) stop("expected 44 WF22 samples")
  a22 <- c22[rowSums(c22 > 0) >= 35, , drop=FALSE]
  cat(lab, ": WF22 abundance genes in >= 35 of 44:", nrow(a22), "\n")
  pc <- cmdscale(vegdist(hell(a22), method="bray"), k=2, eig=TRUE)
  pct <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
  cat(lab, ": PCo1", round(pct[1], 1), "% PCo2", round(pct[2], 1), "%\n")
  df <- data.frame(set=lab, stem=rownames(pc$points), PC1=pc$points[, 1], PC2=pc$points[, 2], pct1=pct[1], pct2=pct[2])
  out <- rbind(out, cbind(df, meta[df$stem, c("animal", "treatment", "egg_mass", "month")]))
}
write.table(out, OUT, sep="\t", quote=FALSE, row.names=FALSE)
cat("wrote", OUT, "\n")
