# Sensitivity for Section 2.6: WF22 PERMANOVA on months 2 and 3 only (30 metagenomes, 15 froglets x 2), where all egg
# masses span the same two sequencing runs. As Step 34e: Hellinger on TMM-CPM with Bray-Curtis (genes > 0 in >= 24 of
# 30 samples), Jaccard on presence >= 5 reads (variable genes); egg mass + treatment with whole-animal permutations,
# month (2 vs 3) within animal, per-animal centroids; run tested at the level where it varies. 999 permutations,
# set.seed(1). Output: results/step34l_m23_permanova.tsv
suppressMessages({library(vegan); library(edgeR); library(permute)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"; OUT <- paste0(W, "/results/step34l_m23_permanova.tsv")
if (file.exists(OUT)) stop("REFUSING: output exists")
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE)
bt <- read.delim(paste0(W, "/results/batch_table_88.tsv"), stringsAsFactors=FALSE)
meta$run <- bt$run[match(meta$stem, bt$stem)]
m <- meta[meta$month %in% c(2, 3), ]; rownames(m) <- m$stem
for (v in c("egg_mass", "month", "treatment", "run", "animal")) m[[v]] <- factor(m[[v]])
cat("samples:", nrow(m), "animals:", nlevels(m$animal), "\n"); print(table(m$egg_mass, m$run))
rpa <- tapply(as.character(m$run), m$animal, function(x) length(unique(x))); cat("animals with > 1 run:", sum(rpa > 1), "\n")
hell <- function(cn) decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(cn)), method="TMM"))), "hellinger")
res <- NULL
add <- function(set, metric, test, r) {
  tt <- as.data.frame(r)
  for (i in seq_len(nrow(tt))) if (!rownames(tt)[i] %in% c("Residual", "Total"))
    res <<- rbind(res, data.frame(set=set, metric=metric, test=test, term=rownames(tt)[i], R2=tt$R2[i], p=tt[["Pr(>F)"]][i]))
}
h1 <- how(nperm=999, plots=Plots(strata=m$animal, type="free"), within=Within(type="none"))
h2 <- how(nperm=999, plots=Plots(strata=m$animal, type="none"), within=Within(type="free"))
for (lab in c("strict_final", "strict_final_EU")) {
  cnt <- read.delim(gzfile(sprintf("%s/results/step34e_counts_%s.tsv.gz", W, lab)), row.names=1, check.names=FALSE)[, m$stem]
  a <- cnt[rowSums(cnt > 0) >= 24, ]; cat(lab, "genes in >= 24 of 30:", nrow(a), "\n")
  pr <- (cnt >= 5) * 1; pr <- pr[rowSums(pr) > 0 & rowSums(pr) < ncol(pr), ]
  for (metric in c("abundance_bray_hellinger", "presence_jaccard_ge5")) {
    if (metric == "abundance_bray_hellinger") { x <- as.matrix(hell(a)); d <- vegdist(x, "bray") } else { x <- t(as.matrix(pr)); d <- vegdist(x, "jaccard", binary=TRUE) }
    set.seed(1); add(lab, metric, "whole_animal", adonis2(d ~ egg_mass + treatment, data=m, permutations=h1, by="terms"))
    set.seed(1); add(lab, metric, "within_animal", adonis2(d ~ month, data=m, permutations=h2))
    if (any(rpa > 1)) { set.seed(1); add(lab, metric, "within_animal_run", adonis2(d ~ run, data=m, permutations=h2))
    } else { set.seed(1); add(lab, metric, "whole_animal_run_first", adonis2(d ~ run + egg_mass + treatment, data=m, permutations=h1, by="terms")) }
    cen <- apply(x, 2, function(v) tapply(v, m$animal, mean)); am <- m[!duplicated(m$animal), ]; rownames(am) <- am$animal; am <- am[rownames(cen), ]
    set.seed(1); add(lab, metric, "centroids", adonis2(vegdist(cen, "bray") ~ egg_mass + treatment, data=am, permutations=999, by="terms"))
  }
}
write.table(res, OUT, sep="\t", quote=FALSE, row.names=FALSE); print(res)
