# Step 34s sensitivity (2026-10-23): does the WF22 month association hold within the two sequencing runs that month 1
# shares with months 2 and 3 (UCD_20240722, UCD_20240726)? Design copied from step34e_hostfactor.R (month within animal,
# 44 samples) and step34l_m23_permanova.R (months 2-3): c_id95 counts (results/step34s/c_id95/step34e_counts_<set>.tsv.gz),
# Hellinger on TMM-CPM computed within the gene set and sample set, Bray-Curtis, genes > 0 in round(0.8 n) samples (= 35 of
# 44 and 24 of 30 as in the sources), 999 permutations within animal (Plots none, Within free), set.seed(1), adonis2.
# A2 enters run before month (adonis2 by="terms", same within-animal permutations). Test B (4 animals x 3 months) has
# 3!^4 = 1,296 within-animal labelings; permute enumerates all of them (below minperm 5,040), so p is exact.
# Run labels from results/batch_table_89.tsv. Abundance only.
# Outputs: results/step34s/c_id95/month_shared_runs.tsv, results/step34s/c_id95/month1_excess_by_animal.tsv
suppressMessages({library(vegan); library(edgeR); library(permute)})
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"; CD <- paste0(W, "/results/step34s/c_id95")
O1 <- paste0(CD, "/month_shared_runs.tsv"); O2 <- paste0(CD, "/month1_excess_by_animal.tsv")
if (file.exists(O1) || file.exists(O2)) stop("REFUSING: outputs exist")
SHARED <- c("UCD_20240722", "UCD_20240726")
meta <- read.delim(paste0(W, "/metadata/wf22_design.tsv"), stringsAsFactors=FALSE); rownames(meta) <- meta$stem
bt <- read.delim(paste0(W, "/results/batch_table_89.tsv"), stringsAsFactors=FALSE)
meta$run <- bt$run[match(meta$stem, bt$stem)]; stopifnot(nrow(meta) == 44, !any(is.na(meta$run)))
for (v in c("egg_mass", "month", "animal", "treatment")) meta[[v]] <- factor(meta[[v]])
hell <- function(cn) decostand(t(cpm(calcNormFactors(DGEList(counts=as.matrix(cn)), method="TMM"))), "hellinger")
cnt <- list(); for (g in c("strict_final", "strict_final_EU"))
  cnt[[g]] <- read.delim(gzfile(sprintf("%s/step34e_counts_%s.tsv.gz", CD, g)), row.names=1, check.names=FALSE)[, meta$stem]
GN <- c(strict_final="all", strict_final_EU="EU")
dist_for <- function(g, m) { c0 <- cnt[[g]][, m$stem]; k <- round(0.8 * nrow(m)); a <- c0[rowSums(c0 > 0) >= k, ]
  list(d=vegdist(hell(a), "bray"), ngenes=nrow(a), k=k) }
res <- NULL
run_test <- function(set, g, m, form) {
  m <- droplevels(m); dd <- dist_for(g, m); d <- dd$d
  h <- how(nperm=999, plots=Plots(strata=m$animal, type="none"), within=Within(type="free"))
  set.seed(1); r <- adonis2(as.formula(paste("d", form)), data=m, permutations=h, by="terms")
  tot <- numPerms(nrow(m), h); np <- if (tot - 1 < getMinperm(h)) tot - 1 else getNperm(h)  # complete enumeration below minperm
  tt <- as.data.frame(r); mn <- table(factor(m$month, levels=c(1, 2, 3)))
  for (t in rownames(tt)) if (!t %in% c("Residual", "Total"))
    res <<- rbind(res, data.frame(set=set, genes=GN[[g]], n_samples=nrow(m), n_animals=nlevels(m$animal),
      months=paste(mn, collapse="/"), terms=paste0(t, " (", sub("^ *~ *", "", form), ")"), R2=round(tt[t, "R2"], 5),
      F=round(tt[t, "F"], 4), p=tt[t, "Pr(>F)"], n_perm=np, genes_used=dd$ngenes, min_samples=dd$k))
  r
}
# ---- reproduce
cat("== reproduction\n")
for (g in names(cnt)) run_test("reproduce_all44", g, meta, "~ month")
m23 <- meta[meta$month %in% c(2, 3), ]; run_test("reproduce_months2_3", "strict_final", m23, "~ month")
rp <- res; print(rp)
exp <- c(0.14112, 0.12574, 0.03186); got <- rp$R2
cat("expected", exp, "| got", got, "\n")
if (any(abs(round(got, 3) - round(exp, 3)) > 0)) { cat("STOP: reproduction differs at the third decimal\n"); quit(status=3) }
cat("reproduction OK\n")
# ---- A1, A2: shared-run set
ms <- meta[meta$run %in% SHARED, ]
cat("\nshared-run set:", nrow(ms), "samples,", length(unique(ms$animal)), "animals\n"); print(table(ms$month, ms$run))
stopifnot(nrow(ms) == 34, sum(ms$month == 1) == 4)
for (g in names(cnt)) { run_test("A1_shared_runs", g, ms, "~ month"); run_test("A2_shared_runs_run_first", g, ms, "~ run + month") }
# ---- B: the 4 animals whose month-1 metagenome is on a shared run
an <- unique(as.character(meta$animal[meta$month == 1 & meta$run %in% SHARED]))
mb <- meta[meta$animal %in% an, ]; em <- table(as.character(mb$egg_mass[mb$month == 1]))
cat("\nB animals:", an, "| month-1 egg masses:", paste(names(em), em, sep=":"), "| samples", nrow(mb), "\n")
if (length(an) != 4 || !identical(as.integer(em[c("1", "2")]), c(1L, 3L)) || nrow(mb) != 12) { cat("STOP: B set is not 4 animals (1 EM1, 3 EM2) x 3\n"); quit(status=4) }
for (g in names(cnt)) print(run_test("B_4_animals_shared_month1", g, mb, "~ month"))
write.table(res, O1, sep="\t", quote=FALSE, row.names=FALSE)
# ---- C: per-animal distances, all genes, all-44 matrix (same filter and transformation as the 44-sample test)
D <- as.matrix(dist_for("strict_final", meta)$d); C <- NULL
for (a in levels(meta$animal)) {
  s <- meta[meta$animal == a, ]; if (!1 %in% s$month) next
  st <- function(mo) { x <- s$stem[s$month == mo]; if (length(x) != 1) stop("animal ", a, " month ", mo, ": ", length(x)); x }
  s1 <- st(1); s2 <- st(2); s3 <- st(3); d12 <- D[s1, s2]; d13 <- D[s1, s3]; d23 <- D[s2, s3]
  C <- rbind(C, data.frame(animal=a, egg_mass=as.character(s$egg_mass[1]), month1_run=s$run[s$month == 1],
    shared=ifelse(s$run[s$month == 1] %in% SHARED, "yes", "no"), d12=round(d12, 4), d13=round(d13, 4), d23=round(d23, 4),
    excess=round((d12 + d13) / 2 - d23, 4)))
}
C <- C[order(C$shared, C$egg_mass, C$animal), ]; write.table(C, O2, sep="\t", quote=FALSE, row.names=FALSE)
cat("\nC: animals with a month-1 metagenome:", nrow(C), "\n"); print(C, row.names=FALSE)
for (s in c("yes", "no")) cat(sprintf("median excess, shared %s (n %d): %.4f\n", s, sum(C$shared == s), median(C$excess[C$shared == s])))
cat("\n"); print(res, row.names=FALSE); cat("wrote", O1, "and", O2, "\n")
