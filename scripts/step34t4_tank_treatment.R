# Step 34t4: tank (4.1, 4.2) by treatment code for the 32 WF24 EM4 froglets of the n = 88 set, exact test of independence
#   (fisher.test, two-sided). Stems from metadata/samples_88_treatment.tsv; tank, treatment, year and egg mass from
#   metadata/wf23_wf24_sequenced_metadata.tsv. UHM586.41010 is not in the n = 88 set; UHM590.41012 is egg mass 2.
# Output (refuses to overwrite): results/step34t/tank_treatment_fisher.txt
W <- "/bigdata/stajichlab/lshad003/wf_protein_space"
out <- file.path(W, "results/step34t/tank_treatment_fisher.txt")
if (file.exists(out)) stop("REFUSING: ", out, " exists")
n88 <- read.delim(file.path(W, "metadata/samples_88_treatment.tsv"), colClasses = "character")$stem
m <- read.delim(file.path(W, "metadata/wf23_wf24_sequenced_metadata.tsv"), colClasses = "character")
d <- m[m$stem %in% n88 & m$wf_year == "WF24" & m$egg_mass == "4", c("stem", "tank", "treatment")]
an <- sub("\\..*", "", d$stem)
if (nrow(d) != 32 || length(unique(an)) != 32) stop("ABORT: expected 32 froglets, got ", nrow(d))
if (any(an %in% c("UHM586", "UHM590"))) stop("ABORT: excluded froglet present")
if (!setequal(unique(d$tank), c("4.1", "4.2"))) stop("ABORT: tanks ", paste(unique(d$tank), collapse = ","))
t <- table(tank = d$tank, treatment = d$treatment)
f <- fisher.test(t)
sink(out)
cat("Step 34t4: WF24 egg mass 4, tank by treatment code, n = 88 set\n")
cat("froglets:", nrow(d), " distinct animals:", length(unique(an)), "\n\n")
print(addmargins(t))
print(f)
cat(sprintf("p-value (4 decimals): %.4f\n", f$p.value))
cat(R.version.string, "\n")
sink()
cat(readLines(out), sep = "\n")
