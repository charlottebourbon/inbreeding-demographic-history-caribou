#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  03 - Kinship (Fkin) estimation and comparison with FROH
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Pipeline (see Methods 2.4):
#    1. Estimate individual kinship with popkin and derive the kinship-based
#       inbreeding coefficient Fkin = 2 * phi_i - 1.
#    2. Merge Fkin with the ROH-based coefficient FROH (from script 02).
#    3. Compare the two metrics:
#         - Standardized Major Axis (SMA) regression, tested against the 1:1 line
#         - per-subpopulation OLS slopes (descriptive)
#         - Bland-Altman agreement (bias and limits of agreement)
#
#  Requires: popkin, BEDMatrix, dplyr, lmodel2
#  Inputs  : caribou_759_autosomes.{bed,bim,fam}  (final dataset, script 01)
#            inbreeding_FROH.txt                   (genome-wide FROH, script 02)
#            population_clusters.txt               (population -> metapopulation)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(popkin)
library(BEDMatrix)
library(dplyr)
library(lmodel2)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
data_dir <- "data"
out_dir  <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

geno_stub <- file.path(data_dir, "caribou_759_autosomes")
n_ind     <- 759
n_snp     <- 33346


# ---------------------------------------------------------------------------
# 1. Estimate kinship and the inbreeding coefficient Fkin
# ---------------------------------------------------------------------------
# Load genotypes (individuals x SNPs) and the accompanying .fam for labels.
X <- BEDMatrix(geno_stub, n = n_ind, p = n_snp, simple_names = TRUE)

fam <- read.table(paste0(geno_stub, ".fam"), header = FALSE)
colnames(fam)[1:2] <- c("group", "id")   # FID = subpopulation, IID = individual

# popkin expects markers x individuals, with subpopulation labels per column.
X <- t(X)
kinship <- popkin(X, subpops = fam$group)

# Individual inbreeding from kinship: Fkin = 2 * phi_i - 1 (popkin::inbr).
# inbr() returns values in individual order, i.e. matching `fam`.
Fkin <- data.frame(
  id    = fam$id,
  group = fam$group,
  Fkin  = inbr(kinship)
)

write.table(Fkin, file.path(out_dir, "inbreeding_Fkin.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Mean Fkin:", mean(Fkin$Fkin),
    "  SD:", sd(Fkin$Fkin), "\n")


# ---------------------------------------------------------------------------
# 2. Merge Fkin with FROH
# ---------------------------------------------------------------------------
# Both are keyed on id + group (from the same .fam), so no string parsing is
# needed. FROH is produced by script 02.
froh <- read.table(file.path(out_dir, "inbreeding_FROH.txt"),
                   header = TRUE, sep = "\t")
# script 02 writes the genome-wide coefficient as Froh_genome
froh <- froh[, c("id", "group", "Froh_genome")]

both <- merge(froh, Fkin, by = c("id", "group"))
cat("Individuals with both metrics:", nrow(both), "\n")

write.table(both, file.path(out_dir, "FROH_vs_Fkin.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---------------------------------------------------------------------------
# 3a. SMA regression of Fkin on FROH, tested against the 1:1 line
# ---------------------------------------------------------------------------
# Standardized Major Axis (SMA) accounts for error in both variables. A slope
# whose 95% CI excludes 1 indicates the two metrics diverge from 1:1.
model_ii <- lmodel2(Fkin ~ Froh_genome, data = both, nperm = 999)
print(model_ii)

sma <- model_ii$regression.results[model_ii$regression.results$Method == "SMA", ]
sma_slope     <- sma$Slope
sma_intercept <- sma$Intercept
sma_ci_lower  <- model_ii$confidence.intervals[
  model_ii$confidence.intervals$Method == "SMA", "2.5%-Slope"]
sma_ci_upper  <- model_ii$confidence.intervals[
  model_ii$confidence.intervals$Method == "SMA", "97.5%-Slope"]

cat(sprintf("SMA slope: %.3f (95%% CI: %.3f - %.3f)\n",
            sma_slope, sma_ci_lower, sma_ci_upper))
cat("R^2:", round(model_ii$rsquare, 3), "\n")
cat("Slope differs from 1:1:",
    ifelse(sma_ci_lower > 1 | sma_ci_upper < 1, "yes", "no"), "\n")


# ---------------------------------------------------------------------------
# 3b. Per-subpopulation OLS slopes (descriptive)
# ---------------------------------------------------------------------------
# SMA is unreliable at small n, so subpopulation-level slopes use OLS.
group_slopes <- both %>%
  group_by(group) %>%
  filter(n() >= 2) %>%
  dplyr::summarise(
    n     = n(),
    slope = coef(lm(Fkin ~ Froh_genome))[2],
    .groups = "drop"
  )

write.table(group_slopes, file.path(out_dir, "FROH_vs_Fkin_group_slopes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---------------------------------------------------------------------------
# 3c. Bland-Altman agreement between Fkin and FROH
# ---------------------------------------------------------------------------
# Plots the difference against the mean of the two measures (computed here as
# summary stats): bias = mean difference; limits of agreement = bias +/- 1.96 SD.
both$mean_value <- (both$Fkin + both$Froh_genome) / 2
both$difference <- both$Fkin - both$Froh_genome

bias      <- mean(both$difference)
sd_diff   <- sd(both$difference)
lower_loa <- bias - 1.96 * sd_diff
upper_loa <- bias + 1.96 * sd_diff

cat(sprintf("Bland-Altman bias (Fkin - FROH): %.4f +/- %.4f (SD)\n",
            bias, sd_diff))
cat(sprintf("Limits of agreement: %.4f to %.4f\n", lower_loa, upper_loa))

# Per-subpopulation bias, for the group-level agreement summary
group_bias <- both %>%
  group_by(group) %>%
  dplyr::summarise(
    n         = n(),
    mean_diff = mean(difference),
    sd_diff   = sd(difference),
    .groups   = "drop"
  ) %>%
  arrange(mean_diff)

write.table(group_bias, file.path(out_dir, "FROH_vs_Fkin_bland_altman_by_group.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Results written to", out_dir, "\n")
