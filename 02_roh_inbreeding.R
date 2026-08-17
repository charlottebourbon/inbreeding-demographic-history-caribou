#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  02 - Runs of homozygosity (ROH) detection and inbreeding (FROH) analysis
#
#  Companion code for:
#    Bourbon et al. (2026) Inbreeding and Demographic History of Caribou
#    (Rangifer tarandus) in Western Canada Inferred From Genome-Wide SNP Data.
#    Evolutionary Applications. doi:10.1111/eva.70311
#
#  Pipeline (see Methods 2.3-2.4):
#    1. Detect ROHs with detectRUNS (consecutive method).
#    2. Per-individual ROH summaries: NROH, LROH, SROH.
#    3. Genome-wide inbreeding coefficient FROH.
#    4. FROH partitioned by ROH length class.
#    5. Group differences (subpopulation and metapopulation):
#       Kruskal-Wallis + epsilon-squared, Dunn's post-hoc, one-sample Wilcoxon.
#
#  ROH detection parameters (chosen after a sensitivity assessment; see
#  Supporting Information S1): consecutive method, minSNP = 20, maxOppRun = 2,
#  maxMissRun = 2, minLengthBps = 300,000. The maximum gap between consecutive
#  SNPs (1 Mb) is the detectRUNS default.
#
#  Requires: detectRUNS, dplyr, tidyr, FSA
#  Sources : froh_inbreeding_class.R  (custom FROH-by-class function)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(detectRUNS)
library(dplyr)
library(tidyr)
library(FSA)

source("froh_inbreeding_class.R")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
data_dir <- "data"      # input files
out_dir  <- "output"    # results written here
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Final autosomal dataset from script 01 (PED/MAP), no sex chromosome.
geno_stub <- file.path(data_dir, "caribou_759_autosomes")
ped_file  <- paste0(geno_stub, ".ped")
map_file  <- paste0(geno_stub, ".map")

# Maps each subpopulation to its metapopulation (structure) cluster.
# Expected columns: population, structure_cluster
cluster <- read.table(file.path(data_dir, "population_clusters.txt"),
                      header = TRUE, sep = "\t")

# Epsilon-squared effect size for a Kruskal-Wallis test.
#   H = KW chi-squared statistic, k = number of groups, n = total sample size
epsilon_squared <- function(H, k, n) (H - k + 1) / (n - k)


# ---------------------------------------------------------------------------
# 1. Detect ROHs (consecutive method)
# ---------------------------------------------------------------------------
runs <- consecutiveRUNS.run(
  genotypeFile = ped_file,
  mapFile      = map_file,
  minSNP       = 20,
  ROHet        = FALSE,       # homozygous runs
  minLengthBps = 300000,
  maxOppRun    = 2,
  maxMissRun   = 2
)

# Standardize subpopulation labels
runs$group <- gsub("_", "-", runs$group)

write.table(runs, file.path(out_dir, "roh_runs.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. Per-individual ROH summaries: NROH, LROH (Mb), SROH (Mb)
# ---------------------------------------------------------------------------
roh_per_ind <- runs %>%
  group_by(id, group) %>%
  dplyr::summarise(
    NROH = n(),
    SROH = sum(lengthBps) / 1e6,      # total ROH length, Mb
    LROH = mean(lengthBps) / 1e6,     # mean ROH length, Mb
    .groups = "drop"
  )

write.table(roh_per_ind, file.path(out_dir, "roh_summary_per_individual.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. Genome-wide inbreeding coefficient (FROH)
# ---------------------------------------------------------------------------
coeff_froh <- Froh_inbreeding(runs, mapFile = map_file)
coeff_froh$group <- gsub(" ", "-", coeff_froh$group)

write.table(coeff_froh, file.path(out_dir, "inbreeding_FROH.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Mean/SD/SE of FROH per subpopulation
froh_by_pop <- coeff_froh %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_FROH = mean(Froh_genome),
    sd_FROH   = sd(Froh_genome),
    se_FROH   = sd(Froh_genome) / sqrt(n()),
    .groups   = "drop"
  )
write.table(froh_by_pop, file.path(out_dir, "FROH_by_subpopulation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Attach metapopulation cluster, then summarize per metapopulation
coeff_froh <- merge(coeff_froh, cluster,
                    by.x = "group", by.y = "population")

froh_by_cluster <- coeff_froh %>%
  group_by(structure_cluster) %>%
  dplyr::summarise(
    mean_FROH = mean(Froh_genome),
    sd_FROH   = sd(Froh_genome),
    se_FROH   = sd(Froh_genome) / sqrt(n()),
    .groups   = "drop"
  )
write.table(froh_by_cluster, file.path(out_dir, "FROH_by_metapopulation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Overall mean FROH:", mean(coeff_froh$Froh_genome),
    "  SD:", sd(coeff_froh$Froh_genome),
    "  range:", paste(round(range(coeff_froh$Froh_genome), 3), collapse = "-"), "\n")


# ---------------------------------------------------------------------------
# 4. FROH by ROH length class
# ---------------------------------------------------------------------------
# Per-individual FROH in each length class (0.3-2, 2-4, 4-6, 6-8, >8 Mb)
coeff_froh_class <- froh_inbreeding_class(runs, map_file = map_file)

write.table(coeff_froh_class, file.path(out_dir, "FROH_by_length_class.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Mean FROH per class per subpopulation
mean_froh_class <- coeff_froh_class %>%
  group_by(group) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop")
write.table(mean_froh_class, file.path(out_dir, "FROH_by_length_class_mean.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Long format: one row per individual x class, for the class-contribution tests
class_long <- coeff_froh_class %>%
  select(group, id, starts_with("Froh_Class_")) %>%
  pivot_longer(cols = starts_with("Froh_Class_"),
               names_to = "length_class",
               names_prefix = "Froh_Class_",
               values_to = "FROH")


# ---------------------------------------------------------------------------
# 5a. Do length classes contribute differently to FROH?
# ---------------------------------------------------------------------------
kw_class <- kruskal.test(FROH ~ length_class, data = class_long)
eps_class <- epsilon_squared(kw_class$statistic,
                             length(unique(class_long$length_class)),
                             nrow(class_long))
print(kw_class)
cat("Effect size (epsilon^2) across length classes:", round(eps_class, 3), "\n")

if (kw_class$p.value < 0.05) {
  dunn_class <- dunnTest(FROH ~ length_class, data = class_long,
                         method = "bonferroni")
  print(dunn_class)
}


# ---------------------------------------------------------------------------
# 5b. FROH differences among subpopulations and metapopulations
# ---------------------------------------------------------------------------
# Among subpopulations
kw_pop  <- kruskal.test(Froh_genome ~ group, data = coeff_froh)
eps_pop <- epsilon_squared(kw_pop$statistic,
                           length(unique(coeff_froh$group)),
                           nrow(coeff_froh))
print(kw_pop)
cat("Effect size (epsilon^2), FROH among subpopulations:", round(eps_pop, 3), "\n")
if (kw_pop$p.value < 0.05) {
  dunn_pop <- dunnTest(Froh_genome ~ group, data = coeff_froh,
                       method = "bonferroni")$res
  write.table(dunn_pop, file.path(out_dir, "FROH_dunn_subpopulation.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Among metapopulations
kw_clu  <- kruskal.test(Froh_genome ~ structure_cluster, data = coeff_froh)
eps_clu <- epsilon_squared(kw_clu$statistic,
                           length(unique(coeff_froh$structure_cluster)),
                           nrow(coeff_froh))
print(kw_clu)
cat("Effect size (epsilon^2), FROH among metapopulations:", round(eps_clu, 3), "\n")
if (kw_clu$p.value < 0.05) {
  dunn_clu <- dunnTest(Froh_genome ~ structure_cluster, data = coeff_froh,
                       method = "bonferroni")$res
  write.table(dunn_clu, file.path(out_dir, "FROH_dunn_metapopulation.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


# ---------------------------------------------------------------------------
# 5c. One-sample Wilcoxon: each group vs the range-wide median FROH
# ---------------------------------------------------------------------------
# Tests whether each subpopulation's FROH is higher (or lower) than the overall
# median. Groups with n < 3 are skipped.
overall_median_froh <- median(coeff_froh$Froh_genome)

wilcox_by_pop <- coeff_froh %>%
  group_by(group) %>%
  filter(n() >= 3) %>%
  dplyr::summarise(
    p_greater = wilcox.test(Froh_genome, mu = overall_median_froh,
                            alternative = "greater")$p.value,
    p_less    = wilcox.test(Froh_genome, mu = overall_median_froh,
                            alternative = "less")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significant_greater = ifelse(p_greater < 0.05, "Yes", "No"),
    significant_less    = ifelse(p_less    < 0.05, "Yes", "No")
  )

write.table(wilcox_by_pop, file.path(out_dir, "FROH_wilcoxon_subpopulation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Same test at the metapopulation level
wilcox_by_cluster <- coeff_froh %>%
  group_by(structure_cluster) %>%
  filter(n() >= 3) %>%
  dplyr::summarise(
    p_greater = wilcox.test(Froh_genome, mu = overall_median_froh,
                            alternative = "greater")$p.value,
    p_less    = wilcox.test(Froh_genome, mu = overall_median_froh,
                            alternative = "less")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significant_greater = ifelse(p_greater < 0.05, "Yes", "No"),
    significant_less    = ifelse(p_less    < 0.05, "Yes", "No")
  )

write.table(wilcox_by_cluster, file.path(out_dir, "FROH_wilcoxon_metapopulation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# ---------------------------------------------------------------------------
# 5d. NROH and LROH differences among subpopulations
# ---------------------------------------------------------------------------
# Same Kruskal-Wallis + epsilon-squared framework applied to ROH number and
# mean ROH length, plus one-sample Wilcoxon vs the overall median.
for (metric in c("NROH", "LROH")) {

  kw <- kruskal.test(roh_per_ind[[metric]], as.factor(roh_per_ind$group))
  eps <- epsilon_squared(kw$statistic,
                         length(unique(roh_per_ind$group)),
                         nrow(roh_per_ind))
  cat("\n", metric, "among subpopulations: KW p =", signif(kw$p.value, 3),
      "  epsilon^2 =", round(eps, 3), "\n")

  overall_median <- median(roh_per_ind[[metric]])
  wilcox_metric <- roh_per_ind %>%
    group_by(group) %>%
    filter(n() >= 3) %>%
    dplyr::summarise(
      p_greater = wilcox.test(get(metric), mu = overall_median,
                              alternative = "greater")$p.value,
      p_less    = wilcox.test(get(metric), mu = overall_median,
                              alternative = "less")$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      significant_greater = ifelse(p_greater < 0.05, "Yes", "No"),
      significant_less    = ifelse(p_less    < 0.05, "Yes", "No")
    )

  write.table(wilcox_metric,
              file.path(out_dir, paste0(metric, "_wilcoxon_subpopulation.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\nDone. Results written to", out_dir, "\n")
